"""Public priors, exact-ID evidence merging and factorized A004 soft demapping.

Observation IDs are processing provenance, not molecular-independence labels.
The demapper deliberately discards joint DNA-plane and alignment correlations.
"""

import itertools
import math

from sync_reference import _sum as _logsum

_NEG_INF = -math.inf
_QUARTETS = tuple(itertools.product((0, 1), repeat=4))


def _watermark_length(watermark):
    if (not isinstance(watermark, tuple) or not 0 < len(watermark) <= 81000
            or len(watermark) % 10
            or any(type(bit) is not int or bit not in (0, 1) for bit in watermark)):
        raise ValueError("watermark must be a bounded binary tuple of length 10*k")
    return len(watermark) // 2


def _validate_rows(rows, maximum):
    if not isinstance(rows, tuple) or len(rows) > maximum:
        raise ValueError("evidence must be a bounded tuple of rows")
    for row in rows:
        if (not isinstance(row, tuple) or len(row) != 4
                or any(type(x) not in (int, float) or x > 0 or x != x for x in row)
                or max(row) != 0):
            raise ValueError("each ATGC log-evidence row must have maximum zero")
        try:
            finite = all(x == _NEG_INF or math.isfinite(x) for x in row)
        except OverflowError as error:
            raise ValueError("log evidence exceeds floating-point range") from error
        if not finite:
            raise ValueError("invalid log-evidence value")


def _log_product(values):
    """Sum log factors without disguising finite overflow as an exact zero."""
    if _NEG_INF in values:
        return _NEG_INF
    try:
        total = math.fsum(values)
    except OverflowError as error:
        raise ValueError("finite log-product overflow") from error
    if not math.isfinite(total):
        raise ValueError("finite log-product overflow")
    return total


def watermark_base_priors(watermark):
    """Return factorized public (A,T,G,C) priors; never inspect payload bits."""
    half = _watermark_length(watermark)
    zeros = tuple(11 / 16 if bit == 0 else 5 / 16 for bit in watermark)
    return tuple((upper * lower, upper * (1 - lower),
                  (1 - upper) * lower, (1 - upper) * (1 - lower))
                 for upper, lower in zip(zeros[:half], zeros[half:]))


def merge_evidence(packets, reference_length):
    """Sum unique observation packets once; reject ID conflicts/contradictions."""
    if type(reference_length) is not int or not 1 <= reference_length <= 40500:
        raise ValueError("reference length must be an integer in [1,40500]")
    if not isinstance(packets, tuple) or len(packets) > 4096:
        raise ValueError("packets must be a tuple of at most 4096 entries")
    unique = {}
    rows_seen = 0
    # Finish whole-batch validation and ID conflict checks before accumulation.
    for packet in packets:
        if not isinstance(packet, tuple) or len(packet) != 3:
            raise ValueError("packet must be (observation_id, start, evidence)")
        identity, start, rows = packet
        if not isinstance(identity, str) or not 0 < len(identity) <= 256:
            raise ValueError("invalid observation ID")
        if type(start) is not int or not 0 <= start <= reference_length:
            raise ValueError("invalid declared packet start")
        _validate_rows(rows, 256)
        if start + len(rows) > reference_length:
            raise ValueError("packet extends outside the reference")
        if identity in unique and unique[identity] != packet:
            raise ValueError("conflicting reuse of an observation ID")
        unique[identity] = packet
        rows_seen += len(rows)

    combined = [[0.0] * 4 for _ in range(reference_length)]
    coverage = [0] * reference_length
    rows_used = 0
    for _, start, rows in unique.values():
        for offset, row in enumerate(rows):
            index = start + offset
            combined[index] = [_log_product((left, right))
                               for left, right in zip(combined[index], row)]
            coverage[index] += 1
        rows_used += len(rows)
    normalized = []
    for row in combined:
        peak = max(row)
        if peak == _NEG_INF:
            raise ValueError("mutually contradictory observations")
        normalized.append(tuple(value - peak for value in row))
    return {"log_evidence": tuple(normalized), "coverage": tuple(coverage),
            "observations_seen": len(packets), "unique_observations": len(unique),
            "duplicate_observations": len(packets) - len(unique),
            "base_rows_seen": rows_seen, "base_rows_used": rows_used}


def demap_codeword(log_evidence, permutation, watermark, *, clip=30.0):
    """Return clipped log(P0/P1) in LDPC order under the fixed factorization."""
    half = _watermark_length(watermark)
    n = half * 8 // 5
    _validate_rows(log_evidence, half)
    if len(log_evidence) != half:
        raise ValueError("evidence length does not match watermark")
    if (not isinstance(permutation, tuple) or len(permutation) != n
            or any(type(i) is not int or not 0 <= i < n for i in permutation)
            or len(set(permutation)) != n):
        raise ValueError("permutation must be a bijection over codeword indices")
    if type(clip) not in (int, float) or not 0 < clip <= 30:
        raise ValueError("clip must be finite and in (0,30]")

    upper, lower = [], []
    log_two = math.log(2)
    for a, t, g, c in log_evidence:
        upper.append((_logsum((a, t)) - log_two, _logsum((g, c)) - log_two))
        lower.append((_logsum((a, g)) - log_two, _logsum((t, c)) - log_two))
    binary_rows = tuple((row[1], row[0]) if mask else row
                        for row, mask in zip(upper + lower, watermark))
    patterns = []
    for quartet in _QUARTETS:
        flag = int(sum(quartet) >= 3)
        patterns.append((flag, *(bit ^ flag for bit in quartet)))
    llrs = [0.0] * n
    for offset in range(0, len(binary_rows), 5):
        group = binary_rows[offset:offset + 5]
        weights = tuple(_log_product(tuple(row[bit] for row, bit in zip(group, pattern)))
                        for pattern in patterns)
        if max(weights) == _NEG_INF:
            raise ValueError("no canonical quartet has nonzero likelihood")
        for bit_index in range(4):
            zero = _logsum(weight for quartet, weight in zip(_QUARTETS, weights)
                           if quartet[bit_index] == 0)
            one = _logsum(weight for quartet, weight in zip(_QUARTETS, weights)
                          if quartet[bit_index] == 1)
            llr = zero - one
            index = 4 * (offset // 5) + bit_index
            llrs[permutation[index]] = max(-clip, min(clip, llr))
    return {"llrs": tuple(llrs),
            "work": {"base_positions": half, "quartet_candidates": 16 * n // 4,
                     "bit_llrs": n}}
