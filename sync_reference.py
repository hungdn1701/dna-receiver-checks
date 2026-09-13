"""Log-domain FBA for the explicit A004 declared-window I/D/T channel.

This is a source adaptation, not a reconstruction of an opaque upstream FBA.
The receiver sees an observation and public priors, never the true payload.
"""

import math

_BASES = "ATGC"
_NEG_INF = -math.inf


def _log(value):
    return math.log(value) if value > 0 else _NEG_INF


def _add(left, right):
    if left == _NEG_INF:
        return right
    if right == _NEG_INF:
        return left
    high, low = max(left, right), min(left, right)
    return high + math.log1p(math.exp(low - high))


def _sum(values):
    result = _NEG_INF
    for value in values:
        result = _add(result, value)
    return result


def _validate(read, priors, rates):
    if (not isinstance(read, str) or len(read) > 384
            or any(base not in _BASES for base in read)):
        raise ValueError("read must be an uppercase ATGC string of length <=384")
    if not isinstance(priors, tuple) or len(priors) > 256:
        raise ValueError("base_priors must be a tuple of at most 256 rows")
    for row in priors:
        if (not isinstance(row, tuple) or len(row) != 4
                or any(type(x) not in (int, float) or not 0 < x <= 1 for x in row)
                or abs(math.fsum(row) - 1) > 1e-12):
            raise ValueError("each prior must be four positive finite probabilities")
    if any(type(x) not in (int, float) or not 0 <= x < 1 for x in rates):
        raise ValueError("channel rates must be finite numbers in [0,1)")


def synchronize_read(read, base_priors, *, p_insert, p_delete, p_substitute):
    """Return likelihood, prior-removed log evidence, posteriors and term counts.

    Evidence rows use (A,T,G,C) order and maximum-zero normalization. Negative
    infinity is a valid zero likelihood. Impossible observations are rejected.
    """
    _validate(read, base_priors, (p_insert, p_delete, p_substitute))
    n, m = len(base_priors), len(read)
    prior_logs = tuple(tuple(math.log(x) for x in row) for row in base_priors)
    stop = math.log1p(-p_insert)
    insertion = _log(p_insert) - math.log(4)
    deletion = stop + _log(p_delete)
    transmit = stop + math.log1p(-p_delete)
    same = math.log1p(-p_substitute)
    different = _log(p_substitute) - math.log(3)
    # Take logs before division: positive subnormal rates must stay nonzero.
    emissions = tuple(tuple(same if base == observed else different
                            for base in _BASES) for observed in read)
    mixtures = [[_sum(p + e for p, e in zip(row, emission))
                 for emission in emissions] for row in prior_logs]

    forward = [[_NEG_INF] * (m + 1) for _ in range(n + 1)]
    forward[0][0] = 0.0
    for i in range(n + 1):
        for j in range(m + 1):
            value = forward[i][j]
            if j < m:
                forward[i][j + 1] = _add(forward[i][j + 1], value + insertion)
            if i < n:
                forward[i + 1][j] = _add(forward[i + 1][j], value + deletion)
            if i < n and j < m:
                forward[i + 1][j + 1] = _add(
                    forward[i + 1][j + 1], value + transmit + mixtures[i][j])
    log_z = forward[n][m] + stop
    if log_z == _NEG_INF:
        raise ValueError("observation has zero probability under the declared channel")

    backward = [[_NEG_INF] * (m + 1) for _ in range(n + 1)]
    backward[n][m] = stop
    for i in range(n, -1, -1):
        for j in range(m, -1, -1):
            if i == n and j == m:
                continue
            value = _NEG_INF
            if j < m:
                value = _add(value, insertion + backward[i][j + 1])
            if i < n:
                value = _add(value, deletion + backward[i + 1][j])
            if i < n and j < m:
                value = _add(value, transmit + mixtures[i][j] + backward[i + 1][j + 1])
            backward[i][j] = value

    evidence, posteriors = [], []
    for i in range(n):
        row = []
        for base in range(4):
            value = _NEG_INF
            for j in range(m + 1):
                value = _add(value, forward[i][j] + deletion + backward[i + 1][j])
                if j < m:
                    value = _add(value, forward[i][j] + transmit
                                 + emissions[j][base] + backward[i + 1][j + 1])
            row.append(value)
        peak = max(row)
        evidence.append(tuple(value - peak for value in row))
        posteriors.append(tuple(math.exp(p + value - log_z)
                                for p, value in zip(prior_logs[i], row)))
    transitions = (n + 1) * m + n * (m + 1) + n * m
    return {"log_likelihood": log_z, "log_evidence": tuple(evidence),
            "posteriors": tuple(posteriors),
            "work": {"lattice_cells": (n + 1) * (m + 1),
                     "forward_transition_terms": transitions,
                     "backward_transition_terms": transitions,
                     "extrinsic_terms": 4 * n * (2 * m + 1)}}
