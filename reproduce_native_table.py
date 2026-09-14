"""Reproduce the five bit-order rows for one pinned external configuration.

Offline: inputs must already be present under their upstream filenames.
No upstream program, read decoder or original-file reconstruction is run.
See NATIVE_REPRODUCTION.md for access instructions and scope.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

from a002_qualification import PchkMatrix, check_syndrome, parse_int32_candidates, parse_pchk
from exact_scoring import parse_bits

UPSTREAM_COMMIT = "ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b"
# These expected pins and row values are assertions, not computed outputs.
ASSET_SPECS = {
    "SequenceLengthALL_FILE001R0667": (324000, "02a34a578244eeea001fb9d82c47848bd09c65f304bce1282e9b25b188c7af4f"),
    "dvb_s2_r2_3.pchk": (950412, "25b638c5e10ba2f436fa6b8d4493b64c91e1fdb28d2b2c7fa845a7f9d9e5b9b1"),
    "encoded_bit.txt": (64800, "bc42f06c5b6ff6a6a569bbcb2575365f9a95003b672ea117413ba5973c2cdfac"),
    "permutation64800": (259200, "42f751c6d95c782ccaa3171aacf354435f77bc4186e433f52e89827995416742"),
    "source_bit.txt": (43200, "0e8a9103ee89bb1689e0d54f263ee562c19be332115d7bdb7fa943b1f753853d"),
}
EXPECTED_ROWS = {
    "stored": {"syndrome_weight": 10988, "suffix_source_mismatches": 21458},
    "undo_public_permutation": {"syndrome_weight": 0, "suffix_source_mismatches": 0},
    "apply_public_permutation": {"syndrome_weight": 10658, "suffix_source_mismatches": 21420},
    "rotate_21600": {"syndrome_weight": 10779, "suffix_source_mismatches": 21666},
    "rotate_43200": {"syndrome_weight": 10805, "suffix_source_mismatches": 21521},
}


def read_pinned_file(path: Path, expected_size: int, expected_sha256: str) -> bytes:
    """Read bounded regular input bytes; never normalize or repair a mismatch."""
    if path.is_symlink() or not path.is_file():
        raise ValueError("input must be a regular, non-symlink file: " + path.name)
    with path.open("rb") as stream:
        blob = stream.read(expected_size + 1)
    if len(blob) != expected_size:
        raise ValueError("input size mismatch: " + path.name)
    if hashlib.sha256(blob).hexdigest() != expected_sha256:
        raise ValueError("input SHA-256 mismatch: " + path.name)
    return blob


def ordering_rows(matrix: PchkMatrix, stored: tuple[int, ...],
                  source: tuple[int, ...], permutation: tuple[int, ...]) -> dict:
    """Compute five fixed candidates; known source is used only for scoring."""
    if not isinstance(matrix, PchkMatrix):
        raise ValueError("expected PchkMatrix")
    check_syndrome(matrix, stored)  # Validate caller-built matrix and word.
    n, parity_length = matrix.cols, matrix.rows
    source_length = n - parity_length
    if (source_length <= 0 or source_length == parity_length
            or not isinstance(source, tuple) or len(source) != source_length
            or any(type(bit) is not int or bit not in (0, 1) for bit in source)):
        raise ValueError("source suffix does not match geometry")
    if (not isinstance(permutation, tuple) or len(permutation) != n
            or any(type(index) is not int for index in permutation)
            or set(permutation) != set(range(n))):
        raise ValueError("expected a complete integer permutation")
    inverse = [0] * n
    for index, destination in enumerate(permutation):
        inverse[destination] = stored[index]
    candidates = {
        "stored": stored,
        "undo_public_permutation": tuple(inverse),
        "apply_public_permutation": tuple(stored[index] for index in permutation),
        "rotate_" + str(parity_length): stored[parity_length:] + stored[:parity_length],
        "rotate_" + str(source_length): stored[source_length:] + stored[:source_length],
    }
    return {
        label: {
            "syndrome_weight": sum(check_syndrome(matrix, word)),
            "suffix_source_mismatches": sum(
                bit != truth for bit, truth in zip(word[-source_length:], source)),
        }
        for label, word in candidates.items()
    }


def verify_expected_rows(rows: dict) -> None:
    """Fail on any changed/missing candidate; do not replace calculated values."""
    if rows != EXPECTED_ROWS:
        raise ValueError("computed ordering rows differ from the published table")


def reproduce(asset_dir: Path) -> dict:
    blobs = {
        name: read_pinned_file(asset_dir / name, size, digest)
        for name, (size, digest) in ASSET_SPECS.items()
    }
    matrix = parse_pchk(blobs["dvb_s2_r2_3.pchk"])
    if (matrix.rows, matrix.cols, len(matrix.entries)) != (21600, 64800, 215999):
        raise ValueError("unexpected matrix geometry")
    stored = parse_bits(blobs["encoded_bit.txt"], 64800)
    source = parse_bits(blobs["source_bit.txt"], 43200)
    candidates = parse_int32_candidates(blobs["permutation64800"], 64800)
    if len(candidates) != 1 or "little" not in candidates:
        raise ValueError("unexpected or ambiguous permutation byte order")
    rows = ordering_rows(matrix, stored, source, candidates["little"])
    verify_expected_rows(rows)
    return {
        "upstream_commit": UPSTREAM_COMMIT,
        "inputs": [
            {"name": name, "bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest()}
            for name, blob in blobs.items()
        ],
        "matrix": {"rows": matrix.rows, "cols": matrix.cols, "entries": len(matrix.entries)},
        "ordering_checks": rows,
        "expected_rows_match": True,
        "scope": "one supplied tuple; five conventions; no sequencing or file recovery",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets", type=Path, required=True,
                        help="directory containing the five unmodified configuration files")
    args = parser.parse_args()
    try:
        result = reproduce(args.assets)
    except (OSError, ValueError) as error:
        print("Native reproduction failed: " + str(error), file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

