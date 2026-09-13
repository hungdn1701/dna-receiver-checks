"""Bounded Neal sparse/pchk readers for the DS1-A002 asset gate.

The pchk envelope is distinct from the embedded sparse matrix body.
Reference: https://github.com/radfordneal/LDPC-codes/blob/master/rcode.c
The inspected source hash/notice is retained in format-clarification-01.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class PchkMatrix:
    rows: int
    cols: int
    entries: tuple[tuple[int, int], ...]
    consumed: int
    trailing: int


def _i32_le(blob: bytes, offset: int) -> int:
    end = offset + 4
    if end > len(blob):
        raise ValueError("truncated int32")
    return int.from_bytes(blob[offset:end], "little", signed=True)


def parse_sparse_matrix(blob: bytes) -> PchkMatrix:
    """Read Neal's sparse body, returning any unconsumed trailing byte count."""
    if not isinstance(blob, bytes) or len(blob) > 2 * 1024 * 1024:
        raise ValueError("sparse body must be bounded bytes")
    rows = _i32_le(blob, 0)
    cols = _i32_le(blob, 4)
    if not (0 < rows <= 64800 and 0 < cols <= 64800):
        raise ValueError("invalid pchk dimensions")
    offset = 8
    row = -1
    entries: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    while True:
        value = _i32_le(blob, offset)
        offset += 4
        if value == 0:
            return PchkMatrix(rows, cols, tuple(entries), offset, len(blob) - offset)
        if value < 0:
            row = -value - 1
            if row >= rows:
                raise ValueError("pchk row out of range")
            continue
        if row < 0:
            raise ValueError("pchk column before row")
        col = value - 1
        if col >= cols:
            raise ValueError("pchk column out of range")
        entry = (row, col)
        if entry in seen:
            raise ValueError("duplicate pchk entry")
        seen.add(entry)
        entries.append(entry)


def parse_pchk(blob: bytes) -> PchkMatrix:
    """Read one full pchk file; reject an unknown marker or trailing bytes."""
    if not isinstance(blob, bytes) or _i32_le(blob, 0) != (ord("P") << 8) + 0x80:
        raise ValueError("invalid pchk file marker")
    matrix = parse_sparse_matrix(blob[4:])
    if matrix.trailing:
        raise ValueError("trailing pchk bytes")
    return PchkMatrix(matrix.rows, matrix.cols, matrix.entries,
                      matrix.consumed + 4, 0)


def parse_int32_candidates(blob: bytes, count: int) -> dict[str, tuple[int, ...]]:
    if (not isinstance(blob, bytes) or type(count) is not int
            or not 0 < count <= 64800 or len(blob) != 4 * count):
        raise ValueError("int32 array has unexpected size")
    candidates: dict[str, tuple[int, ...]] = {}
    for endian in ("little", "big"):
        values = tuple(
            int.from_bytes(blob[offset:offset + 4], endian, signed=True)
            for offset in range(0, len(blob), 4)
        )
        if set(values) == set(range(count)):
            candidates[endian] = values
    return candidates


def check_syndrome(matrix: PchkMatrix, codeword: tuple[int, ...]) -> tuple[int, ...]:
    """Return Hx over GF(2), validating parser- or caller-constructed matrices."""
    if not isinstance(matrix, PchkMatrix):
        raise ValueError("matrix must be a PchkMatrix")
    if (type(matrix.rows) is not int or type(matrix.cols) is not int
            or not 0 < matrix.rows <= 64800 or not 0 < matrix.cols <= 64800):
        raise ValueError("matrix dimensions must be integers in [1, 64800]")
    if not isinstance(matrix.entries, tuple):
        raise ValueError("matrix entries must be a tuple")
    if (not isinstance(codeword, tuple) or len(codeword) != matrix.cols
            or any(type(bit) is not int or bit not in (0, 1) for bit in codeword)):
        raise ValueError("codeword does not match pchk columns")
    syndrome = [0] * matrix.rows
    seen = set()
    for entry in matrix.entries:
        if not isinstance(entry, tuple) or len(entry) != 2:
            raise ValueError("matrix entries must be pairs")
        row, col = entry
        if (type(row) is not int or type(col) is not int
                or not 0 <= row < matrix.rows or not 0 <= col < matrix.cols
                or entry in seen):
            raise ValueError("matrix entry is invalid or duplicated")
        seen.add(entry)
        syndrome[row] ^= codeword[col]
    return tuple(syndrome)
