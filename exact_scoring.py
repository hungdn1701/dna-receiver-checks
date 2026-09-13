"""Independent exact-byte scoring boundary for DS1 qualification."""


def parse_bits(text: bytes, expected_bits: int) -> tuple[int, ...]:
    if not isinstance(text, bytes):
        raise ValueError("text must be bytes")
    if isinstance(expected_bits, bool) or not isinstance(expected_bits, int) or expected_bits <= 0:
        raise ValueError("expected_bits must be a positive integer")
    allowed_ws = b" \t\r\n\v\f"
    bits = []
    for value in text:
        if value in allowed_ws:
            continue
        if value not in (48, 49):
            raise ValueError("non-binary input")
        bits.append(value - 48)
    if len(bits) != expected_bits:
        raise ValueError("bit count mismatch")
    return tuple(bits)


def bits_to_bytes(bits: tuple[int, ...], *, bit_order: str) -> bytes:
    if not isinstance(bits, tuple) or any(
        isinstance(bit, bool) or not isinstance(bit, int) or bit not in (0, 1)
        for bit in bits
    ):
        raise ValueError("bits must be a tuple of binary integers")
    if bit_order not in ("msb", "lsb"):
        raise ValueError("bit_order must be msb or lsb")
    if len(bits) % 8:
        raise ValueError("bit length must be a multiple of eight")
    packed = bytearray()
    for offset in range(0, len(bits), 8):
        octet = bits[offset:offset + 8]
        value = 0
        for index, bit in enumerate(octet):
            shift = 7 - index if bit_order == "msb" else index
            value |= bit << shift
        packed.append(value)
    return bytes(packed)


def score_exact_bytes(candidate: bytes, truth: bytes) -> dict:
    if not isinstance(candidate, bytes) or not isinstance(truth, bytes):
        raise ValueError("candidate and truth must be bytes")
    equal_length = len(candidate) == len(truth)
    mismatches = sum(a != b for a, b in zip(candidate, truth)) if equal_length else None
    return {
        "match": equal_length and candidate == truth,
        "actual_bytes": len(candidate),
        "expected_bytes": len(truth),
        "byte_mismatches": mismatches,
    }
