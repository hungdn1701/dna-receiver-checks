"""Pure noiseless reference mapping for the DS1 A001 qualification."""

_BASES = {(0, 0): "A", (0, 1): "T", (1, 0): "G", (1, 1): "C"}
_PLANES = {base: pair for pair, base in _BASES.items()}


def _bits(value, name):
    if not isinstance(value, tuple) or any(isinstance(bit, bool) or not isinstance(bit, int) or bit not in (0, 1) for bit in value):
        raise ValueError(f"{name} must be a tuple of binary integers")


def _configuration(codeword, permutation, watermark):
    _bits(codeword, "codeword")
    _bits(watermark, "watermark")
    n = len(codeword)
    if n == 0 or n % 8:
        raise ValueError("codeword length must be a positive multiple of eight")
    if not isinstance(permutation, tuple) or len(permutation) != n or any(isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < n for index in permutation) or set(permutation) != set(range(n)):
        raise ValueError("permutation must be a bijection over codeword indices")
    if len(watermark) != 5 * n // 4:
        raise ValueError("watermark length mismatch")
    return n


def encode_reference(codeword: tuple[int, ...], permutation: tuple[int, ...], watermark: tuple[int, ...]) -> str:
    n = _configuration(codeword, permutation, watermark)
    permuted = tuple(codeword[index] for index in permutation)
    sparse = []
    for offset in range(0, n, 4):
        quartet = permuted[offset:offset + 4]
        flag = int(sum(quartet) >= 3)
        sparse.extend((flag, *(bit ^ flag for bit in quartet)))
    sparse = tuple(bit ^ mask for bit, mask in zip(sparse, watermark))
    half = len(sparse) // 2
    return "".join(_BASES[(sparse[i], sparse[half + i])] for i in range(half))


def decode_reference(reference: str, permutation: tuple[int, ...], watermark: tuple[int, ...]) -> tuple[int, ...]:
    if not isinstance(reference, str) or not reference or len(reference) % 5:
        raise ValueError("reference must be a non-empty ACGT string of length divisible by five")
    if any(base not in _PLANES for base in reference):
        raise ValueError("reference contains a non-ACGT base")
    n = len(reference) * 8 // 5
    if not isinstance(permutation, tuple) or len(permutation) != n or any(isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < n for index in permutation) or set(permutation) != set(range(n)):
        raise ValueError("permutation must be a bijection over decoded indices")
    _bits(watermark, "watermark")
    if len(watermark) != 5 * n // 4:
        raise ValueError("watermark length mismatch")
    half = len(reference)
    first_plane = tuple(_PLANES[base][0] for base in reference)
    second_plane = tuple(_PLANES[base][1] for base in reference)
    sparse = first_plane + second_plane
    unmasked = tuple(bit ^ mask for bit, mask in zip(sparse, watermark))
    permuted = []
    # Each encoded quartet occupies five bits, across the full sparse stream.
    for offset in range(0, len(unmasked), 5):
        flag, *encoded = unmasked[offset:offset + 5]
        quartet = tuple(bit ^ flag for bit in encoded)
        if flag != int(sum(quartet) >= 3):
            raise ValueError("non-canonical quartet")
        permuted.extend(quartet)
    codeword = [0] * n
    for index, source in enumerate(permutation):
        codeword[source] = permuted[index]
    return tuple(codeword)
