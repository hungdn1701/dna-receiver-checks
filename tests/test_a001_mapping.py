import unittest

from imeta_mapping import decode_reference, encode_reference


# Explicit truth table, not calls to the implementation under test.
QUARTETS = ("00000", "00001", "00010", "00011", "00100", "00101", "00110", "11000",
            "01000", "01001", "01010", "10100", "01100", "10010", "10001", "10000")


def reference_from_planes(sparse):
    half = len(sparse) // 2
    return "".join("ATGC"[2 * int(a) + int(b)] for a, b in zip(sparse[:half], sparse[half:]))


class MappingTests(unittest.TestCase):
    def test_forward_matches_all_256_quartet_pairs(self):
        for value in range(256):
            bits = tuple(map(int, f"{value:08b}"))
            expected = reference_from_planes(QUARTETS[value >> 4] + QUARTETS[value & 15])
            with self.subTest(value=value):
                self.assertEqual(encode_reference(bits, tuple(range(8)), (0,) * 10), expected)

    def test_inverse_matches_all_256_quartet_pairs(self):
        for value in range(256):
            reference = reference_from_planes(QUARTETS[value >> 4] + QUARTETS[value & 15])
            with self.subTest(value=value):
                self.assertEqual(decode_reference(reference, tuple(range(8)), (0,) * 10),
                                 tuple(map(int, f"{value:08b}")))

    def test_hand_vector_uses_both_planes_and_explicit_permutation_direction(self):
        bits = (0, 1, 1, 1, 1, 0, 1, 0)
        permutation = (2, 0, 1, 4, 5, 7, 3, 6)
        watermark = (0, 1, 1, 0, 0, 1, 1, 0, 1, 0)
        # Permuted B3 -> 10100 00011; XOR -> 11000 11001 -> CCAAT.
        self.assertEqual(encode_reference(bits, permutation, watermark), "CCAAT")
        self.assertEqual(decode_reference("CCAAT", permutation, watermark), bits)

    def test_all_quartets_round_trip_with_identity(self):
        bits = tuple((value >> (3 - index)) & 1 for value in range(16) for index in range(4))
        permutation = tuple(range(64))
        reference = encode_reference(bits, permutation, (0,) * 80)
        self.assertEqual(decode_reference(reference, permutation, (0,) * 80), bits)

    def test_nonzero_watermark_and_noninvolution_permutation(self):
        n = 32
        bits = tuple((index * 3 + 1) % 2 for index in range(n))
        permutation = tuple(list(range(1, n)) + [0])
        watermark = tuple(index % 2 for index in range(5 * n // 4))
        reference = encode_reference(bits, permutation, watermark)
        self.assertEqual(decode_reference(reference, permutation, watermark), bits)

    def test_native_size_geometry(self):
        n = 64800
        bits = tuple(index % 2 for index in range(n))
        reference = encode_reference(bits, tuple(range(n)), (0,) * (5 * n // 4))
        self.assertEqual(len(reference), 40500)
        decoded = decode_reference(reference, tuple(range(n)), (0,) * (5 * n // 4))
        # Keep the equality check; avoid building a huge tuple diff on failure.
        self.assertTrue(decoded == bits, "native-size reference round trip differs")

    def test_native_size_mixed_bits_nonzero_watermark_and_permutation(self):
        n = 64800
        bits = tuple((index // 4 >> (3 - index % 4)) & 1 for index in range(n))
        permutation = tuple((37 * index + 11) % n for index in range(n))
        watermark = tuple((index * index + index // 7) % 2 for index in range(5 * n // 4))
        reference = encode_reference(bits, permutation, watermark)
        self.assertEqual(len(reference), 40500)
        self.assertTrue(decode_reference(reference, permutation, watermark) == bits,
                        "mixed native-size reference round trip differs")

    def test_invalid_configurations_and_bases_rejected(self):
        with self.assertRaises(ValueError):
            encode_reference((0,) * 4, (0, 1, 2, 3), (0,) * 5)
        with self.assertRaises(ValueError):
            encode_reference((0,) * 8, (0, 0, 2, 3, 4, 5, 6, 7), (0,) * 10)
        with self.assertRaises(ValueError):
            decode_reference("AAAAA", tuple(range(8)), (0,) * 9)
        with self.assertRaises(ValueError):
            decode_reference("AAAAX", tuple(range(8)), (0,) * 10)

    def test_noncanonical_quartet_is_rejected(self):
        # 01111 has flag zero but a weight-four quartet: not canonical.
        with self.assertRaises(ValueError):
            decode_reference("AGGGG", tuple(range(8)), (0,) * 10)

    def test_all_invalid_five_bit_words_in_either_half(self):
        invalid = sorted({f"{value:05b}" for value in range(32)} - set(QUARTETS))
        self.assertEqual(len(invalid), 16)
        for word in invalid:
            for sparse in (word + "00000", "00000" + word):
                with self.subTest(sparse=sparse), self.assertRaises(ValueError):
                    decode_reference(reference_from_planes(sparse), tuple(range(8)), (0,) * 10)

    def test_flag_one_zero_tail_is_valid_all_ones_quartet(self):
        self.assertEqual(decode_reference("GAAAA", tuple(range(8)), (0,) * 10),
                         (1, 1, 1, 1, 0, 0, 0, 0))

    def test_invalid_bit_tuples(self):
        for bad in (None, [0] * 8, (True,) * 8, (0.0,) * 8, (2,) * 8, (), (0,) * 7):
            with self.subTest(codeword=bad), self.assertRaises(ValueError):
                encode_reference(bad, tuple(range(8)), (0,) * 10)
        for bad in (None, [0] * 10, (False,) * 10, (0.0,) * 10, (2,) * 10, (0,) * 11):
            with self.subTest(watermark=bad):
                with self.assertRaises(ValueError):
                    encode_reference((0,) * 8, tuple(range(8)), bad)
                with self.assertRaises(ValueError):
                    decode_reference("AAAAA", tuple(range(8)), bad)

    def test_invalid_permutations_in_both_directions(self):
        for bad in (None, list(range(8)), (0,) * 8, tuple(range(7)),
                    (-1,) + tuple(range(1, 8)), tuple(range(1, 9)),
                    (False,) + tuple(range(1, 8)), (0.0,) + tuple(range(1, 8))):
            with self.subTest(permutation=bad):
                with self.assertRaises(ValueError):
                    encode_reference((0,) * 8, bad, (0,) * 10)
                with self.assertRaises(ValueError):
                    decode_reference("AAAAA", bad, (0,) * 10)

    def test_invalid_reference_type_length_and_alphabet(self):
        for reference in (None, b"AAAAA", "", "AAAA", "AAAAAA", "aaaaa", "AAAAN", "AAAA\n"):
            with self.subTest(reference=reference), self.assertRaises(ValueError):
                decode_reference(reference, tuple(range(8)), (0,) * 10)


if __name__ == "__main__":
    unittest.main()
