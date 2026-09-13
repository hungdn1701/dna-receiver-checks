import unittest

from exact_scoring import bits_to_bytes, parse_bits, score_exact_bytes


class ExactScoringTests(unittest.TestCase):
    def test_ascii_whitespace_and_exact_count(self):
        self.assertEqual(parse_bits(b" 0\t1\r0\n1\v\f", 4), (0, 1, 0, 1))

    def test_short_and_extra_bits_are_rejected(self):
        for text in (b"00", b"0001", b"000 1"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_bits(text, 3)

    def test_nonbinary_suffix_is_not_ignored(self):
        for text in (b"000x", b"0N0", b"0\xa00", b"000\x00"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_bits(text, 3)

    def test_invalid_bit_counts_are_rejected(self):
        for count in (0, -1, True, 3.0, "3", None):
            with self.subTest(count=count), self.assertRaises(ValueError):
                parse_bits(b"000", count)

    def test_parser_requires_bytes(self):
        for text in ("000", bytearray(b"000"), None):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_bits(text, 3)

    def test_bit_order_uses_independent_hex_vector(self):
        bits = (1, 0, 0, 1, 0, 1, 1, 0)
        self.assertEqual(bits_to_bytes(bits, bit_order="msb"), b"\x96")
        self.assertEqual(bits_to_bytes(bits, bit_order="lsb"), b"\x69")

    def test_packing_preserves_zero_ff_and_empty(self):
        self.assertEqual(bits_to_bytes((0,) * 8 + (1,) * 8, bit_order="msb"),
                         b"\x00\xff")
        self.assertEqual(bits_to_bytes((), bit_order="msb"), b"")

    def test_packing_requires_explicit_order(self):
        with self.assertRaises(TypeError):
            bits_to_bytes((0,) * 8)

    def test_unknown_order_is_not_inferred(self):
        for order in ("auto", "big", "", None):
            with self.subTest(order=order), self.assertRaises(ValueError):
                bits_to_bytes((0,) * 8, bit_order=order)

    def test_no_implicit_padding_or_nonbinary_values(self):
        for bits in ((0,) * 7, (0,) * 9, (2,) + (0,) * 7,
                     (False,) + (0,) * 7, (0.0,) + (0,) * 7,
                     ("0",) * 8, [0] * 8, None):
            with self.subTest(bits=bits), self.assertRaises(ValueError):
                bits_to_bytes(bits, bit_order="msb")

    def test_full_byte_equality(self):
        for payload in (b"", b"\x00\xff", b"abc\n"):
            with self.subTest(payload=payload):
                self.assertEqual(score_exact_bytes(payload, payload),
                                 {"match": True, "actual_bytes": len(payload),
                                  "expected_bytes": len(payload),
                                  "byte_mismatches": 0})

    def test_same_prefix_different_length_never_matches(self):
        for candidate, truth in ((b"abcX", b"abc"), (b"ab", b"abc"),
                                 (b"", b"a"), (b"a", b"")):
            with self.subTest(candidate=candidate, truth=truth):
                result = score_exact_bytes(candidate, truth)
                self.assertFalse(result["match"])
                self.assertIsNone(result["byte_mismatches"])
                self.assertEqual(result["actual_bytes"], len(candidate))
                self.assertEqual(result["expected_bytes"], len(truth))

    def test_byte_mismatch_is_not_a_bit_ber(self):
        result = score_exact_bytes(b"\xffb", b"\x00b")
        self.assertFalse(result["match"])
        self.assertEqual(result["byte_mismatches"], 1)

    def test_scorer_does_not_coerce_other_types(self):
        for candidate, truth in (("abc", b"abc"), (b"abc", "abc"),
                                 (bytearray(b"x"), b"x"), (None, b"")):
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                score_exact_bytes(candidate, truth)


if __name__ == "__main__":
    unittest.main()
