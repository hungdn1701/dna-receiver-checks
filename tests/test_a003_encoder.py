import itertools
import unittest

from a002_qualification import PchkMatrix
from ldpc_reference import encode_systematic


def matrix_from_rows(rows):
    entries = tuple((r, c) for r, row in enumerate(rows)
                    for c, value in enumerate(row) if value)
    return PchkMatrix(len(rows), len(rows[0]), entries, 0, 0)


class A003EncoderTests(unittest.TestCase):
    def setUp(self):
        self.rows = ((1, 0, 1, 1), (1, 1, 0, 1))
        self.matrix = matrix_from_rows(self.rows)

    def test_encode_matches_four_explicit_words(self):
        expected = {(0, 0): (0, 0, 0, 0), (0, 1): (1, 0, 0, 1),
                    (1, 0): (1, 1, 1, 0), (1, 1): (0, 1, 1, 1)}
        for payload, word in expected.items():
            self.assertEqual(encode_systematic(self.matrix, payload), word)

    def test_independent_dense_parity_and_suffix(self):
        for payload in itertools.product((0, 1), repeat=2):
            word = encode_systematic(self.matrix, payload)
            self.assertIsInstance(word, tuple)
            self.assertEqual(word[2:], payload)
            for row in self.rows:
                self.assertEqual(sum(a * b for a, b in zip(row, word)) % 2, 0)

    def test_does_not_mutate_inputs(self):
        payload = (1, 0)
        original_entries = self.matrix.entries
        encode_systematic(self.matrix, payload)
        self.assertEqual(self.matrix.entries, original_entries)
        self.assertEqual(payload, (1, 0))

    def test_single_parity_row(self):
        matrix = matrix_from_rows(((1, 1, 1),))
        self.assertEqual(encode_systematic(matrix, (1, 1)), (0, 1, 1))
        self.assertEqual(encode_systematic(matrix, (1, 0)), (1, 1, 0))

    def test_missing_parity_diagonal_rejected(self):
        with self.assertRaises(ValueError):
            encode_systematic(matrix_from_rows(((0, 0, 1, 1), (1, 1, 0, 1))), (1, 0))

    def test_non_bidiagonal_parity_block_rejected(self):
        for rows in (((1, 1, 1, 0), (1, 1, 0, 1)),
                     ((1, 0, 1, 1), (0, 1, 0, 1))):
            with self.assertRaises(ValueError):
                encode_systematic(matrix_from_rows(rows), (1, 0))

    def test_duplicate_and_bad_sparse_entries_rejected(self):
        for extra in (((0, 0),), ((-1, 0),), ((0, 4),), ((2, 0),),
                      ((True, 0),), ((0.0, 1),), ((0,),), ("bad",)):
            matrix = PchkMatrix(2, 4, self.matrix.entries + extra, 0, 0)
            with self.assertRaises(ValueError):
                encode_systematic(matrix, (1, 0))

    def test_invalid_dimensions_and_matrix_types_rejected(self):
        bad = (None, PchkMatrix(0, 4, (), 0, 0),
               PchkMatrix(2, 2, (), 0, 0), PchkMatrix(True, 4, (), 0, 0),
               PchkMatrix(2, 64801, (), 0, 0), PchkMatrix(2, 4, [], 0, 0))
        for matrix in bad:
            with self.assertRaises(ValueError):
                encode_systematic(matrix, (1, 0))

    def test_wrong_payload_length_rejected(self):
        for payload in ((), (1,), (1, 0, 1)):
            with self.assertRaises(ValueError):
                encode_systematic(self.matrix, payload)

    def test_nonbinary_or_nontuple_payload_rejected(self):
        for payload in (None, [1, 0], (True, 0), (1.0, 0), (2, 0), ("1", 0)):
            with self.assertRaises(ValueError):
                encode_systematic(self.matrix, payload)

    def test_sparse_entry_order_does_not_change_encoding(self):
        reordered = PchkMatrix(2, 4, self.matrix.entries[::-1], 0, 0)
        self.assertEqual(encode_systematic(reordered, (1, 0)), (1, 1, 1, 0))


if __name__ == "__main__":
    unittest.main()
