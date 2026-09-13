import unittest

from a002_qualification import (
    PchkMatrix,
    check_syndrome,
    parse_int32_candidates,
    parse_pchk as parse_pchk_file,
    parse_sparse_matrix as parse_pchk,
)


def i32(value):
    return int(value).to_bytes(4, "little", signed=True)


class A002QualificationTests(unittest.TestCase):
    def test_pchk_reads_rows_columns_and_one_based_sparse_entries(self):
        # H = [[1, 0, 1], [0, 1, 1]]
        blob = b"".join((i32(2), i32(3), i32(-1), i32(1), i32(3),
                         i32(-2), i32(2), i32(3), i32(0)))
        matrix = parse_pchk(blob)
        self.assertEqual((matrix.rows, matrix.cols), (2, 3))
        self.assertEqual(matrix.entries, ((0, 0), (0, 2), (1, 1), (1, 2)))
        self.assertEqual(matrix.consumed, len(blob))

    def test_pchk_rejects_duplicate_or_out_of_range_entries(self):
        duplicate = b"".join((i32(1), i32(2), i32(-1), i32(1), i32(1), i32(0)))
        with self.assertRaises(ValueError):
            parse_pchk(duplicate)
        bad_column = b"".join((i32(1), i32(2), i32(-1), i32(3), i32(0)))
        with self.assertRaises(ValueError):
            parse_pchk(bad_column)

    def test_pchk_reports_syndrome_over_gf2(self):
        blob = b"".join((i32(2), i32(3), i32(-1), i32(1), i32(3),
                         i32(-2), i32(2), i32(3), i32(0)))
        matrix = parse_pchk(blob)
        self.assertEqual(check_syndrome(matrix, (1, 1, 1)), (0, 0))
        self.assertEqual(check_syndrome(matrix, (1, 0, 0)), (1, 0))

    def test_int32_candidates_accept_only_exact_bijections(self):
        blob = b"".join(i32(value) for value in (2, 0, 1))
        candidates = parse_int32_candidates(blob, 3)
        self.assertEqual(candidates, {"little": (2, 0, 1)})

    def test_full_pchk_envelope_includes_standard_P_marker(self):
        # Independent bytes: marker 0x5080, 1 row, 2 columns, H=[1,1].
        blob = bytes.fromhex("80500000 01000000 02000000 ffffffff "
                             "01000000 02000000 00000000")
        matrix = parse_pchk_file(blob)
        self.assertEqual((matrix.rows, matrix.cols), (1, 2))
        self.assertEqual(matrix.entries, ((0, 0), (0, 1)))
        self.assertEqual(matrix.consumed, 28)
        self.assertEqual(check_syndrome(matrix, (1, 1)), (0,))
        self.assertEqual(check_syndrome(matrix, (1, 0)), (1,))

    def test_full_pchk_rejects_wrong_marker_trailer_and_truncation(self):
        blob = bytes.fromhex("80500000 01000000 02000000 ffffffff "
                             "01000000 02000000 00000000")
        for bad in (None, b"", blob[:3], blob[:-1], blob[4:],
                    b"XXXX" + blob[4:], blob + b"x"):
            with self.subTest(value=bad), self.assertRaises(ValueError):
                parse_pchk_file(bad)

    def test_sparse_rejects_missing_row_terminator_invalid_dims_and_indices(self):
        for values in ((1, 2, 1, 0), (1, 2, -2, 1, 0), (1, 2, -1, 1),
                       (0, 2, 0), (-1, 2, 0), (1, 0, 0),
                       (2147483647, 2, 0)):
            with self.subTest(values=values), self.assertRaises(ValueError):
                parse_pchk(b"".join(i32(x) for x in values))

    def test_sparse_body_reports_trailer_and_preserves_empty_row(self):
        blob = b"".join(i32(x) for x in (2, 3, -2, 2, 0))
        matrix = parse_pchk(blob + b"suffix")
        self.assertEqual(matrix.entries, ((1, 1),))
        self.assertEqual(matrix.trailing, 6)
        self.assertEqual(check_syndrome(matrix, (0, 1, 0)), (0, 1))

    def test_int32_ambiguity_invalid_counts_and_duplicate_values(self):
        self.assertEqual(parse_int32_candidates(i32(0), 1),
                         {"little": (0,), "big": (0,)})
        self.assertEqual(parse_int32_candidates(i32(0) + i32(0), 2), {})
        for count in (0, -1, True, 1.0, None, 64801):
            with self.subTest(count=count), self.assertRaises(ValueError):
                parse_int32_candidates(i32(0), count)
        with self.assertRaises(ValueError):
            parse_int32_candidates(i32(0) + b"x", 1)

    def test_syndrome_rejects_nonbit_and_nontuple_inputs(self):
        matrix = parse_pchk(b"".join(i32(x) for x in (1, 2, -1, 1, 2, 0)))
        for codeword in (None, [0, 1], (True, 0), (0.0, 1), (2, 0), (0,)):
            with self.subTest(codeword=codeword), self.assertRaises(ValueError):
                check_syndrome(matrix, codeword)

    def test_syndrome_rejects_invalid_caller_matrix_dimensions(self):
        dimensions = ((0, 2), (-1, 2), (True, 2), (1.0, 2), (None, 2),
                      (64801, 2), (1, 0), (1, -1), (1, True), (1, 2.0),
                      (1, None), (1, 64801))
        for rows, cols in dimensions:
            matrix = PchkMatrix(rows, cols, (), 0, 0)
            with self.subTest(rows=rows, cols=cols), self.assertRaises(ValueError):
                check_syndrome(matrix, (0, 0))

    def test_syndrome_rejects_negative_and_out_of_range_caller_indices(self):
        for entry in ((-1, 0), (0, -1), (2, 0), (0, 3)):
            matrix = PchkMatrix(2, 3, (entry,), 0, 0)
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                check_syndrome(matrix, (1, 0, 1))

    def test_syndrome_rejects_noninteger_caller_indices(self):
        for entry in ((True, 0), (0, False), (0.0, 0), (0, 1.0),
                      (None, 0), (0, "1")):
            matrix = PchkMatrix(2, 3, (entry,), 0, 0)
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                check_syndrome(matrix, (1, 0, 1))

    def test_syndrome_rejects_malformed_caller_entries(self):
        for entries in (None, [], [(0, 0)], (None,), ([0, 0],),
                        ((0,),), ((0, 0, 0),), ("00",)):
            matrix = PchkMatrix(2, 3, entries, 0, 0)
            with self.subTest(entries=entries), self.assertRaises(ValueError):
                check_syndrome(matrix, (1, 0, 1))

    def test_syndrome_rejects_duplicate_caller_entries(self):
        matrix = PchkMatrix(1, 2, ((0, 0), (0, 0)), 0, 0)
        with self.assertRaises(ValueError):
            check_syndrome(matrix, (1, 0))

    def test_syndrome_preserves_valid_rectangular_empty_and_boundary_cases(self):
        # Independent GF(2) expectations, including rows >= cols.
        cases = (
            (PchkMatrix(2, 3, ((1, 2), (0, 0), (0, 2)), 0, 0),
             (1, 0, 1), (0, 1)),
            (PchkMatrix(2, 2, ((0, 0), (1, 1)), 0, 0), (1, 0), (1, 0)),
            (PchkMatrix(2, 1, ((0, 0), (1, 0)), 0, 0), (1,), (1, 1)),
            (PchkMatrix(2, 3, (), 0, 0), (1, 0, 1), (0, 0)),
            (PchkMatrix(64800, 1, ((64799, 0),), 0, 0),
             (1,), (0,) * 64799 + (1,)),
            (PchkMatrix(1, 64800, ((0, 64799),), 0, 0),
             (0,) * 64799 + (1,), (1,)),
        )
        for matrix, codeword, expected in cases:
            with self.subTest(rows=matrix.rows, cols=matrix.cols):
                self.assertEqual(check_syndrome(matrix, codeword), expected)


if __name__ == "__main__":
    unittest.main()
