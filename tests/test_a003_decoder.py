import inspect
import itertools
import math
import unittest

from a002_qualification import PchkMatrix
from exact_scoring import bits_to_bytes, score_exact_bytes
from ldpc_reference import sum_product_decode


def matrix_from_rows(rows):
    entries = tuple((r, c) for r, row in enumerate(rows)
                    for c, value in enumerate(row) if value)
    return PchkMatrix(len(rows), len(rows[0]), entries, 0, 0)


def exact_tree_posterior(rows, llrs):
    weights = []
    for bits in itertools.product((0, 1), repeat=len(llrs)):
        if any(sum(rows[r][c] * bits[c] for c in range(len(llrs))) % 2
               for r in range(len(rows))):
            continue
        log_weight = sum((1 - 2 * bit) * llr / 2
                         for bit, llr in zip(bits, llrs))
        weights.append((bits, math.exp(log_weight)))
    result = []
    for index in range(len(llrs)):
        p0 = sum(weight for bits, weight in weights if bits[index] == 0)
        p1 = sum(weight for bits, weight in weights if bits[index] == 1)
        result.append(math.log(p0 / p1))
    return tuple(result)


class A003DecoderTests(unittest.TestCase):
    def _invalid_input(self, *args, **kwargs):
        # An unimplemented decoder is not a valid input-rejection result.
        try:
            sum_product_decode(*args, **kwargs)
        except NotImplementedError:
            raise
        except ValueError:
            return
        self.fail("invalid decoder input was accepted")

    def test_one_check_has_analytic_sum_product_messages(self):
        matrix = matrix_from_rows(((1, 1, 1),))
        llrs = (0.6, -0.4, 0.8)
        result = sum_product_decode(matrix, llrs, 1, early_stop=False)
        expected = []
        for index, llr in enumerate(llrs):
            product = 1.0
            for other, other_llr in enumerate(llrs):
                if other != index:
                    product *= math.tanh(other_llr / 2)
            expected.append(llr + 2 * math.atanh(product))
        for actual, wanted in zip(result["posterior_llrs"], expected):
            self.assertAlmostEqual(actual, wanted, places=12)
        self.assertEqual(result["iterations"], 1)

    def test_tree_posterior_matches_exact_enumeration(self):
        rows = ((1, 1, 0), (0, 1, 1))
        llrs = (0.8, -1.1, 0.5)
        result = sum_product_decode(matrix_from_rows(rows), llrs, 2,
                                    early_stop=False)
        expected = exact_tree_posterior(rows, llrs)
        for actual, wanted in zip(result["posterior_llrs"], expected):
            self.assertAlmostEqual(actual, wanted, places=11)

    def test_nonzero_codeword_is_decoded_without_truth_argument(self):
        matrix = matrix_from_rows(((1, 1, 0), (0, 1, 1)))
        result = sum_product_decode(matrix, (-8.0, 8.0, -8.0), 5,
                                    early_stop=False)
        self.assertEqual(result["hard_bits"], (1, 1, 1))
        self.assertTrue(result["parity_satisfied"])

    def test_zero_llr_has_declared_bit_zero_tie(self):
        result = sum_product_decode(matrix_from_rows(((1, 1, 1),)),
                                    (0.0, 0.0, 0.0), 2, early_stop=False)
        self.assertEqual(result["hard_bits"], (0, 0, 0))
        self.assertEqual(result["iterations"], 2)
        self.assertTrue(all(value == 0.0 for value in result["posterior_llrs"]))

    def test_large_llrs_are_finite_and_clipped(self):
        result = sum_product_decode(matrix_from_rows(((1, 1, 1),)),
                                    (1000.0, -1000.0, 1000.0), 2,
                                    clip=10.0, early_stop=False)
        self.assertTrue(all(math.isfinite(value)
                            for value in result["posterior_llrs"]))
        self.assertTrue(all(abs(value) <= 10.0 for value in result["posterior_llrs"]))

    def test_early_stop_reports_initial_valid_word(self):
        result = sum_product_decode(matrix_from_rows(((1, 1),)),
                                    (5.0, 5.0), 7, early_stop=True)
        self.assertEqual(result["iterations"], 0)
        self.assertEqual(result["hard_bits"], (0, 0))
        self.assertEqual(result["syndrome_weight"], 0)
        self.assertEqual(result["syndrome_checks"], 1)

    def test_iteration_cap_is_obeyed_when_early_stop_disabled(self):
        result = sum_product_decode(matrix_from_rows(((1, 1, 1),)),
                                    (0.4, -0.2, 0.3), 3, early_stop=False)
        self.assertEqual(result["iterations"], 3)
        self.assertEqual(result["check_updates"], 9)
        self.assertEqual(result["syndrome_checks"], 4)

    def test_parity_valid_does_not_mean_exact_recovery(self):
        result = sum_product_decode(matrix_from_rows(((1,) * 8,)),
                                    (0.0,) * 8, 0, early_stop=True)
        self.assertTrue(result["parity_satisfied"])
        self.assertEqual(result["hard_bits"], (0,) * 8)
        candidate = bits_to_bytes(result["hard_bits"], bit_order="msb")
        self.assertFalse(score_exact_bytes(candidate, b"\xff")["match"])
        self.assertNotIn("truth", result)

    def test_truth_and_scorer_are_not_decoder_inputs(self):
        parameters = tuple(inspect.signature(sum_product_decode).parameters)
        self.assertEqual(parameters,
                         ("matrix", "llrs", "max_iterations", "clip", "early_stop"))
        sum_product_decode(matrix_from_rows(((1, 1),)), (1.0, 1.0), 1)

    def test_invalid_matrix_and_llr_inputs_are_rejected(self):
        matrix = matrix_from_rows(((1, 1),))
        self._invalid_input(None, (1.0, 1.0), 1)
        self._invalid_input(matrix, (1.0,), 1)
        self._invalid_input(matrix, (True, 1.0), 1)
        self._invalid_input(matrix, (1.0, 1.0), -1)
        self._invalid_input(matrix, (1.0, 1.0), 1, clip=0.0)
        self._invalid_input(matrix, (1.0, 1.0), 1, early_stop=1)
        for llrs in (None, [1.0, 1.0], (float("nan"), 0.0),
                     (float("inf"), 0.0), ("1", 0.0)):
            self._invalid_input(matrix, llrs, 1)
        for iterations in (True, 1.5, None):
            self._invalid_input(matrix, (1.0, 1.0), iterations)
        for clip in (True, 31.0, float("inf"), float("nan")):
            self._invalid_input(matrix, (1.0, 1.0), 1, clip=clip)

    def test_malformed_sparse_entries_are_rejected(self):
        duplicate = PchkMatrix(1, 2, ((0, 0), (0, 0), (0, 1)), 0, 0)
        self._invalid_input(duplicate, (1.0, 1.0), 1)
        for entries in ([], ((0,),), ((True, 0),), ((0, 2),), ((-1, 0),)):
            self._invalid_input(PchkMatrix(1, 2, entries, 0, 0), (1.0, 1.0), 1)
        for rows, cols in ((0, 2), (2, 2), (True, 2), (1, 64801)):
            self._invalid_input(PchkMatrix(rows, cols, (), 0, 0), (1.0, 1.0), 1)

    def test_inputs_are_not_mutated(self):
        matrix = matrix_from_rows(((1, 1, 0), (0, 1, 1)))
        entries = matrix.entries
        llrs = [1.0, -1.0, 1.0]
        try:
            sum_product_decode(matrix, tuple(llrs), 1)
        except NotImplementedError:
            raise
        self.assertEqual(matrix.entries, entries)
        self.assertEqual(llrs, [1.0, -1.0, 1.0])


if __name__ == "__main__":
    unittest.main()
