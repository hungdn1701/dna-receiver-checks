"""Synthetic tests for evidence merging and watermark demapping."""

import inspect
import itertools
import math
import unittest

from a002_qualification import PchkMatrix
from imeta_mapping import encode_reference
from ldpc_reference import encode_systematic, sum_product_decode
from soft_evidence import demap_codeword, merge_evidence, watermark_base_priors
from sync_reference import synchronize_read

BASES = "ATGC"
PERM = (3, 5, 0, 7, 4, 2, 1, 6)
WATERMARK = (0, 1, 0, 0, 1, 1, 0, 1, 1, 0)
NEUTRAL = (0.0, 0.0, 0.0, 0.0)
ROW = (0.0, -1.0, -2.0, -3.0)


def hard_rows(reference, wrong=-math.inf):
    return tuple(tuple(0.0 if base == observed else wrong for base in BASES)
                 for observed in reference)


def all_word_llrs(rows):
    """Independent eight-bit enumeration, using only the qualified forward mapper."""
    joint = []
    for word in itertools.product((0, 1), repeat=8):
        reference = encode_reference(word, PERM, WATERMARK)
        weight = math.exp(sum(rows[i][BASES.index(base)]
                              for i, base in enumerate(reference)))
        joint.append((word, weight))
    return tuple(math.log(sum(w for bits, w in joint if bits[i] == 0)
                          / sum(w for bits, w in joint if bits[i] == 1))
                 for i in range(8))


class A004SoftTests(unittest.TestCase):
    def reject(self, function, *args, **kwargs):
        try:
            function(*args, **kwargs)
        except NotImplementedError:
            raise
        except ValueError:
            return
        self.fail("invalid soft-evidence input was accepted")

    def test_public_prior_matches_all_256_input_words(self):
        priors = watermark_base_priors(WATERMARK)
        counts = [[0] * 4 for _ in range(5)]
        for word in itertools.product((0, 1), repeat=8):
            for i, base in enumerate(encode_reference(word, PERM, WATERMARK)):
                counts[i][BASES.index(base)] += 1
        self.assertEqual(priors, tuple(tuple(value / 256 for value in row)
                                      for row in counts))

    def test_public_prior_api_cannot_accept_payload(self):
        watermark_base_priors(WATERMARK)
        self.assertEqual(tuple(inspect.signature(watermark_base_priors).parameters),
                         ("watermark",))

    def test_merge_overlaps_and_uncovered_positions(self):
        other = (-2.0, 0.0, -1.0, -3.0)
        result = merge_evidence((("one", 0, (ROW, ROW)), ("two", 1, (other,))), 3)
        self.assertEqual(result["log_evidence"],
                         (ROW, (-1.0, 0.0, -2.0, -5.0), NEUTRAL))
        self.assertEqual(result["coverage"], (1, 2, 0))
        self.assertEqual(result["base_rows_seen"], 3)
        self.assertEqual(result["base_rows_used"], 3)

    def test_identical_observation_id_is_counted_once(self):
        packet = ("one", 0, (ROW,))
        result = merge_evidence((packet, packet), 1)
        self.assertEqual(result["log_evidence"], (ROW,))
        self.assertEqual(result["coverage"], (1,))
        self.assertEqual(result["observations_seen"], 2)
        self.assertEqual(result["unique_observations"], 1)
        self.assertEqual(result["duplicate_observations"], 1)
        self.assertEqual(result["base_rows_seen"], 2)
        self.assertEqual(result["base_rows_used"], 1)

    def test_distinct_ids_with_identical_evidence_both_count(self):
        result = merge_evidence((("one", 0, (ROW,)), ("two", 0, (ROW,))), 1)
        self.assertEqual(result["log_evidence"], ((0.0, -2.0, -4.0, -6.0),))
        self.assertEqual(result["unique_observations"], 2)
        self.assertEqual(result["duplicate_observations"], 0)
        self.assertEqual(result["coverage"], (2,))

    def test_conflicting_reuse_of_an_id_is_rejected(self):
        one = ("one", 0, (ROW,))
        self.reject(merge_evidence, (one, ("one", 1, (ROW,))), 2)
        self.reject(merge_evidence, (one, ("one", 0, (NEUTRAL,))), 2)

    def test_contradictory_exact_evidence_is_not_fabricated_as_erasure(self):
        a, t = hard_rows("AT")
        self.reject(merge_evidence, (("a", 0, (a,)), ("t", 0, (t,))), 1)

    def test_erasure_demap_is_exactly_zero_llr(self):
        result = demap_codeword((NEUTRAL,) * 5, PERM, WATERMARK)
        self.assertEqual(result["llrs"], (0.0,) * 8)
        self.assertEqual(result["work"],
                         {"base_positions": 5, "quartet_candidates": 32, "bit_llrs": 8})

    def test_noiseless_inverse_for_all_256_words(self):
        for word in itertools.product((0, 1), repeat=8):
            rows = hard_rows(encode_reference(word, PERM, WATERMARK))
            result = demap_codeword(rows, PERM, WATERMARK)
            self.assertEqual(tuple(int(llr < 0) for llr in result["llrs"]), word)
            self.assertTrue(all(abs(llr) == 30 for llr in result["llrs"]))

    def test_factorizable_base_evidence_matches_direct_word_oracle(self):
        rows = []
        for i in range(5):
            # Product of two independent binary factors, so exact-word and
            # declared factorized demapping have the same marginal posterior.
            upper, lower = (2 + i, 7 - i), (3 + i, 5)
            weights = (upper[0] * lower[0], upper[0] * lower[1],
                       upper[1] * lower[0], upper[1] * lower[1])
            rows.append(tuple(math.log(x / max(weights)) for x in weights))
        rows = tuple(rows)
        result = demap_codeword(rows, PERM, WATERMARK)
        for actual, expected in zip(result["llrs"], all_word_llrs(rows)):
            self.assertAlmostEqual(actual, expected, places=11)

    def test_joint_base_correlations_are_explicitly_discarded(self):
        first = tuple(math.log(x / 4) for x in (4, 2, 1, 3))
        second = tuple(math.log(x / 3) for x in (3, 3, 2, 2))
        a = demap_codeword((first,) * 5, PERM, WATERMARK)
        b = demap_codeword((second,) * 5, PERM, WATERMARK)
        for left, right in zip(a["llrs"], b["llrs"]):
            self.assertAlmostEqual(left, right, places=11)

    def test_final_llr_clipping_preserves_sign_for_finite_evidence(self):
        word = (1, 0, 1, 1, 0, 0, 1, 0)
        rows = hard_rows(encode_reference(word, PERM, WATERMARK), wrong=-100.0)
        result = demap_codeword(rows, PERM, WATERMARK, clip=2.5)
        self.assertEqual(result["llrs"], tuple(-2.5 if bit else 2.5 for bit in word))

    def test_noncanonical_exact_quartet_is_rejected(self):
        # Upper sparse group 11100 has weight three, so it is not canonical.
        self.reject(demap_codeword, hard_rows("GGGAA"), tuple(range(8)), (0,) * 10)

    def test_watermark_values_and_geometry_are_validated(self):
        for bad in (None, list(WATERMARK), (), (0,) * 9, (0,) * 81010,
                    (True,) + WATERMARK[1:], (2,) + WATERMARK[1:]):
            self.reject(watermark_base_priors, bad)

    def test_message_rows_and_numerical_overflow_are_validated(self):
        for bad in ([0, 0, 0, 0], (0, 0, 0), (0, 0, True, 0),
                    (0, math.nan, 0, 0), (0, math.inf, 0, 0),
                    (0, 1, 0, 0), (-1, -1, -1, -1), (-math.inf,) * 4):
            self.reject(merge_evidence, (("one", 0, (bad,)),), 1)
        huge = (0.0, -1e308, -1e308, -1e308)
        self.reject(merge_evidence, (("one", 0, (huge,)), ("two", 0, (huge,))), 1)

    def test_packet_shape_bounds_and_id_are_validated(self):
        packet = ("one", 0, (ROW,))
        for packets, length in ((None, 1), ([packet], 1), ((packet,) * 4097, 1),
                                ((("", 0, (ROW,)),), 1), ((("x" * 257, 0, (ROW,)),), 1),
                                (((2, 0, (ROW,)),), 1), ((("one", True, (ROW,)),), 1),
                                ((("one", -1, (ROW,)),), 1), ((("one", 1, (ROW,)),), 1),
                                ((("one", 0, (ROW,) * 257),), 400),
                                ((("one", 0, [ROW]),), 1), ((("one", 0),), 1),
                                ((packet,), True), ((packet,), 0), ((packet,), 40501)):
            self.reject(merge_evidence, packets, length)

    def test_demap_geometry_permutation_and_clip_are_validated(self):
        self.reject(demap_codeword, (NEUTRAL,) * 4, PERM, WATERMARK)
        for bad in (None, list(PERM), PERM[:-1], (0,) * 8,
                    (True,) + PERM[1:], (8,) + PERM[1:]):
            self.reject(demap_codeword, (NEUTRAL,) * 5, bad, WATERMARK)
        for bad in (0, -1, 31, True, math.nan, math.inf, "3"):
            self.reject(demap_codeword, (NEUTRAL,) * 5, PERM, WATERMARK, clip=bad)

    def test_pure_interfaces_and_immutable_inputs(self):
        packets = (("one", 0, (ROW,)),)
        result = merge_evidence(packets, 1)
        self.assertEqual(packets, (("one", 0, (ROW,)),))
        self.assertEqual(result, merge_evidence(packets, 1))
        self.assertIsInstance(result["log_evidence"], tuple)
        self.assertIsInstance(result["coverage"], tuple)
        self.assertEqual(tuple(inspect.signature(merge_evidence).parameters),
                         ("packets", "reference_length"))
        self.assertEqual(tuple(inspect.signature(demap_codeword).parameters),
                         ("log_evidence", "permutation", "watermark", "clip"))

    def test_tiny_fba_to_ldpc_path_uses_public_priors_and_dedup(self):
        priors = watermark_base_priors(WATERMARK)
        entries = tuple((r, c) for r in range(4) for c in range(8)
                        if c == r or (r > 0 and c == r - 1) or c == 4 + r)
        matrix = PchkMatrix(4, 8, entries, 0, 0)
        word = encode_systematic(matrix, (1, 0, 1, 1))
        observed = encode_reference(word, PERM, WATERMARK)
        sync = synchronize_read(observed, priors, p_insert=0, p_delete=0,
                                p_substitute=0.000003)
        packet = ("synthetic-one", 0, sync["log_evidence"])
        merged = merge_evidence((packet, packet), 5)
        self.assertEqual(merged["unique_observations"], 1)
        soft = demap_codeword(merged["log_evidence"], PERM, WATERMARK)
        decoded = sum_product_decode(matrix, soft["llrs"], 4)
        self.assertEqual(decoded["hard_bits"], word)
        self.assertTrue(decoded["parity_satisfied"])

    def test_empty_batch_reports_no_observation_or_work(self):
        result = merge_evidence((), 3)
        self.assertEqual(result, {"log_evidence": (NEUTRAL,) * 3,
                                  "coverage": (0, 0, 0),
                                  "observations_seen": 0, "unique_observations": 0,
                                  "duplicate_observations": 0,
                                  "base_rows_seen": 0, "base_rows_used": 0})
