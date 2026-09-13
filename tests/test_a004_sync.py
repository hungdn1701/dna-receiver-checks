"""Finite-oracle tests for declared-window synchronization."""

from fractions import Fraction
import inspect
import itertools
import math
import unittest

from sync_reference import synchronize_read

BASES = "ATGC"
UNIFORM = (0.25, 0.25, 0.25, 0.25)
RATES = dict(p_insert=0.2, p_delete=0.1, p_substitute=0.15)
PRIORS = ((0.1, 0.2, 0.3, 0.4), (0.4, 0.1, 0.2, 0.3))


def enumerate_joint(read, priors, rates):
    """Enumerate actual source words and edit paths with exact rational weights.

    No forward/backward arrays, log sums, caching or production helpers.
    Intended only for at most two source bases and three observed bases.
    """
    rows = tuple(tuple(Fraction(str(x)) for x in row) for row in priors)
    pi, pd, ps = (Fraction(str(rates[name])) for name in
                  ("p_insert", "p_delete", "p_substitute"))
    total = Fraction(0)
    joint = [[Fraction(0) for _ in BASES] for _ in rows]
    for word in itertools.product(range(4), repeat=len(rows)):
        path_weights = []

        def walk(i, j, weight):
            if i == len(word) and j == len(read):
                path_weights.append(weight * (1 - pi))
                return
            if j < len(read):
                walk(i, j + 1, weight * pi / 4)
            if i < len(word):
                walk(i + 1, j, weight * (1 - pi) * pd)
            if i < len(word) and j < len(read):
                emit = 1 - ps if BASES[word[i]] == read[j] else ps / 3
                walk(i + 1, j + 1, weight * (1 - pi) * (1 - pd) * emit)

        walk(0, 0, Fraction(1))
        weight = sum(path_weights) * math.prod(rows[i][b] for i, b in enumerate(word))
        total += weight
        for i, base in enumerate(word):
            joint[i][base] += weight
    return total, tuple(tuple(row) for row in joint)


class A004SyncTests(unittest.TestCase):
    def reject(self, read, priors, **rates):
        try:
            synchronize_read(read, priors, **(rates or RATES))
        except NotImplementedError:
            raise
        except ValueError:
            return
        self.fail("invalid or impossible synchronization input was accepted")

    def test_empty_source_and_read_stops_once(self):
        result = synchronize_read("", (), **RATES)
        self.assertAlmostEqual(result["log_likelihood"], math.log(0.8), places=13)
        self.assertEqual(result["posteriors"], ())
        self.assertEqual(result["log_evidence"], ())

    def test_insertion_only_has_trailing_stop_weight(self):
        result = synchronize_read("ATG", (), **RATES)
        expected = math.log(0.8) + 3 * math.log(0.2 / 4)
        self.assertAlmostEqual(result["log_likelihood"], expected, places=12)
        self.assertEqual(result["log_evidence"], ())

    def test_deleted_positions_have_neutral_extrinsic_evidence(self):
        result = synchronize_read("", PRIORS, **RATES)
        self.assertAlmostEqual(result["log_likelihood"],
                               math.log(0.8) + 2 * math.log(0.8 * 0.1), places=12)
        for actual, prior, evidence in zip(result["posteriors"], PRIORS,
                                           result["log_evidence"]):
            for value, wanted in zip(actual, prior):
                self.assertAlmostEqual(value, wanted, places=12)
            self.assertEqual(evidence, (0.0, 0.0, 0.0, 0.0))

    def test_noiseless_observation_has_exact_degenerate_posteriors(self):
        result = synchronize_read("AG", (UNIFORM,) * 2,
                                  p_insert=0, p_delete=0, p_substitute=0)
        self.assertAlmostEqual(result["log_likelihood"], math.log(1 / 16), places=12)
        for index, base in enumerate("AG"):
            wanted = tuple(float(candidate == base) for candidate in BASES)
            self.assertEqual(result["posteriors"][index], wanted)
            self.assertEqual(result["log_evidence"][index],
                             tuple(0.0 if candidate == base else -math.inf
                                   for candidate in BASES))

    def test_memoryless_evidence_removes_the_local_prior(self):
        results = [synchronize_read("T", (prior,), p_insert=0, p_delete=0,
                                    p_substitute=0.12) for prior in PRIORS]
        wanted = tuple(0.0 if base == "T" else math.log(0.04 / 0.88)
                       for base in BASES)
        for result, prior in zip(results, PRIORS):
            z = sum(p * (0.88 if b == "T" else 0.04) for b, p in zip(BASES, prior))
            for actual, expected in zip(result["log_evidence"][0], wanted):
                self.assertAlmostEqual(actual, expected, places=12)
            self.assertAlmostEqual(result["posteriors"][0][1],
                                   prior[1] * 0.88 / z, places=12)
        self.assertNotEqual(results[0]["posteriors"], results[1]["posteriors"])

    def test_fba_matches_exact_source_and_path_enumeration(self):
        for read in ("", "A", "AT", "CGA"):
            result = synchronize_read(read, PRIORS, **RATES)
            z, joint = enumerate_joint(read, PRIORS, RATES)
            self.assertAlmostEqual(result["log_likelihood"], math.log(float(z)),
                                   places=11)
            for i, row in enumerate(joint):
                evidence = [float(value / Fraction(str(PRIORS[i][b])))
                            for b, value in enumerate(row)]
                for b, value in enumerate(row):
                    self.assertAlmostEqual(result["posteriors"][i][b],
                                           float(value / z), places=11)
                    expected = math.log(evidence[b] / max(evidence))
                    self.assertAlmostEqual(result["log_evidence"][i][b],
                                           expected, places=11)

    def test_one_base_includes_both_insert_delete_orderings(self):
        result = synchronize_read("C", PRIORS[:1], **RATES)
        likelihoods = [2 * (0.2 / 4) * 0.8 * 0.1
                       + 0.8 * 0.9 * (0.85 if base == "C" else 0.05)
                       for base in BASES]
        z = 0.8 * sum(p * value for p, value in zip(PRIORS[0], likelihoods))
        self.assertAlmostEqual(result["log_likelihood"], math.log(z), places=12)
        for actual, value in zip(result["log_evidence"][0], likelihoods):
            self.assertAlmostEqual(actual, math.log(value / max(likelihoods)),
                                   places=12)

    def test_log_likelihood_survives_probability_underflow(self):
        prior = (1 - 3e-6, 1e-6, 1e-6, 1e-6)
        result = synchronize_read("C" * 256, (prior,) * 256,
                                  p_insert=0, p_delete=0, p_substitute=0)
        self.assertAlmostEqual(result["log_likelihood"], 256 * math.log(1e-6),
                               places=8)
        self.assertTrue(all(abs(row[3] - 1.0) < 1e-9
                            for row in result["posteriors"]))

    def test_uninformative_substitutions_leave_priors_unchanged(self):
        result = synchronize_read("ATG", PRIORS, p_insert=0.2, p_delete=0.1,
                                  p_substitute=0.75)
        for posterior, prior, evidence in zip(result["posteriors"], PRIORS,
                                              result["log_evidence"]):
            for value, wanted in zip(posterior, prior):
                self.assertAlmostEqual(value, wanted, places=12)
            self.assertTrue(all(abs(value) < 1e-12 for value in evidence))

    def test_work_counts_include_every_structural_term(self):
        result = synchronize_read("AGC", PRIORS, **RATES)
        self.assertEqual(result["work"],
                         {"lattice_cells": 12, "forward_transition_terms": 23,
                          "backward_transition_terms": 23, "extrinsic_terms": 56})

    def test_window_and_read_limits_are_rejected_before_allocation(self):
        self.reject("A", (UNIFORM,) * 257)
        self.reject("A" * 385, (UNIFORM,))

    def test_invalid_channel_rates_are_rejected(self):
        for name in RATES:
            for bad in (True, None, "0.1", -0.1, 1, math.inf, math.nan, 10**400):
                self.reject("A", (UNIFORM,), **dict(RATES, **{name: bad}))

    def test_invalid_prior_shape_or_mass_is_rejected(self):
        for bad in (None, [UNIFORM], ([0.25] * 4,), ((0.5, 0.5),),
                    ((0, 0.2, 0.3, 0.5),), ((-0.1, 0.3, 0.3, 0.5),),
                    ((True, 0.2, 0.3, 0.5),), ((math.nan, 0.2, 0.3, 0.5),),
                    ((math.inf, 0.2, 0.3, 0.5),), ((0.1, 0.2, 0.3, 0.3),),
                    (("0.1", 0.2, 0.3, 0.4),), ((10**400, 0.2, 0.3, 0.4),)):
            self.reject("A", bad)

    def test_invalid_read_alphabet_or_type_is_rejected(self):
        for bad in (None, b"AT", ["A"], "a", "N", "AT\n"):
            self.reject(bad, (UNIFORM,))

    def test_api_has_no_payload_alignment_or_scorer_parameter(self):
        synchronize_read("A", (UNIFORM,), **RATES)
        params = inspect.signature(synchronize_read).parameters
        self.assertEqual(tuple(params),
                         ("read", "base_priors", "p_insert", "p_delete", "p_substitute"))
        self.assertTrue(all(params[key].kind == inspect.Parameter.KEYWORD_ONLY
                            for key in RATES))

    def test_inputs_are_immutable_and_outputs_are_independent(self):
        snapshot = tuple(tuple(row) for row in PRIORS)
        first = synchronize_read("AT", PRIORS, **RATES)
        second = synchronize_read("AT", PRIORS, **RATES)
        self.assertEqual(PRIORS, snapshot)
        self.assertEqual(first, second)
        self.assertIsNot(first, second)
        self.assertIsInstance(first["log_evidence"], tuple)
        self.assertTrue(all(isinstance(row, tuple) for row in first["posteriors"]))

    def test_subnormal_substitution_evidence_is_not_silently_zero(self):
        rate = 1e-323
        result = synchronize_read("T", (UNIFORM,), p_insert=0, p_delete=0,
                                  p_substitute=rate)
        wanted = math.log(rate) - math.log(3)
        self.assertAlmostEqual(result["log_evidence"][0][0], wanted, places=11)
        self.assertTrue(math.isfinite(result["log_evidence"][0][0]))

    def test_zero_probability_observations_fail_closed(self):
        self.reject("AT", (UNIFORM,), p_insert=0, p_delete=0, p_substitute=0)
        self.reject("", (UNIFORM,), p_insert=0, p_delete=0, p_substitute=0)
