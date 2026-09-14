"""Synthetic regression checks for the native-table command; no native inputs."""

import copy
import hashlib
import importlib
from pathlib import Path
import tempfile
import unittest

from a002_qualification import PchkMatrix


class NativeTableTests(unittest.TestCase):
    def setUp(self):
        self.driver = importlib.import_module("reproduce_native_table")

    def test_five_orders_use_inverse_direction_and_source_suffix(self):
        matrix = PchkMatrix(2, 6, ((0, 0), (0, 1), (1, 1), (1, 5)), 0, 0)
        rows = self.driver.ordering_rows(
            matrix, (1, 0, 1, 1, 0, 0), (0, 0, 0, 1), (1, 2, 0, 5, 3, 4))
        expected = {
            "stored": (1, 3),
            "undo_public_permutation": (0, 0),
            "apply_public_permutation": (2, 3),
            "rotate_2": (1, 2),
            "rotate_4": (1, 2),
        }
        self.assertEqual(
            {key: (value["syndrome_weight"], value["suffix_source_mismatches"])
             for key, value in rows.items()}, expected)

    def test_duplicate_permutation_is_rejected(self):
        matrix = PchkMatrix(2, 6, ((0, 0), (1, 1)), 0, 0)
        with self.assertRaises(ValueError):
            self.driver.ordering_rows(
                matrix, (0,) * 6, (0,) * 4, (0, 0, 2, 3, 4, 5))

    def test_wrong_source_geometry_is_rejected(self):
        matrix = PchkMatrix(2, 6, ((0, 0), (1, 1)), 0, 0)
        with self.assertRaises(ValueError):
            self.driver.ordering_rows(matrix, (0,) * 6, (0,) * 3, tuple(range(6)))

    def test_exact_file_is_read_without_modification(self):
        blob = b"0101"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            path.write_bytes(blob)
            self.assertEqual(self.driver.read_pinned_file(
                path, len(blob), hashlib.sha256(blob).hexdigest()), blob)
            self.assertEqual(path.read_bytes(), blob)

    def test_short_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            path.write_bytes(b"010")
            with self.assertRaisesRegex(ValueError, "size"):
                self.driver.read_pinned_file(path, 4, hashlib.sha256(b"0101").hexdigest())

    def test_trailing_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            path.write_bytes(b"01010")
            with self.assertRaisesRegex(ValueError, "size"):
                self.driver.read_pinned_file(path, 4, hashlib.sha256(b"0101").hexdigest())

    def test_same_size_wrong_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input"
            path.write_bytes(b"1010")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                self.driver.read_pinned_file(path, 4, hashlib.sha256(b"0101").hexdigest())

    def test_missing_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises((OSError, ValueError)):
                self.driver.read_pinned_file(Path(directory) / "missing", 4, "0" * 64)

    def test_expected_values_do_not_replace_computed_values(self):
        changed = copy.deepcopy(self.driver.EXPECTED_ROWS)
        changed["undo_public_permutation"]["syndrome_weight"] = 1
        preserved = copy.deepcopy(changed)
        with self.assertRaises(ValueError):
            self.driver.verify_expected_rows(changed)
        self.assertEqual(changed, preserved)

    def test_missing_candidate_is_not_an_exact_table_match(self):
        changed = copy.deepcopy(self.driver.EXPECTED_ROWS)
        del changed["stored"]
        with self.assertRaises(ValueError):
            self.driver.verify_expected_rows(changed)


if __name__ == "__main__":
    unittest.main()

