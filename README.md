# DNA-storage receiver checks

Reference Python components for checking bit ordering, output validation
and soft-evidence interfaces in a watermark-assisted DNA-storage receiver.

## Quick start

Use Python 3.11; the qualified reference environment is Python 3.11.9 on
Linux. No third-party package, network connection or external dataset is
needed to run the synthetic tests. From the repository root:

```sh
python -B -m unittest discover -s tests -v
```

The released source and seven test modules are byte-identical to the
qualified package, which passed 105 tests with no failures, errors or skips.
No new runtime test was performed solely for release packaging. A successful
run prints `OK`. Support on other interpreter versions and platforms has
not been established.

Read the [user guide](USER_GUIDE.md) for three worked examples, API contracts,
expected outcomes and resource limits.

## Components

| Module | Purpose |
|---|---|
| `exact_scoring.py` | Strict bit parsing, explicit byte packing, complete-byte comparison |
| `imeta_mapping.py` | Permutation, sparse mapping, watermark and DNA-base mapping |
| `a002_qualification.py` | Parity-check files, integer permutations and syndrome checks |
| `ldpc_reference.py` | Structured systematic encoding and flooding sum-product decoding |
| `sync_reference.py` | Bounded forward--backward synchronization |
| `soft_evidence.py` | Public priors, observation-ID merging and approximate soft demapping |

The module filenames retain the identifiers used during development so that
the tested source files and imports remain unchanged.

## Scope

This is reference software for interface qualification and method development.
The synchronization component needs a caller-specified reference window.
The tests use deterministic synthetic fixtures; they do not measure recovery
from sequencing reads or establish a performance gain. Native configuration
inputs used in the associated bit-order case study are not bundled.

The [input provenance](INPUT_PROVENANCE.md) identifies those external inputs
by repository revision, paths and checksums. The synthetic test suite is
self-contained and is not a reproduction of the native case-study table.

## Version, license and support

The initial release is `v0.1.0`. See [CHANGELOG.md](CHANGELOG.md) and
[CITATION.cff](CITATION.cff). Interfaces may change before version 1.0.

The original software is distributed under the [MIT License](LICENSE).
An identical [Licence.txt](Licence.txt) is provided for SoftwareX metadata
requirements. See [third-party attribution](THIRD_PARTY_NOTICES.md) for the
published mapping and input sources. The MIT grant for this repository does
not relicense external data or dependencies.

Maintainer: Hung N. Dang, Posts and Telecommunications Institute of Technology,
Hanoi, Vietnam. Contact: [hungdn@ptit.edu.vn](mailto:hungdn@ptit.edu.vn).

Code development and documentation used AI assistance. The software makes no
generative-model calls at runtime.
