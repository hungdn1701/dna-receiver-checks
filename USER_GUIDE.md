# DNA-storage receiver checks: user guide

This collection provides exact output scoring, representation checks and
reference soft-receiver components for a watermark-assisted DNA-storage design.
It is intended for researchers checking interfaces before comparing decoders.
It is not a sequencing pipeline or a reproduction of an upstream executable.

## Requirements and running the examples

The reference environment is Python 3.11.9. The six modules use only the Python
standard library. No package installation, network access or external dataset
is needed for the examples below. Run commands from the directory containing
the six `.py` modules and `tests/`. Use `python3` instead of `python` if that is
the name of your Python 3.11 interpreter.

```sh
python -B -m unittest discover -s tests -p 'test_a00*.py' -v
```

The pattern selects the seven receiver test modules, including when the
development repository contains other intake tests. A successful run ends with `OK`;
the test names identify each checked behavior. Failures and errors are not
successful decoding outcomes. Platform support beyond the recorded test
environment has not been established.

## Example 1: exact outputs and explicit byte order

```sh
python -B -m unittest discover -s tests -p test_a001_scoring.py -v
```

`parse_bits(text, expected_bits)` accepts bytes containing binary symbols and
ASCII whitespace. It requires the complete stream to contain exactly the
declared number of bits. The input `b"0001"` at length three raises `ValueError`,
as do a short stream and a nonbinary suffix. This constructed example explains
why prefix agreement differs from complete-stream agreement; it is not an
observed failure of another decoder.

`bits_to_bytes(bits, bit_order=...)` requires an immutable tuple of integer
bits and an explicit `"msb"` or `"lsb"` convention. The vector
`(1, 0, 0, 1, 0, 1, 1, 0)` produces hexadecimal `96` in most-significant-bit-first
order and `69` in least-significant-bit-first order. Non-byte-aligned inputs
are rejected, rather than padded.

`score_exact_bytes(candidate, truth)` compares the complete byte strings and
their lengths. Comparing `b"abcX"` with `b"abc"` returns `match=False`,
`actual_bytes=4`, `expected_bytes=3` and `byte_mismatches=None`. For equal-length
inputs the last field counts differing bytes, not differing bits. A valid
original-file assessment also requires a documented source-bit extraction,
packing and padding convention. The scorer cannot infer these from a codeword.

## Example 2: matrix and representation checks

```sh
python -B -m unittest discover -s tests -p test_a002_qualification.py -v
python -B -m unittest discover -s tests -p test_a001_mapping.py -v
```

`parse_pchk` reads the complete little-endian parity-check file, including
the `0x5080` marker and the sparse matrix body. Unknown markers, trailing
bytes and duplicate entries are rejected. `parse_sparse_matrix` is the
lower-level body reader; it reports trailing bytes and must not be confused
with the complete-file validator. Its size limit is 2 MiB for the sparse
body; callers should bound file reads before creating the input byte string.
The small fixture with matrix `H=[1,1]`
gives syndrome `(0,)` for codeword `(1,1)` and `(1,)` for `(1,0)`.

`check_syndrome` accepts parser-produced or directly constructed `PchkMatrix`
values. It validates integer dimensions in 1..64,800, tuple entries containing
integer row/column pairs, index bounds and uniqueness. Invalid matrices raise
`ValueError`; negative indices do not wrap. Unlike the encoder and decoder,
the syndrome checker also accepts square and tall matrices. Empty rows have
zero syndrome. Construction alone does not validate a `PchkMatrix`.

`parse_int32_candidates(blob, count)` returns only byte-order interpretations
that are complete permutations of `0,...,count-1`. An empty result means no
valid interpretation; two results mean ambiguity, not permission to choose
arbitrarily. The one-element zero array deliberately demonstrates ambiguity.

`encode_reference(codeword, permutation, watermark)` serializes bits as
`stored[i] = codeword[permutation[i]]`, applies the four-to-five sparse map
and the binary watermark, then maps two bit planes to DNA bases. The reverse
function `decode_reference` restores `codeword[permutation[i]] = stored[i]`
and checks sparse-pattern validity. Mapping tests include a non-involutive
permutation, nonzero watermark and independent quartet-table expectations.
A zero syndrome establishes code membership; without the source constraint
it does not establish which message was encoded.

### Relation to the published-size case study

The article's native tuple uses a supplied parity-check matrix, serialized
codeword, permutation, source bits and watermark. Those third-party input
bodies are not included in these synthetic examples. Reproducing that table
requires the permitted, pinned inputs and the separate native qualification
procedure. The small fixtures above do not reproduce the native table.

## Example 3: observation identity and soft evidence

```sh
python -B -m unittest discover -s tests -p test_a004_soft.py -v
```

`merge_evidence(packets, reference_length)` consumes tuples of
`(observation_id, start, rows)`. Starts are zero-based positions in a declared
reference. Each row contains log-likelihood evidence in **A, T, G, C** order,
with maximum zero. An uncovered position has row `(0,0,0,0)` and coverage zero.
Coverage counts included observation packets, not bases sequenced or molecules.
An empty packet counts as an observation but adds no evidence or coverage.

For row `(0,-1,-2,-3)`, sending the same packet twice under the same identifier
keeps that row unchanged and reports two observations seen, one unique
observation and one duplicate. Sending it under two distinct identifiers
produces `(0,-2,-4,-6)` with coverage two. Reusing an identifier with a different
start or row raises `ValueError`. The caller must preserve observation
identities consistently: identical sequence strings and distinct identifiers
do not establish molecular independence.

The tests also check overlapping evidence, conflicting exact observations,
and an eight-bit synthetic path connecting the receiver components.
This path uses a known reference span, public watermark priors, zero
insertion/deletion rates and one unique observation supplied twice. Its hard
decisions already satisfy parity, so the decoder exits without message-passing
iterations. It checks interface composition and duplicate handling. The
original payload is used to construct and assess the fixture, not as a
receiver-function argument. Passing the fixture is not a read-error-rate
measurement or proof of process-level truth isolation.

## Soft-receiver interfaces and limits

`synchronize_read(read, base_priors, p_insert=..., p_delete=...,
p_substitute=...)` uses a declared insertion/deletion/substitution channel.
The read contains uppercase A, T, G or C. Priors are positive probability
rows in A, T, G, C order, each summing to one within an absolute tolerance
of 1e-12. The returned `log_evidence` excludes the explicit
target-symbol prior; `posteriors` include that prior. Use `log_evidence`,
not the posterior rows, as input to evidence merging.

Synchronization accepts at most 256 reference positions and 384 observed
bases. It models consumption of the complete declared reference window and
observation, rather than locating an arbitrary read within a long reference.
The caller must choose the corresponding window; missing reference positions
are treated as deletions. Evidence
merging accepts at most 4,096 packets, 256 rows per packet and 40,500 reference
positions. Its checks reject contradictory exact evidence and finite
arithmetic overflow rather than treating them as neutral evidence.

Adding per-read evidence after marginalizing over neighbouring symbols is
an approximation. Shared sequence context can couple these messages even
for distinct observations; probability calibration under indels has not
been established. Observation deduplication alone does not resolve this.

`watermark_base_priors(watermark)` constructs factorized public priors.
`demap_codeword(log_evidence, permutation, watermark, *, clip=30.0)` returns
log-likelihood ratios with sign convention `log(P(bit=0)/P(bit=1))` in LDPC
order. The demapper separates the two DNA bit planes, so correlations between
them are discarded. Clipping must be positive and at most 30.

`encode_systematic(matrix, payload)` supports the documented unit
lower-bidiagonal parity block, not arbitrary LDPC generator construction.
`sum_product_decode(matrix, llrs, max_iterations, clip=30.0, early_stop=True)`
uses flooding updates. Its `parity_satisfied` result is a code-membership
check, not a guarantee of source or file recovery.
Both LDPC functions require integer dimensions satisfying
`0 < rows < cols <= 64800`; this restriction does not apply to `check_syndrome`.

## Interpretation and reuse

The exact scorer can be used separately from the DNA-specific components.
The matrix and mapping checks are useful when porting a serializer. The
evidence tests specify an interface for a separate alignment or localization
method, whose validity must be assessed independently. The package does not
automatically inspect a wrapper for source-dependent stopping, supply a
matched performance baseline, reconstruct an original file or establish
independent biological replicates. No external adoption is claimed.

## Distribution and input access

Version 0.1.0 is distributed under the MIT License. The maintainer is
Hung N. Dang (hungdn@ptit.edu.vn). See README.md for version information,
THIRD_PARTY_NOTICES.md for attribution and INPUT_PROVENANCE.md for the
external configuration tuple. The examples need only the included synthetic
fixtures; external native-input bodies are not redistributed.
