# Reproduce the native bit-order table

The synthetic test suite is self-contained. This separate command reproduces
the five rows of the native configuration comparison from the supplied input
tuple. It does not run the upstream receiver, acquire sequencing reads or
reconstruct an original file.

## 1. Obtain the exact inputs

Use the original filenames in a new directory. The five files, sizes,
checksums and upstream terms are documented in [INPUT_PROVENANCE.md](INPUT_PROVENANCE.md).
The selected revision is
`ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b`, not a moving branch.

If access and reuse are permitted under the applicable source terms, the
following shell commands download only those five files. They stop if the
destination directory already exists or a request fails; they do not
execute anything from the upstream repository.

```sh
set -e
mkdir native-inputs
base='https://raw.githubusercontent.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/configure'
for name in SequenceLengthALL_FILE001R0667 dvb_s2_r2_3.pchk encoded_bit.txt permutation64800 source_bit.txt
do
    curl --fail --silent --show-error --max-time 60 "$base/$name" --output "native-inputs/$name"
done
```

Alternatively, download those files from the pinned GitHub links and retain
their original bytes. Do not copy text from the rendered GitHub preview or
normalize line endings. The runner validates byte counts and SHA-256 before
parsing. The download snippet requires curl; the Python runner does not.

## 2. Run the offline comparison

From the software repository root, using Python 3.11:

```sh
python -B reproduce_native_table.py --assets native-inputs
```

The five input files are read, not changed. The command prints JSON to
standard output and exits zero only if all input pins, geometry checks and
computed rows match. Missing, modified, short or extended input files cause
a nonzero exit. Supply ordinary files, not symlinks.

## 3. Inspect the expected result

The JSON `ordering_checks` object must contain:

| Candidate | Syndrome weight | Source-suffix mismatches |
|---|---:|---:|
| stored | 10988 | 21458 |
| undo_public_permutation | 0 | 0 |
| apply_public_permutation | 10658 | 21420 |
| rotate_21600 | 10779 | 21666 |
| rotate_43200 | 10805 | 21521 |

The output also contains `expected_rows_match: true`, the input checksums,
and matrix dimensions 21600 by 64800 with 215999 nonzero entries.
The command computes the rows from the input bytes before comparing them
with the expected table. It does not print stored expectations as results.

For stored bits s and permutation p, the inverse candidate is defined by
x[p[i]] = s[i]. The forward candidate is x[i] = s[p[i]]. Rotations are left
cyclic shifts of the stored word. Each distance compares the last 43200 bits
with the supplied source; each syndrome is Hx over GF(2).

Four files enter that arithmetic. The watermark file is checksum-checked to
identify the complete original five-file tuple; this command does not run
the watermark mapper or the native-geometry soft-receiver fixture.

## Interpretation

Only the inverse convention satisfies both constraints among these five
inspected alternatives. This does not establish uniqueness over all
permutations. Source-suffix agreement and parity satisfaction do not prove
complete original-file recovery. The command is a deterministic
reproduction check, not a performance benchmark or new sequencing result.

The [source provenance note](SOURCE_PROVENANCE.md) separately documents the
static observations about the upstream scorer and stage policy. Those
observations are not outputs of this command.

