# Inspected upstream source

This note supports the observations about output scoring and stage progression
in the accompanying software article. It identifies twelve repository files
from one later snapshot of the watermark-assisted source system. It does not
assess the complete original release or establish equivalence to the code
archived with the publication.

Repository: dna-storage-lab/Bootstrap-readout-using-hidden-references.
Revision: ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b.
The [machine-readable inventory](inspected-source-files.json) records each
path, byte count, SHA-256 and Git blob identifier. The inventory includes
LICENSE and README.md; it is not a claim that twelve source programs were read.

## Prefix scoring

The Hamming-distance function iterates until its declared bit count is
reached, checks that each consumed character is binary and rejects a short
prefix. It then writes the distance and closes the files without checking
for a remaining suffix. Thus equality of the consumed prefixes does not
verify equality of the complete streams.

- [Scorer, lines 18-47][scorer].
- The three stage scripts pass 43,200 as the compared count and provide the
  supplied source bits: [stage I][stage1], [stage II][stage2], [stage III][stage3].

These are static observations of the named functions and invocations, not
results of executing the upstream programs. File extraction, byte order and
padding require their own contracts.

## Source-dependent stage progression

The wrapper reads each stage's reported distance. A zero first-stage distance
selects that result; otherwise it invokes the second stage. The second-stage
condition similarly controls whether the third stage is invoked.

- [Stage I decision and transition][stop1].
- [Stage II decision and transition][stop2].
- [Stage III result classification][stop3].

The visible LDPC invocations in the stage scripts above do not pass source
truth as an explicit decoder argument. That statement does not certify
hidden file access or the behavior of uninspected binaries.

This orchestration can serve as a known-answer evaluation harness. The
observation neither quantifies a recovery advantage nor implies misconduct.
A fixed-candidate Boolean comparison of any-stage and terminal-stage success
is an identity, not a replay of these programs or a measured performance gap.

## Complete inspected set

- [LICENSE](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/LICENSE)
- [README.md](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/README.md)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/src/decode_feedback_align.cpp](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/src/decode_feedback_align.cpp)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/src/CalBitError.c](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/src/CalBitError.c)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/build.sh](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/build.sh)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/recover.sh](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/recover.sh)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step1.sh](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step1.sh)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step2.sh](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step2.sh)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step3.sh](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step3.sh)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/src/post_dec_hamming_dis.c](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/src/post_dec_hamming_dis.c)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/src/Seq_chunker.c](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/src/Seq_chunker.c)
- [Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/src/SubsampleFastqRandom.c](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/src/SubsampleFastqRandom.c)

## Access and reuse

Links point to the pinned upstream files; their bodies are not copied here.
Consult the original source's terms before reuse. The repository-level MIT
notice does not establish permission for every externally hosted dataset.
The five configuration inputs used for the bit-order table are listed
separately in [INPUT_PROVENANCE.md](INPUT_PROVENANCE.md).

[scorer]: https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/src/post_dec_hamming_dis.c#L18-L47
[stage1]: https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step1.sh#L126-L145
[stage2]: https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step2.sh#L91-L105
[stage3]: https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/R0.67_bootstrap_recovery_step3.sh#L87-L101
[stop1]: https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/recover.sh#L95-L114
[stop2]: https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/recover.sh#L116-L137
[stop3]: https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/blob/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/recover.sh#L139-L154

