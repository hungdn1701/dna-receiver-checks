# External configuration inputs

The synthetic test suite runs without external files. The native bit-order
case study instead uses five configuration and known-answer assets from
[the pinned upstream revision](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references/tree/ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b/Recovery_code/Figure5/bootstrap_recovery_TypeI%2BII%2BIIIReads_FBA/R0.67/configure).

Repository: `dna-storage-lab/Bootstrap-readout-using-hidden-references`.
Revision: `ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b`.
Path prefix:
`Recovery_code/Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/configure/`.

These are not sequencing reads. Their bodies are not redistributed here.
The upstream root has an MIT notice (retained in THIRD_PARTY_NOTICES.md);
users should check the applicable source terms before obtaining or reusing
any external material.

## Exact inputs

| File | Bytes | Role |
|---|---:|---|
| `SequenceLengthALL_FILE001R0667` | 324000 | Binary watermark |
| `dvb_s2_r2_3.pchk` | 950412 | Parity-check matrix |
| `encoded_bit.txt` | 64800 | Serialized known-answer word |
| `permutation64800` | 259200 | Permutation |
| `source_bit.txt` | 43200 | Known source suffix |

SHA-256 values of the original bytes:

```text
02a34a578244eeea001fb9d82c47848bd09c65f304bce1282e9b25b188c7af4f  SequenceLengthALL_FILE001R0667
25b638c5e10ba2f436fa6b8d4493b64c91e1fdb28d2b2c7fa845a7f9d9e5b9b1  dvb_s2_r2_3.pchk
bc42f06c5b6ff6a6a569bbcb2575365f9a95003b672ea117413ba5973c2cdfac  encoded_bit.txt
42f751c6d95c782ccaa3171aacf354435f77bc4186e433f52e89827995416742  permutation64800
0e8a9103ee89bb1689e0d54f263ee562c19be332115d7bdb7fa943b1f753853d  source_bit.txt
```

The pinned path and checksums identify the inspected tuple, not every
possible input of the upstream system. Later-snapshot equivalence to the
publication-time archive has not been established.

## Relation to the case study

The case study compares the stored word, inverse-permuted word,
forward-permuted word, and left cyclic shifts by 21,600 and 43,200 bits.
For each word it evaluates the syndrome weight over GF(2) and the Hamming
distance of the last 43,200 bits from the supplied source. The inverse
permutation is defined by `codeword[permutation[i]] = stored[i]`.

The library provides the parsing, mapping and syndrome primitives for such
a comparison. This initial release does not include a native-data acquisition
or case-study runner. The synthetic tests demonstrate the interfaces and
must not be presented as reproducing that native table.
