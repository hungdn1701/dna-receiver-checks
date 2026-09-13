# Source attribution and third-party notices

## Watermark-assisted mapping

The permutation, four-to-five sparse mapping, watermark XOR and two-plane
base mapping in `imeta_mapping.py` are an independently written implementation
of the mapping described in `process_r083` in the following source:

- Repository: [Bootstrap-readout-using-hidden-references](https://github.com/dna-storage-lab/Bootstrap-readout-using-hidden-references).
- Inspected revision: `ddf71a0e6db8c5faf328f7a990d1020aba9cbc0b`.
- Source file: `decode_feedback_align.cpp`, Git blob
  `0d70cb9de75c5b6126e18f06e3dce8b0da4572e9`.
- Related publication: [Chen et al., DOI 10.1002/imt2.70105](https://doi.org/10.1002/imt2.70105).

No upstream binary, third-party library or complete receiver implementation
is included. The inverse mapping, strict scorer and validation tests are
independently implemented in this package. The inspected snapshot has not
been established to be identical to the publication-time archive.

The upstream MIT notice is retained below as attribution. It does not assert
that all third-party material in the upstream repository has the same license.

## Other implementation sources

The sparse-file marker was checked against Radford Neal's LDPC-codes
`rcode.c`; no source body from that file is included or executed.
The encoder, decoder, synchronizer and evidence interfaces are independent
reference implementations. No external package is required at runtime.

## External input material

Native configuration inputs and sequencing data are excluded from this
distribution. [INPUT_PROVENANCE.md](INPUT_PROVENANCE.md) records the inspected
input locations and hashes. Obtain external materials from their source and
review applicable terms; this repository does not grant rights to them.

## Retained upstream license notice

```text
MIT License

Copyright (c) 2025 dna-storage-lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
