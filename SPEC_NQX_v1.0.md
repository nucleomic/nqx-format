# NQX v1.0 Hybrid File Format Specification

**Status:** Draft / Reference Specification  
**Version:** 1.0  
**Author:** İlteriş Eren AMİL  
**Maintainer:** İlteriş Eren AMİL  
**Contact:** contact@nucleomic.com

---

## 1. Introduction

NQX is a binary container format designed to store sequencing reads that are
logically equivalent to FASTQ data but in a more compact, structured, and
versioned layout.

NQX v1.0 (“hybrid”) consists of:

1. A fixed-size **binary pre-header**
2. A structured **JSON metadata header**
3. A binary **global header**
4. A sequence of **blocks**, each containing:
   - A binary block header
   - Four per-stream compressed payloads:
     - ID
     - PLUS
     - SEQ
     - QUAL

All multi-byte integers in NQX v1.0 are **little-endian**.

The design goal is full round-trip fidelity:

> FASTQ → NQX → FASTQ MUST reproduce the original FASTQ byte-for-byte.

NQX v1.0 is specified for FASTQ (Phred+33, Illumina-style) reads.

---

## 2. File Extension and Versioning

Recommended extension:  
- `.nqx` (raw NQX v1.0 file)  
- `.nqx.gz` or `.nqx.zst` (optional external compression)

Versioning fields appear in both:
- The JSON header
- The binary global header

Reader rules:
- If pre-header `magic != "NQJ1"` → not hybrid NQX v1.0  
- If global header `magic != "NQX1"` → unsupported payload version

---

## 3. Overall File Layout

The logical structure of an NQX v1.0 file is:

```
+---------------------------+ offset 0
| PRE_HEADER (12 bytes)     |
+---------------------------+
| JSON_HEADER (UTF-8)       |
+---------------------------+
| GLOBAL_HEADER (binary)    |
+---------------------------+
| BLOCK_HEADER (block 0)    |
| BLOCK_STREAMS (block 0)   |
+---------------------------+
| BLOCK_HEADER (block 1)    |
| BLOCK_STREAMS (block 1)   |
+---------------------------+
| ...                       |
```

---

## 4. Binary Pre-header (12 bytes)

C-like structure:

```c
struct PreHeader {
    char   magic[4];      // MUST be "NQJ1"
    uint32 json_length;   // length in bytes of the JSON header
    uint32 flags;         // reserved, MUST be 0 in v1.0
};
```

### Field definitions

- **magic**  
  MUST be the ASCII bytes:
  ```
  4E 51 4A 31   ("NQJ1")
  ```

- **json_length**  
  Byte-length of the following UTF-8 JSON header.  
  MUST be > 0.

- **flags**  
  Reserved for future use.  
  For NQX v1.0, MUST be `0`.

If the pre-header cannot be fully read, or if the fields are invalid,  
a reader MUST reject the file.

---

## 5. JSON Metadata Header

Immediately following the pre-header, the file contains a UTF-8 JSON object
with exactly `json_length` bytes.

### Required JSON structure

```
{
  "format": "NQX",
  "variant": "hybrid",
  "version": "1.0",
  "source": "fastq",
  "platform": "illumina",
  "quality_encoding": "phred33",
  "block_scheme": {
    "type": "4_stream",
    "streams": [ "id", "plus", "seq", "qual" ],
    "codec": "zstd",
    "block_read_target": 5000
  }
}
```

### Required Fields

- `"format"`: MUST be `"NQX"`
- `"variant"`: `"hybrid"` for NQX v1.0
- `"version"`: MUST be `"1.0"`
- `"source"`: MUST be `"fastq"` for v1.0 reference encoder
- `"platform"`: `"illumina"`
- `"quality_encoding"`: MUST be `"phred33"`

### Block Scheme

- `"type"` MUST be `"4_stream"`
- `"streams"` MUST be `[ "id", "plus", "seq", "qual" ]`
- `"codec"` MUST be `"zstd"` for all streams in v1.0 reference implementation
- `"block_read_target"`: target number of reads per block

Readers MUST fail for malformed JSON.

---

## 6. Global Binary Header

Immediately after the JSON header:

```c
struct GlobalHeader {
    char   magic[4];          // "NQX1"
    uint8  version_major;     // MUST be 1
    uint8  version_minor;     // MUST be 0
    uint8  platform;          // 0=unknown, 1=illumina
    uint8  quality_encoding;  // 0=phred33
    uint64 total_reads;       // total FASTQ reads
    uint64 total_bases;       // sum of SEQ lengths
    uint32 block_read_target; // same as JSON.block_scheme.block_read_target
};
```

### Field definitions

- **magic**  
  MUST be `"NQX1"`.

- **version_major / version_minor**  
  MUST be `1.0` for this spec.

- **platform**  
  - `0` = unknown  
  - `1` = Illumina  
  (future revisions may add more)

- **quality_encoding**  
  - `0` = Phred+33

- **total_reads**  
  Sum of all block_read_count.

- **total_bases**  
  Sum of lengths of all SEQ lines.

- **block_read_target**  
  Recommended reads per block.

If the global header cannot be fully read or contains invalid fields,
reader MUST reject the file.

---

## 7. Block Structure

Each block consists of:

```
[ BLOCK_HEADER ] [ COMP_ID ] [ COMP_PLUS ] [ COMP_SEQ ] [ COMP_QUAL ]
```

The block header describes sizes and codec IDs.  
Each compressed stream follows immediately after the header.

---

## 8. Block Header Format

C-like layout:

```c
struct BlockHeader {
    uint32 block_id;
    uint32 block_read_count;

    uint32 uncomp_id_bytes;
    uint32 uncomp_plus_bytes;
    uint32 uncomp_seq_bytes;
    uint32 uncomp_qual_bytes;

    uint32 comp_id_bytes;
    uint32 comp_plus_bytes;
    uint32 comp_seq_bytes;
    uint32 comp_qual_bytes;

    uint8  codec_id_id;
    uint8  codec_id_plus;
    uint8  codec_id_seq;
    uint8  codec_id_qual;
};
```

### Field definitions

- **block_id**  
  Zero-based block index.

- **block_read_count**  
  Number of FASTQ records in this block.

- **uncomp_*_bytes**  
  Uncompressed byte lengths of each stream.

- **comp_*_bytes**  
  Compressed byte lengths.

- **codec_id_\***  
  - `0` = raw (no compression)  
  - `1` = Zstandard (zstd)  
  v1.0 reference encoder uses **1** for all streams.

---

## 9. Compressed Streams

Four streams appear in the exact order:

1. COMP_ID  
2. COMP_PLUS  
3. COMP_SEQ  
4. COMP_QUAL

A reader MUST:

1. Read exactly `comp_*_bytes` for each stream  
2. Decompress using `codec_id_*`  
3. Validate that decompressed size equals `uncomp_*_bytes`

Mismatch ⇒ file is corrupt.

---

## 10. FASTQ Reconstruction Rules

After decompression:

```
id_lines   = splitlines(ID_STREAM)
plus_lines = splitlines(PLUS_STREAM)
seq_lines  = splitlines(SEQ_STREAM)
qual_lines = splitlines(QUAL_STREAM)
```

Constraints:

- All four arrays MUST be equal length.  
- Length MUST equal `block_read_count`.

Reconstruction:

For each index `i`:

```
id_lines[i]   + "\n"
seq_lines[i]  + "\n"
plus_lines[i] + "\n"
qual_lines[i] + "\n"
```

Writer must preserve content exactly (no trimming, no reformatting).

---

## 11. Compression Rules

### Codec IDs

- `0` = RAW  
- `1` = Zstandard (zstd)

NQX v1.0 reference encoder uses:

```
codec_id_id   = 1
codec_id_plus = 1
codec_id_seq  = 1
codec_id_qual = 1
```

Compression level is implementation-defined and not stored.

Readers MUST support both codec IDs.

---

## 12. Validation Rules (strict)

A reader MUST reject a file if:

### Pre-header
- magic != "NQJ1"  
- json_length == 0  
- premature EOF

### JSON
- not valid UTF-8  
- missing required fields  
- version != "1.0"

### Global header
- magic != "NQX1"  
- version_major != 1  

### Blocks
- premature EOF  
- mismatched uncompressed sizes  
- codec_id not recognized  
- number of lines != block_read_count  

Readers MAY additionally validate:

- sum(block_read_count) == total_reads  
- sum(SEQ line lengths) == total_bases  

---

## 13. Reference Implementation

The following Python files represent the official reference implementation for NQX v1.0 hybrid:

- `encoder.py`  
- `decoder.py`  
- `nqx-info.py` (metadata inspector)

These tools implement exactly the format defined in this document.

---

## 14. Test Vectors

Round-trip test:

```
python encoder.py test.fastq test.nqx
python decoder.py test.nqx test_roundtrip.fastq
```

FASTQ files MUST match byte-for-byte.

External compression (optional):

```
gzip -9 test.fastq
gzip -9 test.nqx
```

NQX.gz SHOULD be comparable or smaller, depending on data.

---

## 15. Future Directions (NQX v2+)

Possible future revisions:

- Additional layout types  
- Domain-specific compression  
- Parallel block encoding metadata  
- Extended JSON fields (sample metadata, read groups)  
- Support for multi-file NQX containers  

Compatibility rules:

- Pre-header and JSON `"version"` MUST accurately describe new variants  
- Readers MUST explicitly check version compatibility  

---

## 16. License

This specification and the reference implementation are released under the MIT License.

```text
Copyright (c) 2025
İlteriş Eren AMİL
```

See `LICENSE` file for full terms.

---

## 17. Citation

If you use the NQX v1.0 format or its reference implementation, please cite:

İlteriş Eren AMİL.  
**NQX v1.0 Hybrid File Format Specification.**  
Zenodo. 2025.  
https://doi.org/10.5281/zenodo.17772957

---