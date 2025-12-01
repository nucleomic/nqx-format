# NQX v1.0 Hybrid Format  
**Reference Specification and Implementation**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17772957.svg)](https://doi.org/10.5281/zenodo.17772957)


NQX is a compact, lossless, block-based container format for FASTQ sequencing data.  
It combines a lightweight JSON header with a structured binary payload to provide:

- Full FASTQ round-trip fidelity  
- Independent compression of ID, PLUS, SEQ and QUAL streams  
- Block-based layout enabling future advanced compression models  
- Clear versioning and extensible metadata  

This repository contains the **official NQX v1.0 hybrid specification** and the  
**reference Python encoder/decoder**.

---

## Features

- Hybrid structure: JSON metadata + binary block payload  
- Fully lossless (byte-identical FASTQ → NQX → FASTQ)  
- Zstandard-compressed 4-stream architecture  
- Integer-safe, little-endian binary headers  
- Extensible for future compression schemes (NQX v2+)  
- Simple, readable, versioned metadata via JSON header  

---

## File Structure (High-Level)

```
+--------------------------+
| PRE_HEADER (12 bytes)    |  -- binary: "NQJ1", json_length, flags
+--------------------------+
| JSON_HEADER (UTF-8)      |  -- metadata & layout description
+--------------------------+
| GLOBAL_HEADER (binary)   |  -- version, counts, block size, etc.
+--------------------------+
| BLOCK HEADER 0           |
| BLOCK STREAMS 0          |
| BLOCK HEADER 1           |
| BLOCK STREAMS 1          |
| ...                      |
+--------------------------+
```

Full byte-level format description is available in:

**`SPEC_NQX_v1.0.md`**

---

## Reference Tools

This repository includes:

| File            | Description                                     |
|-----------------|-------------------------------------------------|
| `encoder.py`    | Reference encoder (FASTQ → NQX v1.0)            |
| `decoder.py`    | Reference decoder (NQX v1.0 → FASTQ)            |
| `nqx-info.py`   | Metadata inspector for NQX files                |
| `SPEC_NQX_v1.0.md` | Official file format specification          |

All tools require **Python 3.8+** and:

```
pip install zstandard
```

---

## Usage

### Encode FASTQ → NQX

```
python encoder.py input.fastq output.nqx
```

Optional block size (default 5000 reads):

```
python encoder.py input.fastq output.nqx 8000
```

---

### Decode NQX → FASTQ

```
python decoder.py input.nqx output.fastq
```

---

### Inspect metadata (JSON + PRE_HEADER + GLOBAL_HEADER)

```
python nqx-info.py file.nqx
```

Example output:

```
PRE HEADER
  magic      : b'NQJ1'
  json_len   : 231
  flags      : 0

JSON HEADER
{
  ... metadata ...
}

GLOBAL HEADER
  magic          : b'NQX1'
  version        : 1.0
  total_reads    : 128000
  total_bases    : 38400000
  block_target   : 5000
```

---

## Validation (Round-Trip Test)

```
python encoder.py test.fastq test.nqx
python decoder.py test.nqx test_roundtrip.fastq
```

Compare:

- Windows:

```
fc test.fastq test_roundtrip.fastq
```

- Linux/Mac:

```
diff test.fastq test_roundtrip.fastq
```

**The files MUST be identical.**

---

## Example Dataset

The `examples/` directory may include:

- `test_R1.fastq`  
- `test_R1.nqx`  
- `test_R1_roundtrip.fastq`  

These serve as official test vectors for verifying implementations.

---

## Versioning

- This repository defines **NQX v1.0 hybrid**.  
- Future versions (NQX v2+) may introduce new block layouts and codec stacks.  
- Backward compatibility is guaranteed through the metadata version fields in the JSON header.
- NQX file format version: 1.0 (SPEC_NQX_v1.0.md)
- GitHub/Zenodo software release: v1.0.1 (DOI: 10.5281/zenodo.17772957)


---

## Developed By

**İlteriş Eren AMİL**  
Contact: contact@nucleomic.com  
Website: https://nucleomic.com


---

## License

MIT License

---

## Citation

If you use the NQX format or the reference implementation in your work, please cite:

İlteriş Eren AMİL.  
**NQX v1.0 Hybrid File Format Specification.**  
Zenodo. 2025.  
https://doi.org/10.5281/zenodo.17772957

---