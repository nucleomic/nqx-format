
from __future__ import annotations

import json
import struct
from pathlib import Path
from typing import Iterator, Tuple, BinaryIO

import zstandard as zstd  

PRE_HEADER_STRUCT = struct.Struct("<4sII")
PRE_MAGIC = b"NQJ1"


# =========================
# Binary NQX payload sabitleri
# =========================

# Global header: "<4sBBBBQQI"
# 4s : MAGIC = b"NQX1"
# B  : version_major
# B  : version_minor
# B  : platform  (0=unknown, 1=illumina, 2=ont, 3=pb, ...)
# B  : quality_encoding (0=phred33, 1=phred64)
# Q  : total_reads
# Q  : total_bases
# I  : block_read_target
GLOBAL_HEADER_STRUCT = struct.Struct("<4sBBBBQQI")

MAGIC = b"NQX1"
VERSION_MAJOR = 1
VERSION_MINOR = 0

PLATFORM_UNKNOWN = 0
PLATFORM_ILLUMINA = 1

QUAL_PHRED33 = 0
QUAL_PHRED64 = 1  
# Block header: "<IIIIIIIIIIBBBB"
# I : block_id
# I : block_read_count
# I : uncomp_id_bytes
# I : uncomp_plus_bytes
# I : uncomp_seq_bytes
# I : uncomp_qual_bytes
# I : comp_id_bytes
# I : comp_plus_bytes
# I : comp_seq_bytes
# I : comp_qual_bytes
# B : codec_id_id
# B : codec_id_plus
# B : codec_id_seq
# B : codec_id_qual
BLOCK_HEADER_STRUCT = struct.Struct("<IIIIIIIIIIBBBB")

CODEC_ID_RAW = 0
CODEC_ID_ZSTD = 1


# =========================
# FASTQ yardımcıları
# =========================

def iter_fastq_reads(f: BinaryIO) -> Iterator[Tuple[str, str, str, str]]:
    while True:
        id_line = f.readline()
        if not id_line:
            return  

        seq_line = f.readline()
        plus_line = f.readline()
        qual_line = f.readline()

        if not seq_line or not plus_line or not qual_line:
            raise ValueError("FASTQ dosyası bozuk: 4 satırlık blok tamamlanamadı.")

        id_str = id_line.decode("utf-8").rstrip("\r\n")
        seq_str = seq_line.decode("utf-8").rstrip("\r\n")
        plus_str = plus_line.decode("utf-8").rstrip("\r\n")
        qual_str = qual_line.decode("utf-8").rstrip("\r\n")

        if not id_str.startswith("@"):
            raise ValueError(f"FASTQ ID satırı '@' ile başlamıyor: {id_str!r}")

        if len(seq_str) != len(qual_str):
            raise ValueError(
                f"SEQ ve QUAL uzunlukları uyuşmuyor (len(seq)={len(seq_str)}, len(qual)={len(qual_str)})"
            )

        yield id_str, seq_str, plus_str, qual_str


# =========================
# Encoder ana fonksiyonu
# =========================

def encode_fastq_to_nqx(
    fastq_path: str | Path,
    nqx_path: str | Path,
    block_size: int = 5000,
    platform: int = PLATFORM_ILLUMINA,
    quality_encoding: int = QUAL_PHRED33,
) -> None:
    """
    Verilen FASTQ dosyasını NQX v1.0 hibrit formatına dönüştürür.

    Dosya yapısı:
    [ PRE_HEADER (12B) ]
    [ JSON_HEADER (json_length B, UTF-8) ]
    [ GLOBAL_HEADER (binary) ]
    [ BLOK(lar): BLOCK_HEADER + 4 adet zstd sıkıştırılmış stream ]
    """
    fastq_path = Path(fastq_path)
    nqx_path = Path(nqx_path)

    if block_size <= 0:
        raise ValueError("block_size pozitif bir tamsayı olmalıdır.")

   
    header_dict = {
        "format": "NQX",
        "variant": "hybrid",
        "version": "1.0",
        "source": "fastq",
        "platform": "illumina",
        "quality_encoding": "phred33",
        "block_scheme": {
            "type": "4_stream",
            "streams": ["id", "plus", "seq", "qual"],
            "codec": "zstd",
            "block_read_target": block_size,
        },
    }

    json_bytes = json.dumps(header_dict, separators=(",", ":")).encode("utf-8")
    json_length = len(json_bytes)

   
    zstd_level = 3
    zctx = zstd.ZstdCompressor(level=zstd_level)

    total_reads = 0
    total_bases = 0
    block_id = 0

    with fastq_path.open("rb") as fin, nqx_path.open("wb") as fout:
        
        pre_header = PRE_HEADER_STRUCT.pack(PRE_MAGIC, json_length, 0)
        fout.write(pre_header)
        fout.write(json_bytes)

       
        global_header_pos = fout.tell()
        dummy_global_header = GLOBAL_HEADER_STRUCT.pack(
            MAGIC,
            VERSION_MAJOR,
            VERSION_MINOR,
            platform,
            quality_encoding,
            0, 
            0, 
            block_size,
        )
        fout.write(dummy_global_header)

        
        reader = iter_fastq_reads(fin)

        while True:
            ids: list[str] = []
            pluses: list[str] = []
            seqs: list[str] = []
            quals: list[str] = []

            
            for _ in range(block_size):
                try:
                    id_str, seq_str, plus_str, qual_str = next(reader)
                except StopIteration:
                    break

                ids.append(id_str)
                pluses.append(plus_str)
                seqs.append(seq_str)
                quals.append(qual_str)

                total_reads += 1
                total_bases += len(seq_str)

            if not ids:
                
                break

            block_read_count = len(ids)

            
            id_bytes = ("\n".join(ids) + "\n").encode("utf-8")
            plus_bytes = ("\n".join(pluses) + "\n").encode("utf-8")
            seq_bytes = ("\n".join(seqs) + "\n").encode("utf-8")
            qual_bytes = ("\n".join(quals) + "\n").encode("utf-8")

            uncomp_id_bytes = len(id_bytes)
            uncomp_plus_bytes = len(plus_bytes)
            uncomp_seq_bytes = len(seq_bytes)
            uncomp_qual_bytes = len(qual_bytes)

            
            comp_id = zctx.compress(id_bytes)
            comp_plus = zctx.compress(plus_bytes)
            comp_seq = zctx.compress(seq_bytes)
            comp_qual = zctx.compress(qual_bytes)

            comp_id_bytes = len(comp_id)
            comp_plus_bytes = len(comp_plus)
            comp_seq_bytes = len(comp_seq)
            comp_qual_bytes = len(comp_qual)

            
            block_header = BLOCK_HEADER_STRUCT.pack(
                block_id,
                block_read_count,
                uncomp_id_bytes,
                uncomp_plus_bytes,
                uncomp_seq_bytes,
                uncomp_qual_bytes,
                comp_id_bytes,
                comp_plus_bytes,
                comp_seq_bytes,
                comp_qual_bytes,
                CODEC_ID_ZSTD,
                CODEC_ID_ZSTD,
                CODEC_ID_ZSTD,
                CODEC_ID_ZSTD,
            )
            fout.write(block_header)

            
            fout.write(comp_id)
            fout.write(comp_plus)
            fout.write(comp_seq)
            fout.write(comp_qual)

            block_id += 1

       
        end_pos = fout.tell()
        fout.seek(global_header_pos)
        final_global_header = GLOBAL_HEADER_STRUCT.pack(
            MAGIC,
            VERSION_MAJOR,
            VERSION_MINOR,
            platform,
            quality_encoding,
            total_reads,
            total_bases,
            block_size,
        )
        fout.write(final_global_header)
        fout.seek(end_pos)


# =========================
# CLI
# =========================

if __name__ == "__main__":
    import sys

    if not (3 <= len(sys.argv) <= 4):
        print("Kullanım: python encoder.py <girdi.fastq> <cikti.nqx> [block_size]")
        sys.exit(1)

    in_fastq = sys.argv[1]
    out_nqx = sys.argv[2]
    if len(sys.argv) == 4:
        try:
            block_size_val = int(sys.argv[3])
        except ValueError:
            print("block_size tamsayı olmalıdır.")
            sys.exit(1)
    else:
        block_size_val = 5000

    encode_fastq_to_nqx(in_fastq, out_nqx, block_size=block_size_val)
    print("NQX v1.0 hibrit encode tamamlandı.")
