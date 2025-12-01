"""
NQX v1.0 hybrid decoder
Tek NQX (hybrid) -> tek FASTQ

Dosya yapısı:
[ PRE_HEADER (12B) ]
[ JSON_HEADER (json_length B, UTF-8) ]
[ GLOBAL_HEADER (binary) ]
[ BLOK(lar): BLOCK_HEADER + 4 adet zstd sıkıştırılmış stream ]

Bu dosya, encoder.py'nin NQX v1.0 hibrit sürümü ile üretilmiş .nqx dosyalarını
FASTQ'a geri çevirmek için referans decoder olarak tasarlanmıştır.
"""

from __future__ import annotations

import json
import struct
from typing import BinaryIO

import zstandard as zstd


# =========================
# Hibrit pre-header + JSON
# =========================

# Pre-header: "<4sII"
# 4s : PRE_MAGIC = b"NQJ1"
# I  : json_length (UTF-8 JSON header byte sayısı)
# I  : flags (şimdilik 0, gelecekte kullanım için ayrılmış)
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

PLATFORM_UNKNOWN = 0
PLATFORM_ILLUMINA = 1

QUAL_PHRED33 = 0
QUAL_PHRED64 = 1  # geleceğe yönelik, şimdilik kullanılmıyor

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
# Yardımcı fonksiyonlar
# =========================

def read_exact(f: BinaryIO, n: int) -> bytes:
    """Dosyadan tam n byte okur, eksikse hata fırlatır."""
    data = f.read(n)
    if len(data) != n:
        raise ValueError(f"Beklenen {n} byte, okunan {len(data)} byte (dosya bozuk veya eksik).")
    return data


def read_pre_and_json_header(f: BinaryIO) -> dict:
    """
    Hibrit NQX v1.0 pre-header + JSON header'ı okur ve JSON dict döner.
    Dosya pozisyonu bu fonksiyon sonrası GLOBAL_HEADER'ın başında olur.
    """
    pre_raw = read_exact(f, PRE_HEADER_STRUCT.size)
    magic, json_length, flags = PRE_HEADER_STRUCT.unpack(pre_raw)

    if magic != PRE_MAGIC:
        raise ValueError("PRE_MAGIC uyuşmuyor. Bu dosya hibrit NQX v1.0 formatında değil.")

    if json_length <= 0:
        raise ValueError("JSON header uzunluğu geçersiz.")

    json_bytes = read_exact(f, json_length)
    header_dict = json.loads(json_bytes.decode("utf-8"))

    # Basit doğrulamalar
    if header_dict.get("format") != "NQX":
        raise ValueError("JSON header 'format' alanı NQX değil.")

    if header_dict.get("version") != "1.0":
        raise ValueError("NQX version 1.0 bekleniyordu.")

    return header_dict


def read_global_header(f: BinaryIO) -> dict:
    """
    GLOBAL_HEADER_STRUCT'a göre global header'ı okur ve alanları dict olarak döner.
    Dosya pozisyonu bu fonksiyon sonrası ilk block header'ın başında olur.
    """
    header_bytes = read_exact(f, GLOBAL_HEADER_STRUCT.size)
    (
        magic,
        ver_major,
        ver_minor,
        platform,
        qual_enc,
        total_reads,
        total_bases,
        block_read_target,
    ) = GLOBAL_HEADER_STRUCT.unpack(header_bytes)

    if magic != MAGIC:
        raise ValueError("NQX global header magic uyuşmuyor (NQX1 bekleniyordu).")

    if ver_major != 1:
        raise ValueError(f"Desteklenmeyen NQX major versiyonu: {ver_major}")

    return {
        "version_major": ver_major,
        "version_minor": ver_minor,
        "platform": platform,
        "quality_encoding": qual_enc,
        "total_reads": total_reads,
        "total_bases": total_bases,
        "block_read_target": block_read_target,
    }


def decompress_stream(codec_id: int, comp_data: bytes) -> bytes:
    """
    Verilen codec_id'ye göre stream'i açar.
    Şimdilik:
      - CODEC_ID_RAW  -> no-op
      - CODEC_ID_ZSTD -> zstandard
    """
    if codec_id == CODEC_ID_RAW:
        return comp_data
    elif codec_id == CODEC_ID_ZSTD:
        zctx = zstd.ZstdDecompressor()
        return zctx.decompress(comp_data)
    else:
        raise ValueError(f"Desteklenmeyen codec_id: {codec_id}")


# =========================
# Ana decoder fonksiyonu
# =========================

def decode_nqx_to_fastq(nqx_path: str, fastq_path: str) -> None:
    """
    Hibrit NQX v1.0 dosyasını FASTQ formatına decode eder.

    Girdi:
      nqx_path  : encoder.py hibrit v1.0 tarafından üretilmiş .nqx dosyası
    Çıktı:
      fastq_path: FASTQ dosyası (round-trip için orijinal FASTQ ile birebir aynı olmalı)
    """
    with open(nqx_path, "rb") as fin, open(fastq_path, "wb") as fout:
        # 1) PRE_HEADER + JSON_HEADER
        json_header = read_pre_and_json_header(fin)

        # İstersen burada json_header içeriğini loglayabilir,
        # block_scheme'e göre davranış ayarlayabilirsin.
        # Şimdilik sadece doğrulama amaçlı okuyoruz.

        # 2) GLOBAL_HEADER (binary NQX payload başlangıcı)
        global_header = read_global_header(fin)

        # 3) Blokları sırayla oku
        while True:
            # Bir block header okumaya çalış
            block_header_bytes = fin.read(BLOCK_HEADER_STRUCT.size)
            if not block_header_bytes:
                # EOF (blok kalmadı)
                break

            if len(block_header_bytes) != BLOCK_HEADER_STRUCT.size:
                raise ValueError("Blok header eksik veya bozuk.")

            (
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
                codec_id_id,
                codec_id_plus,
                codec_id_seq,
                codec_id_qual,
            ) = BLOCK_HEADER_STRUCT.unpack(block_header_bytes)

            # Sıkıştırılmış stream'leri oku
            comp_id = read_exact(fin, comp_id_bytes)
            comp_plus = read_exact(fin, comp_plus_bytes)
            comp_seq = read_exact(fin, comp_seq_bytes)
            comp_qual = read_exact(fin, comp_qual_bytes)

            # Codec'e göre decompress et
            id_bytes = decompress_stream(codec_id_id, comp_id)
            plus_bytes = decompress_stream(codec_id_plus, comp_plus)
            seq_bytes = decompress_stream(codec_id_seq, comp_seq)
            qual_bytes = decompress_stream(codec_id_qual, comp_qual)

            # Uncompressed uzunlukları kontrol et
            if len(id_bytes) != uncomp_id_bytes:
                raise ValueError(
                    f"ID stream uncompressed size uyuşmuyor (beklenen={uncomp_id_bytes}, gerçek={len(id_bytes)})"
                )
            if len(plus_bytes) != uncomp_plus_bytes:
                raise ValueError(
                    f"PLUS stream uncompressed size uyuşmuyor (beklenen={uncomp_plus_bytes}, gerçek={len(plus_bytes)})"
                )
            if len(seq_bytes) != uncomp_seq_bytes:
                raise ValueError(
                    f"SEQ stream uncompressed size uyuşmuyor (beklenen={uncomp_seq_bytes}, gerçek={len(seq_bytes)})"
                )
            if len(qual_bytes) != uncomp_qual_bytes:
                raise ValueError(
                    f"QUAL stream uncompressed size uyuşmuyor (beklenen={uncomp_qual_bytes}, gerçek={len(qual_bytes)})"
                )

            # Newline bazlı read list'lerine ayrıştır
            id_lines = id_bytes.decode("utf-8").splitlines()
            plus_lines = plus_bytes.decode("utf-8").splitlines()
            seq_lines = seq_bytes.decode("utf-8").splitlines()
            qual_lines = qual_bytes.decode("utf-8").splitlines()

            if not (
                len(id_lines)
                == len(plus_lines)
                == len(seq_lines)
                == len(qual_lines)
                == block_read_count
            ):
                raise ValueError(
                    f"Blok {block_id} için satır sayıları uyuşmuyor: "
                    f"id={len(id_lines)}, plus={len(plus_lines)}, "
                    f"seq={len(seq_lines)}, qual={len(qual_lines)}, "
                    f"block_read_count={block_read_count}"
                )

            # FASTQ'ya yaz
            for i in range(block_read_count):
                # Orijinal encoder, PLUS satırını da saklıyordu,
                # bu yüzden round-trip'te aynen geri yazıyoruz.
                fout.write((id_lines[i] + "\n").encode("utf-8"))
                fout.write((seq_lines[i] + "\n").encode("utf-8"))
                fout.write((plus_lines[i] + "\n").encode("utf-8"))
                fout.write((qual_lines[i] + "\n").encode("utf-8"))


# =========================
# CLI
# =========================

if __name__ == "__main__":
    # python decoder.py test_R1_2.nqx decoded_R1.fastq
    import sys

    if len(sys.argv) != 3:
        print("Kullanım: python decoder.py <girdi.nqx> <cikti.fastq>")
        sys.exit(1)

    decode_nqx_to_fastq(sys.argv[1], sys.argv[2])
    print("NQX v1.0 hibrit decode tamamlandı.")
