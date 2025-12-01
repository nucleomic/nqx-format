"""
nqx-info.py
NQX v1.0 Hybrid Format – Metadata Inspector
"""

import json
import struct
from pathlib import Path

PRE_HEADER_STRUCT = struct.Struct("<4sII")
PRE_MAGIC = b"NQJ1"

GLOBAL_HEADER_STRUCT = struct.Struct("<4sBBBBQQI")
MAGIC = b"NQX1"


def read_exact(f, n):
    data = f.read(n)
    if len(data) != n:
        raise ValueError(f"Expected {n} bytes, got {len(data)} bytes.")
    return data


def print_nqx_info(path: str):
    p = Path(path)
    with p.open("rb") as f:

        # --- PRE HEADER ---
        pre = read_exact(f, PRE_HEADER_STRUCT.size)
        magic, json_len, flags = PRE_HEADER_STRUCT.unpack(pre)

        print("PRE HEADER")
        print("  magic      :", magic)
        print("  json_len   :", json_len)
        print("  flags      :", flags)
        print()

        if magic != PRE_MAGIC:
            raise ValueError("Invalid PRE_MAGIC (not NQX v1.0 hybrid).")

        # --- JSON HEADER ---
        json_bytes = read_exact(f, json_len)
        j = json.loads(json_bytes.decode("utf-8"))

        print("JSON HEADER")
        print(json.dumps(j, indent=2))
        print()

        # --- GLOBAL HEADER ---
        gh_raw = read_exact(f, GLOBAL_HEADER_STRUCT.size)
        (
            magic2,
            ver_major,
            ver_minor,
            platform,
            qual_enc,
            total_reads,
            total_bases,
            block_target,
        ) = GLOBAL_HEADER_STRUCT.unpack(gh_raw)

        print("GLOBAL HEADER")
        print("  magic          :", magic2)
        print("  version        :", f"{ver_major}.{ver_minor}")
        print("  platform       :", platform)
        print("  quality_enc    :", qual_enc)
        print("  total_reads    :", total_reads)
        print("  total_bases    :", total_bases)
        print("  block_target   :", block_target)
        print()


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python nqx-info.py <file.nqx>")
        raise SystemExit
    print_nqx_info(sys.argv[1])
