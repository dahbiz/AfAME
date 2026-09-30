#!/usr/bin/env python3
"""Extract every length-eight record from Danielsen's published GF(9) database.

Input is the unmodified selfdualcodes3.txt.bz2 from
https://www.codetables.de/larsed/nonbinary/selfdualcodes3.txt.bz2
No source code from that site is executed.
"""
import argparse
import bz2
import hashlib
import json
from pathlib import Path

URL = "https://www.codetables.de/larsed/nonbinary/selfdualcodes3.txt.bz2"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    selected, positions = [], []
    with bz2.open(args.archive, "rt", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            if line.split("\t", 1)[0] == "8":
                selected.append(line.rstrip("\r\n") + "\n")
                positions.append(line_number)
    if len(selected) != 817:
        raise ValueError(f"Expected all 817 length-eight classes, found {len(selected)}")
    payload = "".join(selected).encode("ascii")
    (args.output / "danielsen_length8.tsv").write_bytes(payload)
    manifest = {
        "source_url": URL,
        "source_archive_sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest(),
        "extracted_records_sha256": hashlib.sha256(payload).hexdigest(),
        "records": len(selected),
        "source_line_numbers": positions,
        "completeness_source": "Danielsen, Adv. Math. Commun. 3, 329-348 (2009); IEEE Trans. Inf. Theory 58, 5500-5511 (2012)",
        "record_format": "length, D/I, minimum distance, weight enumerator, automorphism group order, weighted edges; tab-separated",
    }
    (args.output / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in manifest.items() if k != "source_line_numbers"}, indent=2))


if __name__ == "__main__":
    main()
