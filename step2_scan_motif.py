import re
import os
from multiprocessing import Pool, cpu_count
from itertools import islice

FASTA = "data/promoters_500bp.fa"
OUTFILE = "data/nrf1_dreg_hitss.dreg"
MOTIF = re.compile(r"GCGC..GCGC", re.IGNORECASE)
WORKERS = max(1, cpu_count() - 1)
CHUNK = 5000

def parse_fasta_chunks(path, chunk_size):
    chunk = []
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    chunk.append((header, "".join(seq_parts)))
                    if len(chunk) >= chunk_size:
                        yield chunk
                        chunk = []
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            chunk.append((header, "".join(seq_parts)))
    if chunk:
        yield chunk

def scan_chunk(records):
    hits = []
    for header, seq in records:
        for m in MOTIF.finditer(seq):
            hits.append(f"{header}\t{m.start()+1}\t{m.end()}\t{m.group()}\t+")
    return hits

if __name__ == "__main__":
    print(f"Scanning {FASTA} with {WORKERS} workers...")
    print("Counting sequences...", end=" ", flush=True)
    total = sum(1 for l in open(FASTA) if l.startswith(">"))
    print(f"{total:,} sequences found")

    done = 0
    all_hits = []

    with Pool(WORKERS) as pool:
        for hits in pool.imap_unordered(scan_chunk, parse_fasta_chunks(FASTA, CHUNK)):
            all_hits.extend(hits)
            done += CHUNK
            pct = min(100, done / total * 100)
            print(f"  {min(done, total):,}/{total:,} sequences ({pct:.1f}%) — {len(all_hits)} hits so far", flush=True)

    with open(OUTFILE, "w") as out:
        for h in all_hits:
            out.write(h + "\n")

    print(f"\nDone. {len(all_hits)} total hits written to {OUTFILE}")