import gzip
import csv

input_file = "data/human_gene_annotation.tsv.gz"
tss_bed = "data/tss.bed"
chrom_sizes_file = "data/hg38.chrom.sizes"
promoter_bed = "data/promoters_500bp.bed"
SLOP = 500

standard_chroms = {f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]}

# Load chrom sizes for clamping
chrom_sizes = {}
with open(chrom_sizes_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            chrom_sizes[parts[0]] = int(parts[1])

seen = set()
rows = []

with gzip.open(input_file, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        chrom_raw = row["chromosome_name"].strip()
        chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"
        if chrom not in standard_chroms:
            continue
        if chrom not in chrom_sizes:
            continue
        tss = int(row["transcription_start_site"])
        gene = row["external_gene_name"].strip()
        strand_raw = row["strand"].strip()
        strand = "+" if strand_raw == "1" else "-"
        start = tss
        end = tss + 1
        name = f"{chrom}@{start}-{end}|{gene}"
        key = (chrom, start, gene)
        if key in seen:
            continue
        seen.add(key)
        rows.append((chrom, start, end, name, ".", strand, chrom_sizes[chrom]))

rows.sort(key=lambda x: (x[0], x[1]))

with open(tss_bed, "w") as out:
    for r in rows:
        out.write("\t".join(str(x) for x in r[:6]) + "\n")

print(f"Written {len(rows)} TSS entries to {tss_bed}")

# Build promoter bed directly with clamping (no bedtools slop needed)
skipped = 0
written = 0
with open(promoter_bed, "w") as out:
    for chrom, start, end, name, dot, strand, chromlen in rows:
        if strand == "+":
            pstart = max(0, start - SLOP)
            pend = end
        else:
            pstart = start
            pend = min(chromlen, end + SLOP)
        if pstart >= pend:
            skipped += 1
            continue
        out.write(f"{chrom}\t{pstart}\t{pend}\t{name}\t.\t{strand}\n")
        written += 1

print(f"Written {written} promoter entries to {promoter_bed} ({skipped} skipped)")
print("\nNow run: bash step2_extract_and_dreg.sh")