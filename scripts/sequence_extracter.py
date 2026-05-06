import csv 
import os

binding_site_path = "../data/dosr_binding_sites.csv"
genome_path = "../data/GCF_000195955.2_ASM19595v2_genomic.fna"
output_path = "../data/sequences.fasta"

binding_sites = []
with open(binding_site_path, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        binding_sites.append((
            row["sequence"],
            int(row["start"]),
            int(row["end"]),
            row["strand"]
        ))


genome = ""
with open(genome_path, "r") as f:
    for line in f:
        if line.startswith(">"):
            pass
        else:
            genome += line.strip()


def negative_strand(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[b] for b in reversed(seq))


sequences = []
for i, (seq, start, end , strand) in enumerate(binding_sites):
    window_start = (start -1)- 250
    window_end = end + 250

    window_start = max(0, window_start)
    window_end = min(len(genome), window_end)

    window_sequence = genome[window_start:window_end]
    if strand == "-":
        window_sequence = negative_strand(window_sequence)
    sequences.append((seq, window_sequence))
    print(f"Site {i+1}: {seq} -> {len(window_sequence)}bp extracted")

with open(output_path, "w")as f:
    for i, (seq, window) in enumerate(sequences):
        f.write(window + "\n")
        
