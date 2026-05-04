import csv
"""
Collected the data from:  http://www.collectf.org/browse/view_motif_reports_by_TF_and_species/4/25/
Since data not formatted in a downloadable way, had to manually create a dataset.

"""

binding_sites = [
    ("TTAGGGCCGGAAGTCCCCAA", 1965468, 1965487, "+"),
    ("TGGCGGACGAATGACCCCAG", 1965508, 1965527, "+"),
    ("GCCGGGACTTCAGGCCCTAT", 1965530, 1965549, "+"),
    ("TGGGGACCATTGACCCTG", 2279061, 2279078, "+"),
    ("GCGGGCCATTTGTCCGCG", 2279026, 2279043, "+"),
    ("TAGGGACAGAAGTCCCCG", 2279004, 2279021, "+"),
    ("TCGGGGACTTCTGTCCCTAG", 2279003, 2279022, "-"),
    ("ACAGGGTCAATGGTCCCCAA", 2279060, 2279079, "-" ),
    ("GTGGGGACCAACGCCCCTGG", 3500819, 3500838, "-"),
    ("ATAAGGACTAACGGCCCTCA", 3500840, 3500859, "-"),
    ("ACCTGGACGAGCCACCCGTG",3500911, 3500930, "-"),
    ("ACGGGATGTATCCGCCCCAG",3500761, 3500780, "-"),
    ("GGTCGGCCTTATGCCCCGTG",3500348, 3500367, "-" ),
    ("TGCGGGTGGATCGGGCCATC",3500084, 3500103, "-" )
]

output_file = "dosr_binding_sites.csv"

with open(output_file, "w", newline="")as f:
    writer = csv.writer(f)
    writer.writerow(["sequence", "start", "end", "strand"])
    for site in binding_sites:
        writer.writerow(site)

for seq, start, end, strand in binding_sites:
    print(f"{seq:<25} {start:<12}{end:<12}{strand}  ")