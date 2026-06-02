from utils import Score
from entropyGibbsSample import entropyScore
known_sequences = [
        "TTAGGGCCGGAAGTCCCCAA",
        "TGGCGGACGAATGACCCCAG",
        "GCCGGGACTTCAGGCCCTAT",
        "TGGGGACCATTGACCCTG",
        "GCGGGCCATTTGTCCGCG",
        "TAGGGACAGAAGTCCCCG",
        "TCGGGGACTTCTGTCCCTAG",
        "ACAGGGTCAATGGTCCCCAA",
        "GTGGGGACCAACGCCCCTGG",
        "ATAAGGACTAACGGCCCTCA",
        "ACCTGGACGAGCCACCCGTG",
        "ACGGGATGTATCCGCCCCAG",
        "GGTCGGCCTTATGCCCCGTG",
        "TGCGGGTGGATCGGGCCATC"]
clean = [s for s in known_sequences if len(s) == 20]
print(f"Ground Truth Score: {Score(clean)}")
print(f"Ground Truth Entropy Score: {entropyScore(clean)}")