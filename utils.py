import sys
import numpy as np
import math
from scipy.optimize import linear_sum_assignment

def ProfileBuilder(motifs, k):
    """
    This is from our rosalind assignments, reimplemented for practice
    """
    DNA = "ACGT"
    profile = np.ones((4,k))
    for i in range(k):
        for motif in motifs:
            row = DNA.find(motif[i])
            profile[row][i] +=1        
    profile = profile / (len(motifs) + 4)
    return profile


def Score(motifs):
    """
    This is from our rosalind assignments, reimplemented for practice
    """
    counts = {}
    score = 0
    for col in range(len(motifs[0])):
        counts = {"A" :0, "C" :0, "G": 0, "T":0}

        for motif in motifs:
            counts[motif[col]] +=1
        maxCount = max(counts.values())
        score += (len(motifs) - maxCount)
    return score



def profileRandomKmer(sequence, profile, k):
    """
    This is from our rosalind assignments, reimplemented for practice
    """
    probs = []
    DNA = "ACGT"
    for i in range(len(sequence) -k +1):
        kmer = sequence[i:i+k]
        prob = 1.0
        for j, letter in enumerate(kmer):
            prob *= profile[DNA.find(letter)][j]
        probs.append(prob)

    probs = np.array(probs)
    probs /= np.sum(probs)

    index = np.digitize(np.random.random(), np.cumsum(probs))
    return sequence[index: index+k]

def readSequences(filename):
    sequences = []
    with open(filename,"r") as f:
        for line in f:
            sequences.append(line.strip())
    return sequences



def compareToGroundTruth(foundMotifs):
    """
    This method allows us to compare our found motifs to the ground truth, since score
    is not super relevant. We care about the actual binding sites from CollectTf
    """

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
    
    matrix = np.array([
        [sum(f==k for f,k in zip(found, known)) / max(len(found),len(known)) * 100
        for known in known_sequences]
        for found in foundMotifs])

    row_ind, col_ind = linear_sum_assignment(-matrix)  

    print(f"\n{'Site':<6} {'Known':<25} {'Best Match':<25} {'Match %'}")
    print("-" * 70)
    total = 0
    for r, c in zip(row_ind, col_ind):
        pct = matrix[r, c]
        total += pct
        print(f"{c+1:<6} {known_sequences[c]:<25} {foundMotifs[r]:<25} {pct:.1f}%")

    print(f"\nOverall match: {total/len(known_sequences):.1f}%")


