import sys
import numpy as np

def ProfileBuilder(motifs, k):
    DNA = "ACGT"
    profile = np.ones((4,k))

    for i in range(k):
        for motif in motifs:
            row = DNA.find(motif[i])
            profile[row][i] +=1        
    profile = profile / (len(motifs) + 4)

    return profile


def Score(motifs):
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
        "TGCGGGTGGATCGGGCCATC",
    ]
    
    print(f"\n{'Site':<6} {'Found':<25} {'Known':<25} {'Match %'}")
    print("-" * 70)
    
    total_match = 0
    for i, (found, known) in enumerate(zip(foundMotifs, known_sequences)):

        matches = sum(f == k for f, k in zip(found, known))
        percentage = matches / max(len(found), len(known)) * 100
        total_match += percentage
        print(f"{i+1:<6} {found:<25} {known:<25} {percentage:.1f}%")
    print(f"\nOverall match: {total_match/len(foundMotifs):.1f}%")
