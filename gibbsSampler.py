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


def gibbsSampler(Dna, k, t, N):
    motifs = []
    for sequence in Dna:
        start = np.random.randint(0, len(sequence) -k + 1)
        motifs.append(sequence[start: start + k])
    bestMotifs = motifs[:]
    for j in range(1, N):
        i = np.random.randint(0, t)
        profile = ProfileBuilder(motifs[:i] + motifs[i+1:], k)
        motifi = profileRandomKmer(Dna[i], profile ,k )
        motifs[i] =motifi
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs[:]
    return bestMotifs

def readSequences(filename):
    sequences = []
    with open(filename,"r") as f:
        for line in f:
            sequences.append(line.strip())
    return sequences

if __name__ == "__main__":
    k = 20
    N = 500
    sequences = readSequences("data/sequences.fasta")
    t= len(sequences)

    bestMotifs = gibbsSampler(sequences, k,t,N)

    print(f"\nBest Motifs Found:")
    for i, motif in enumerate(bestMotifs):
        print(f"  Site {i+1}: {motif}")
    print(f"\nScore: {Score(bestMotifs)}")



