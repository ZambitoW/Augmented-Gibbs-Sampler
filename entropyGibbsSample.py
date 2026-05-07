import numpy as np
from utils import ProfileBuilder, Score, profileRandomKmer, readSequences, compareToGroundTruth
import math

def motifEntropy(motif,profile, k):
    DNA = "ACGT"
    entropy = 0.0
    for j in range(k):
        for bp in DNA:
            p = profile[DNA.find(bp)][j]
            if p>0 :
                entropy -= p * math.log2(p)
    return entropy



def entropyGibbsSampler(Dna, k, t, N):
    motifs = []
    for sequence in Dna:
        start = np.random.randint(0, len(sequence) -k + 1)
        motifs.append(sequence[start: start + k])
    bestMotifs = motifs[:]

    for j in range(1, N):
        entropies = []
        for index in range(t):
            profile_index = ProfileBuilder(motifs[:index] + motifs[index+1:], k)
            entropy = motifEntropy(motifs[index], profile_index, k)
            entropies.append(entropy) 
        
        entropies = np.array(entropies)
        max_entropy = np.max(entropies)
        weights = np.exp(entropies - max_entropy)
        weights /= weights.sum()
        i=np.random.choice(t,p=weights)

        profile = ProfileBuilder(motifs[:i] + motifs[i+1:], k)
        motifi = profileRandomKmer(Dna[i], profile ,k )
        motifs[i] =motifi
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs[:]
    return bestMotifs



if __name__ == "__main__":
    k = 20
    N = 2000
    sequences = readSequences("data/sequences.fasta")
    t = len(sequences)
 
    runs = 30
    scores = []
    overallBestMotifs = None
    overallBestScore = float("inf")
 
    for r in range(runs):
        bestMotifs = entropyGibbsSampler(sequences, k, t, N)
        score = Score(bestMotifs)
        scores.append(score)
        print(f"Run {r+1}: Score = {score}")
 
        if score < overallBestScore:
            overallBestScore = score
            overallBestMotifs = bestMotifs[:]
 
    avg_score = sum(scores) / runs
    print(f"\nAverage Score : {avg_score:.2f}")
    print(f"Best Score    : {overallBestScore}")
 
    print("\nBest Motifs Found:")
    for i, motif in enumerate(overallBestMotifs):
        print(f"  Site {i+1}: {motif}")
 
    compareToGroundTruth(overallBestMotifs)