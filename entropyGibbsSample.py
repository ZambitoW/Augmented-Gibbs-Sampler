import numpy as np
from utils import ProfileBuilder, profileRandomKmer, readSequences, compareToGroundTruth
import math
"""
Entropy-based Gibbs Sampler
Uses entropy as the scoring function instead of Hamming distance.
Lower entropy = more conserved motif set = better.
"""
def entropyScore(motifs):
    DNA = "ACGT" 
    k = len(motifs[0])
    t = len(motifs)
    total_entropy = 0.0

    for col in range(k):
        counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        for motif in motifs:
            counts[motif[col]] +=1
        for base in DNA:
            p = counts[base]/t
            if p>0:
                total_entropy -= p * math.log2(p)
    return total_entropy


def entropyGibbsSampler(Dna, k, t, N):
    motifs = []
    for sequence in Dna:
        start = np.random.randint(0, len(sequence) - k + 1)
        motifs.append(sequence[start: start + k])
    bestMotifs = motifs[:]
    bestScore = entropyScore(bestMotifs)

    for j in range(1, N):
        i = np.random.randint(0, t)
        profile = ProfileBuilder(motifs[:i] + motifs[i+1:], k)
        motifi = profileRandomKmer(Dna[i], profile, k)
        motifs[i] = motifi
        current_score = entropyScore(motifs)
        if current_score < bestScore:
            bestMotifs = motifs[:]
            bestScore = current_score

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
        score = entropyScore(bestMotifs)
        scores.append(score)
        print(f"Run {r+1}: Entropy Score = {score:.4f}")

        if score < overallBestScore:
            overallBestScore = score
            overallBestMotifs = bestMotifs[:]

    avg_score = sum(scores) / runs
    print(f"\nAverage Entropy Score : {avg_score:.4f}")
    print(f"Best Entropy Score    : {overallBestScore:.4f}")

    print("\nBest Motifs Found:")
    for i, motif in enumerate(overallBestMotifs):
        print(f"  Site {i+1}: {motif}")

    compareToGroundTruth(overallBestMotifs)