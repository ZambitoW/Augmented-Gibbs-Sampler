import numpy as np
from utils import ProfileBuilder, Score, profileRandomKmer, readSequences, compareToGroundTruth

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
        bestMotifs = gibbsSampler(sequences, k, t, N)
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