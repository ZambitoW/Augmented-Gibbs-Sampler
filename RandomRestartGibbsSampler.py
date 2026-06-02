import numpy as np
from utils import ProfileBuilder, Score, profileRandomKmer, readSequences

def gibbsSampler(Dna, k, t, N):
    motifs = []
    for sequence in Dna:
        start = np.random.randint(0, len(sequence) - k + 1)
        motifs.append(sequence[start:start + k])
    bestMotifs = motifs[:]
    for j in range(1, N):
        i = np.random.randint(0, t)
        profile = ProfileBuilder(motifs[:i] + motifs[i+1:], k)
        motifi = profileRandomKmer(Dna[i], profile, k)
        motifs[i] = motifi
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs[:]
    return bestMotifs

def randomRestartsSampler(Dna, k, t, N, runs):
    overallBestMotifs = None
    overallBestScore = float("inf")
    scores = []

    for r in range(runs):
        bestMotifs = gibbsSampler(Dna, k, t, N)
        score = Score(bestMotifs)
        scores.append(score)
        print(f"Run {r+1}: Score = {score}")

        if score < overallBestScore:
            overallBestMotifs = bestMotifs[:]
            overallBestScore = score
    avgScore = sum(scores) / runs
    return overallBestMotifs, overallBestScore, avgScore

if __name__ == "__main__":
    k = 20
    N = 2000
    runs = 30
    sequences = readSequences("data/sequences.fasta")
    t = len(sequences)
    bestMotifs, bestScore, avgScore = randomRestartsSampler(
        sequences, k, t, N, runs
    )
    print("\nBest Motifs Across All Runs:")
    for i, motif in enumerate(bestMotifs):
        print(f"  Site {i+1}: {motif}")
    print(f"\nBest Score: {bestScore}")
    print(f"Average Score over {runs} runs: {avgScore:.2f}")
