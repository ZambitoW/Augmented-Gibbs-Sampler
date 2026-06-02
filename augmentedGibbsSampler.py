import numpy as np
import math
from utils import ProfileBuilder, profileRandomKmer, Score, readSequences, compareToGroundTruth

"""
This model is the augmented one from our proposal that combines both the entropy
and simulated annealing models. 
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

def augmentedGibbsSampler(Dna,k ,t, N):
    motifs = []
    for sequence in Dna:
        start = np.random.randint(0, len(sequence) - k + 1)
        motifs.append(sequence[start: start + k])
    bestMotifs = motifs[:]
    bestScore = Score(bestMotifs)

    T = 5.0 * k
    cooling_rate = math.exp(math.log(0.01) / N)
    exploration_rate = 0.2

    for j in range(1, N):
        i = np.random.randint(0,t)
        profile = ProfileBuilder(motifs[:i] + motifs[i+1:], k)
        if np.random.random() < exploration_rate:
            start = np.random.randint(0, len(Dna[i]) - k + 1)
            motifi= Dna[i][start:start+k]
        else:
            motifi= profileRandomKmer(Dna[i], profile, k)
    
    #MH acceptance
        candidate_motifs = motifs[:]
        candidate_motifs[i] =motifi
        delta = entropyScore(candidate_motifs) - entropyScore(motifs)
        if delta < 0 or np.random.random() < math.exp(-delta/T):
            motifs[i] = motifi
    
        current_score = entropyScore(motifs)
        if current_score< bestScore:
            bestMotifs = motifs[:]
            bestScore = current_score

        T *= cooling_rate

    return bestMotifs
        

def runAugmented(Dna, k, t, N, runs):
    overallBestMotifs = None
    overallBestScore = float("inf")
    scores =[]

    for r in range(runs):
        motifs = augmentedGibbsSampler(Dna, k , t, N)
        score = entropyScore(motifs)
        scores.append(score)

        print(f"Run {r+1}: Score = {score}")
        if score< overallBestScore:
            overallBestScore = score
            overallBestMotifs = motifs
        
    avgScore = sum(scores) / runs
    return overallBestMotifs, overallBestScore, avgScore

if __name__ == "__main__":
    k = 20
    N= 2000
    runs = 30
    sequences = sequences = readSequences("data/sequences.fasta")
    t = len(sequences)

    bestMotifs, bestScore, avgScore = runAugmented(sequences, k,t,N, runs)

    print("\nBest Motifs Across All Runs:")
    for i, motif in enumerate(bestMotifs):
        print (f"  Site{i+1}: {motif}")
    print(f"\n Best Score: {bestScore}")
    print(f"Avg Score: {avgScore:.2f}")

    compareToGroundTruth(bestMotifs)