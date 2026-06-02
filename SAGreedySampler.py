import numpy as np
import math
from utils import ProfileBuilder, profileRandomKmer, Score, readSequences, compareToGroundTruth

"""
This Algorithm starts with GreedyMotifSearch as implemented on rosalind and then runs GibbsSampling on it.
It also has Simulated Annealing as this was show to be the most beneficial tool to improve Gibbs. 
Only uses 1 run because Greedy always give you same starting motif, so for time reasons only once. 
"""

def profileMostProbable(sequence, k, profile):
    DNA = "ACGT"
    bestProb = -1
    bestKmer = sequence[:k]
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        prob = 1.0
        for j, base in enumerate(kmer):
            prob *= profile[DNA.find(base)][j]
        if prob > bestProb:
            bestProb = prob
            bestKmer = kmer
    return bestKmer
 
 
def greedyMotifSearch(Dna, k, t):
    bestMotifs = [dna[:k] for dna in Dna]
    for i in range(len(Dna[0]) - k + 1):
        firstMotif = Dna[0][i:i+k]
        motifs = [firstMotif]
        for j in range(1, t):
            profile = ProfileBuilder(motifs, k)
            currentBestMotif = profileMostProbable(Dna[j], k, profile)
            motifs.append(currentBestMotif)
 
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

def SAGibbsSampler(Dna, k, t, N):
    motifs = greedyMotifSearch(Dna, k, t)
    bestMotifs = motifs[:]
    bestScore = Score(bestMotifs)

    T = 5.0*k
    cooling_rate = math.exp(math.log(0.01) / N)
    exploration_rate = 0.20 
    for j in range(1, N):
        i = np.random.randint(0, t)
        profile = ProfileBuilder(motifs[:i] + motifs[i+1:], k)
        if np.random.random() < exploration_rate:
            start = np.random.randint(0, len(Dna[i]) - k +1)
            motifi = Dna[i][start: start+k]
        else: 
            motifi = profileRandomKmer(Dna[i], profile ,k )

        candidate_motifs= motifs[:]
        candidate_motifs[i]=motifi
        delta = Score(candidate_motifs) - Score(motifs)

        if delta<0 or np.random.random() < math.exp(-delta/T):
            motifs[i]= motifi
        current_score = Score(motifs)
        if current_score< bestScore:
            bestMotifs = motifs[:]
            bestScore = current_score
        
        T *= cooling_rate
    return bestMotifs



if __name__ == "__main__":
    k = 20
    N = 2000
    sequences = readSequences("data/sequences.fasta")
    t = len(sequences)
 
    runs = 1
    scores =[]
    overallBestMotifs =None
    overallBestScore = float("inf")
 
    for r in range(runs):
        bestMotifs = SAGibbsSampler(sequences, k, t,N)
        score = Score(bestMotifs)
        scores.append(score)
        print(f"Run {r+1}: Score = {score}")
 
        if score < overallBestScore:
            overallBestScore = score
            overallBestMotifs = bestMotifs[:]
 
    avg_score = sum(scores) / runs
    print(f"\nAverage Score: {avg_score:.2f}")
    print(f"Best Score: {overallBestScore}")
 
    print("\nBest Motifs Found:")
    for i, motif in enumerate(overallBestMotifs):
        print(f"Site {i+1}: {motif}")
 
    compareToGroundTruth(overallBestMotifs)