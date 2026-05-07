import numpy as np
import math
from utils import ProfileBuilder, profileRandomKmer, Score, readSequences, compareToGroundTruth

def motifEntropy(motif, profile, k):
    DNA = "ACGT"
    entropy = 0.0
    for j in range(k):
        for bp in DNA:
            p = profile[DNA.find(bp)][j]
            if p>0:
                entropy -= p *math.log2(p)

    return entropy

def augmentedGibbsSampler(Dna,k ,t, N):
    motifs = []
    for sequence in Dna:
        start = np.random.randint(0, len(sequence) - k +1)
        motifs.append(sequence[start: start + k])
    bestMotifs = motifs[:]
    bestScore = Score(bestMotifs)

    T = 5.0 * k
    cooling_rate = math.exp(math.log(0.01) / N)
    exploration_rate = 0.2

    for j in range(1, N):
        entropies = []
        for index in range(t):
            profile_index = ProfileBuilder(motifs[:index] + motifs[index+1:], k)
            ent = motifEntropy(motifs[index], profile_index,k)
            entropies.append(ent)
        entropies = np.array(entropies)
        weights = np.exp(entropies - np.max(entropies))
        weights /= weights.sum()
        i = np.random.choice(t, p = weights)

        profile = ProfileBuilder(motifs[:i] + motifs[i +1 :], k)
        if np.random.random()< exploration_rate:
            start = np.random.randint(0, len(Dna[i]) - k +1)
            motifi = Dna[i][start:start+k]
        else:
            motifi = profileRandomKmer(Dna[i], profile, k)

        candidate_motifs = motifs[:]
        candidate_motifs[i] = motifi
        delta = Score(candidate_motifs) - Score(motifs)

        if delta < 0 or np.random.random() < math.exp(-delta/T):
            motifs[i] = motifi

        current_score = Score(motifs)
        if current_score < bestScore:
            bestMotifs = motifs[:]
            bestScore = current_score
        T*= cooling_rate
    return bestMotifs

def runAugmented(Dna, k, t, N, runs):
    overallBestMotifs = None
    overallBestScore = float("inf")
    scores =[]

    for r in range(runs):
        motifs = augmentedGibbsSampler(Dna, k , t, N)
        score = Score(motifs)
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