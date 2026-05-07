"""
benchmark.py
-------------
Runs all Gibbs sampler versions and compares results.

Algorithms compared:
    1. Baseline Gibbs Sampler
    2. SA Gibbs Sampler (simulated annealing + exploration)
    3. Entropy Gibbs Sampler (entropy-guided selection)
    4. Augmented Gibbs Sampler (SA + entropy + exploration)

Usage:
    python scripts/benchmark.py
"""

import time
from utils import Score, readSequences, compareToGroundTruth
from gibbsSampler import gibbsSampler
from SAGibbsSampler import SAGibbsSampler
from entropyGibbsSample import entropyGibbsSampler
from augmentedGibbsSampler import augmentedGibbsSampler

# ============================================================
# PARAMETERS
# ============================================================

K    = 20
N    = 2000
RUNS = 30

# ============================================================
# RUNNER — wraps any sampler with multiple restarts
# ============================================================

def runSampler(name, samplerFn, sequences, k, t, N, runs):
    print(f"\n{'='*55}")
    print(f"  {name}")
    print(f"{'='*55}")

    scores = []
    allMotifs = []
    start_time = time.time()

    for r in range(runs):
        motifs = samplerFn(sequences, k, t, N)
        score = Score(motifs)
        scores.append(score)
        allMotifs.append(motifs)
        print(f"  Run {r+1:2d}: Score = {score}")

    elapsed = time.time() - start_time
    avg = sum(scores) / runs
    best = min(scores)
    worst = max(scores)
    bestMotifs = allMotifs[scores.index(best)]

    print(f"\n  Avg Score  : {avg:.2f}")
    print(f"  Best Score : {best}")
    print(f"  Worst Score: {worst}")
    print(f"  Time       : {elapsed:.1f}s")

    print(f"\n  Best Motifs Found:")
    for i, motif in enumerate(bestMotifs):
        print(f"    Site {i+1:2d}: {motif}")

    groundTruthPct = compareToGroundTruth(bestMotifs)

    return {
        "name"        : name,
        "avg"         : avg,
        "best"        : best,
        "worst"       : worst,
        "scores"      : scores,
        "bestMotifs"  : bestMotifs,
        "groundTruth" : groundTruthPct,
        "time"        : elapsed,
    }

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    sequences = readSequences("data/sequences.fasta")
    t = len(sequences)

    print(f"Dataset   : DosR M. tuberculosis H37Rv")
    print(f"Sequences : {t}")
    print(f"Motif k   : {K}")
    print(f"Iterations: {N}")
    print(f"Runs      : {RUNS}")

    results = []
    results.append(runSampler("Baseline Gibbs",   gibbsSampler,          sequences, K, t, N, RUNS))
    results.append(runSampler("SA Gibbs",         SAGibbsSampler,        sequences, K, t, N, RUNS))
    results.append(runSampler("Entropy Gibbs",    entropyGibbsSampler,   sequences, K, t, N, RUNS))
    results.append(runSampler("Augmented Gibbs",  augmentedGibbsSampler, sequences, K, t, N, RUNS))

    # ============================================================
    # SUMMARY TABLE
    # ============================================================

    print(f"\n{'='*65}")
    print(f"FINAL SUMMARY")
    print(f"{'='*65}")
    print(f"{'Algorithm':<22} {'Avg':>6} {'Best':>6} {'Worst':>6} {'GT Match':>10} {'Time':>8}")
    print(f"{'-'*65}")
    for r in results:
        print(f"{r['name']:<22} {r['avg']:>6.2f} {r['best']:>6} {r['worst']:>6} {r['groundTruth']:>9.1f}% {r['time']:>7.1f}s")

    print(f"\nLower score = better | Higher GT Match % = better")

    best_score_result = min(results, key=lambda x: x["avg"])
    best_truth_result = max(results, key=lambda x: x["groundTruth"])

    print(f"\nBest avg score    : {best_score_result['name']} ({best_score_result['avg']:.2f})")
    print(f"Best ground truth : {best_truth_result['name']} ({best_truth_result['groundTruth']:.1f}%)")