# Enhanced Gibbs Sampling for Transcription Factor Binding Site Discovery

**CS 321 · Middlebury College · Spring 2025**

Implemented and evaluated four augmented versions of the Gibbs Sampler algorithm for DNA motif finding, applied to real *M. tuberculosis* genomic data. The goal was to reduce convergence to local optima and evaluate which enhancements improve biological accuracy.

> Built with Owen Thornton

---

## Background

Motif finding identifies short, recurring DNA patterns such as transcription factor binding sites. The classic **Gibbs Sampler** is a randomized algorithm that iteratively updates motifs using a probability profile — but it frequently gets stuck in local optima.

This project implements and compares four enhancements:

| Algorithm | Key Idea |
|---|---|
| Baseline Gibbs | Standard Rosalind implementation — random init, unconditional acceptance |
| SA Gibbs | Simulated Annealing with Metropolis-Hastings acceptance + random k-mer injection |
| Greedy + SA | Exhaustive greedy warm start followed by SA with lower temperature |
| Entropy | Shannon entropy score instead of Hamming distance |
| Augmented (SA + Entropy) | Combines SA, entropy scoring, and random k-mer injection |

---

## Dataset

- **14 DosR binding sites** from *M. tuberculosis* H37Rv genome (CollecTF + NCBI)
- **14 × 520bp** extracted input sequences
- **20bp ground truth** binding sites
- Parameters: `k=20`, `N=2000`, `30 restarts`

---

## Results

| Algorithm | Avg Score (↓) | Avg Entropy (↓) | GT Match % (↑) |
|---|---|---|---|
| Ground Truth | 81 | 23.65 | 100% |
| Baseline Gibbs | 76.10 | — | 23.1% |
| SA Gibbs | 71.20 | — | 61.4% |
| Greedy + SA | **63** | — | 33.9% |
| Entropy | — | **15.15** | 63.2% |
| Augmented (SA + Entropy) | — | 15.51 | **63.2%** |

**Key findings:**
- SA Gibbs improved ground truth match from 23.1% → 61.4% over baseline
- Greedy + SA achieved the best Hamming score (63) but worst GT match (33.9%) — greedy initialization biased search toward early local optima
- Entropy and Augmented models peaked at 63.2% GT match but over-conserved flanking positions, converging on a more specific k-mer than the true binding sites
- Conservation-based objectives cannot reliably distinguish the true motif from other highly conserved sequences in 520bp windows — a position-specific scoring matrix anchored to known sites would likely improve accuracy

---

## File Structure

```
Augmented-Gibbs-Sampler/
├── gibbsSampler.py            # Baseline Gibbs Sampler
├── SAGibbsSampler.py          # Simulated Annealing variant
├── SAGreedySampler.py         # Greedy warm start + SA
├── entropyGibbsSample.py      # Entropy scoring variant
├── augmentedGibbsSampler.py   # Full augmented (SA + Entropy)
├── RandomRestartGibbsSampler.py
├── groundTruthScore.py        # GT match evaluation
├── utils.py                   # Shared utilities
├── data/                      # M. tuberculosis genome + binding sites
├── scripts/                   # Dataset construction scripts
└── CS321_Final.pdf            # Conference-style poster
```

---

## Setup

```bash
pip install -r requirements.txt
```

**Run the augmented sampler:**
```bash
python augmentedGibbsSampler.py
```

---

## References

1. Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*, Vol. 1
2. NCBI: *M. tuberculosis* genome [NC_018143.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_018143.1)
3. CollecTF Database: DosR binding site data [EXPSITE_000013c0](http://www.collectf.org/EXPSITE_000013c0)
