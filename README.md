# two-row-sliding-dominoes

Exhaustive verification code for the paper

> L. Dong and Z. Li, *Exact component diameters and weighted distances in two-row sliding-domino puzzles* (submitted to Discrete Applied Mathematics).

The paper studies horizontal and vertical dominoes that slide without rotating on a board with two rows and an even number of vacant cells. It proves a complete reachability invariant, the exact unit, straight and finger diameters of every component, closed-form distances for two vacancies, weighted diameters of labeled chambers, and exact metrics for rectangular pieces on boards of arbitrary height when the vacancy counts fix their rows. This repository contains one Python program that enumerates every configuration of small boards, builds the three move graphs, runs breadth-first search, and checks the closed forms of the paper against the search results.

## Requirements

Python 3.8 or later. The program uses only the standard library and contains no random numbers, so every run produces the same output.

## Quick start

```bash
python verify_two_row.py --out results/verification_log.txt
```

This command runs all checks of Section 10 of the paper (about two and a half minutes on a desktop computer) and stops with an assertion error if any check fails. The file `results/verification_log.txt` is the output reported in the paper.

The option `--heuristic-only` computes supplementary statistics of the label-erasure heuristic in the labeled chamber (`results/heuristic_log.txt`); these statistics are not part of the paper.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--nmax-q1` | 10 | largest board width for the two-vacancy checks |
| `--nmax-gen` | 8 | largest board width for the signature, count and extremum checks |
| `--nmax-comp` | 10 | largest board width for the diameter of every individual component |
| `--kmax-lab` | 4 | largest half-width `k` of the labeled chamber |
| `--heuristic-only` | off | compute only the supplementary heuristic statistics |
| `--out` | `verification_log.txt` | output file |

## What is checked

| Statement in the paper | Range | Check |
|---|---|---|
| Unit, straight and finger diameters of every component (Theorem 3.8, Example 3.10) | all inventories, width 10 or less (1460 components) | diameters from breadth-first search equal `EN+A` and the stated macro diameters |
| Components without an extra vacancy pair are partial cubes of isometric dimension `A` (Corollary 6.2, Example 6.3) | width 10 or less (545 components) | Djokovic-Winkler relation |
| Chamber signature is a complete reachability invariant; component count (Theorem 3.4); extrema over inventories (Corollary 3.9) | width 8 or less, every vacancy count | components of the unit-move graph match the signatures, the count formula and the maximum diameters |
| Grid-chain and triangular components, distance formulas (Theorems 4.1, 4.2, 5.1, 5.2) | two vacancies, width 10 or less | formulas compared with search on all 32,612 ordered pairs of mutually reachable states, under each move convention |
| Partial cubes with two vacancies (Theorem 6.1) and the four-vacancy counterexample (Example 6.4) | two vacancies, width 10 or less | Djokovic-Winkler relation; isometric dimension equals diameter |
| Component and state counts, generating functions, extremal size and diameter (Theorem 7.1, Corollary 7.2) | two vacancies, width 10 or less | exact counts per inventory |
| Reduction of the labeled chamber to a weighted 2 x k tile puzzle (Theorem 8.1) and its two components (Corollary 8.2) | k = 2, 3, 4 | unit distance equals the (2,1)-weighted tile distance and straight distance equals the tile distance on aligned states |
| Labeled diameters (Corollary 8.6) | k = 2, 3, 4 | exact diameters for k = 2, 3; distances of a reflection pair for k = 4 |
| Rectangles on a three-row board (Theorem 9.2, Example 9.3) | one 3 x 15 board | 42 states and diameter 12 under each move convention |

## Repository layout

```
verify_two_row.py              the verification program
results/verification_log.txt   output of the checks reported in the paper
results/heuristic_log.txt      supplementary heuristic statistics
CITATION.cff                   citation metadata
LICENSE                        MIT licence
```

## Citation

If you use this code, please cite the paper. Citation metadata are in `CITATION.cff`; the entry will be updated with the journal reference after publication.

## License

MIT. See `LICENSE`.