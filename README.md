# two-row-sliding-dominoes

Exhaustive verification code for the paper

> L. Dong and Z. Li, *Reachability, exact diameters, and label-erasure heuristics in two-row sliding-domino puzzles* (submitted to the ICGA Journal).

The paper studies horizontal and vertical dominoes that slide without rotating on a board with two rows and an even number of vacant cells. It proves a complete reachability invariant, exact maximum distances under unit, straight and finger moves, closed-form distances for two vacancies, and the distance lost when a pattern-database-style abstraction forgets the labels of the dominoes. This repository contains one Python program that enumerates every configuration of small boards, builds the three move graphs, runs breadth-first search, and checks the closed forms of the paper against the search results.

## Requirements

Python 3.8 or later. The program uses only the standard library and contains no random numbers, so every run produces the same output.

## Quick start

```bash
python verify_two_row.py --out results/verification_log.txt
python verify_two_row.py --heuristic-only --out results/heuristic_log.txt
```

The first command runs all closed-form checks (about 30 seconds on a desktop computer) and stops with an assertion error if any check fails. The second command computes the heuristic statistics of Table 2 (about 10 seconds). The files in `results/` are the outputs reported in the paper.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--nmax-q1` | 10 | largest board width for the two-vacancy checks |
| `--nmax-gen` | 8 | largest board width for the checks with any number of vacancies |
| `--kmax-lab` | 4 | largest half-width `k` of the labeled chamber |
| `--heuristic-only` | off | compute only the label-erasure heuristic statistics |
| `--out` | `verification_log.txt` | output file |

## What is checked

| Statement in the paper | Range | Check |
|---|---|---|
| Chamber signature is a complete reachability invariant; component count (Theorem 3.4) | all inventories, width 8 or less | components of the unit-move graph match the signatures and the count formula |
| Exact diameters `D_U = QN`, `D_S = min{n-1, 2N}`, `D_F = N` (Theorems 3.7 and 3.8, Table 3) | width 8 or less, every vacancy count | maximum component diameters from breadth-first search |
| Grid-chain and triangular components, distance formulas (Theorems 4.1, 4.2, 5.1, 5.2) | two vacancies, width 10 or less | formulas compared with search on all 32,612 ordered pairs of mutually reachable states, under each move convention |
| Partial cubes (Theorem 6.1) and the four-vacancy counterexample (Example 6.2) | two vacancies, width 10 or less | Djokovic-Winkler relation; isometric dimension equals diameter |
| Component and state counts, generating functions, extremal size and diameter (Theorem 7.1, Corollary 7.2) | two vacancies, width 10 or less | exact counts per inventory |
| Reduction of the labeled chamber to a weighted 2 x k tile puzzle (Theorem 8.1) and its two components (Corollary 8.2) | k = 2, 3, 4 | unit distance equals the (2,1)-weighted tile distance and straight distance equals the tile distance on aligned states |
| Labeled diameters and a pair with adjacent unlabeled images (Section 10) | k = 2, 3, 4 | exact diameters for k = 2, 3; distances of the reflection pair |
| Label-erasure heuristic quality (Table 2) | k = 2, 3, 4 | all ordered pairs for k = 2, 3; the first 24 states of the sorted state list as sources for k = 4 |
| Mixed rectangular pieces (Example 9.3) | one 2 x 15 board | 17 states and diameter 9 under each move convention |

## Repository layout

```
verify_two_row.py          the verification program
results/verification_log.txt   output of the closed-form checks
results/heuristic_log.txt      output used for Table 2
CITATION.cff               citation metadata
LICENSE                    MIT licence
```

## Citation

If you use this code, please cite the paper. Citation metadata are in `CITATION.cff`; the entry will be updated with the journal reference after publication.

## License

MIT. See `LICENSE`.
