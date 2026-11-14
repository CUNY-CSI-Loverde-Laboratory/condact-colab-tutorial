# About CONDACT

## What it measures

Many methods for analysing simulations ask whether two residues tend to be in related *positions* (for example mutual information or principal component analysis). CONDACT asks a different question: whether two residues tend to *change* at related *times*.

Each residue is reduced to a single angle:

- for amino acids, the first side-chain dihedral angle, χ1 (alanine and glycine have none and are skipped);
- for nucleotides, the dihedral angle between the sugar and the base.

That angle usually sits near one of a few preferred values, called **states**. A jump from one state to another is a **transition**. CONDACT records the time of every transition of every selected residue, and all later quantities are calculated from those lists of times.

## Definitions

For two residues X and Y:

| Quantity | Meaning |
|---|---|
| Transition times | The times at which a residue changes state. |
| Waiting time | The time between two consecutive transitions of the same residue. |
| Persistence time, ⟨τ<sub>X</sub>⟩ | How long, on average, X waits before its next transition, measured from a random moment. |
| Exchange time, ⟨τ<sub>X\|Y</sub>⟩ | How long, on average, X waits before its next transition, measured from a transition of Y. |

The **conditional activity** compares the two waiting times:

$$A[X][Y] = -\ln\left(\frac{\langle\tau_{X|Y}\rangle}{\langle\tau_X\rangle}\right)$$

| Value | Interpretation |
|---|---|
| `A[X][Y] > 0` | X changes sooner than usual after Y changes. Y promotes X. |
| `A[X][Y] = 0` | The transitions of X and Y are independent. |
| `A[X][Y] < 0` | X changes later than usual after Y changes. Y suppresses X. |
| `A[X][X]` | The **dynamical memory** of X. A positive value means the transitions of X come in bursts. |

The measure is directional: `A[X][Y]` is in general not equal to `A[Y][X]`. In the matrix that CONDACT returns, the row is the responding residue X and the column is the residue Y that changes first.

## What you can do with the result

- **Dynamical memory**: the diagonal of the matrix.
- **Directional coupling**: the off-diagonal entries, tested for statistical significance with a permutation test and a false discovery rate correction ({doc}`tutorials/03-significance`).
- **Dynamically connected domains**: groups of mutually coupled residues, from the principal eigenvector of the matrix ({doc}`tutorials/04-domains-routes`).
- **Communication routes**: chains of coupled residues that are also close in space ({doc}`tutorials/04-domains-routes`).

## Limitations

:::{important}
- **Long trajectories are required.** A residue needs at least 10 transitions to be analysed. Residues with fewer are reported as `NaN`.
- **States are a coarse-graining.** Splitting a continuous angle into a few states can hide small conformational changes, and badly placed boundaries distort the result ({doc}`tutorials/02-state-boundaries`).
- **Cost grows with the square of the number of residues**, because every pair is compared.
:::

## Authors

CONDACT was written by Augustine C. Onyema, assisted by Chukwuebuka Dikeocha, in the Loverde Laboratory at the College of Staten Island and the Graduate Center, City University of New York. It is built on [MDAnalysis](https://www.mdanalysis.org) as an MDAKit and is released under the GNU General Public License, version 2.
