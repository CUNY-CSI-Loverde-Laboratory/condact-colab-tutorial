# Frequently asked questions

## Most of my matrix is `NaN`. What went wrong?

Probably nothing in the code. A residue is only analysed if it made at least 10 transitions, and every entry involving a residue with fewer is reported as `NaN`. Check the number of transitions per residue (see {doc}`tutorials/01-quickstart`). If the counts are low, you need a longer trajectory, or a simulation in which those residues adopt multiple transitions. If they are unexpectedly low or high, check your state boundaries ({doc}`tutorials/02-state-boundaries`).

## How long does my simulation need to be?

Long enough for the residues you care about to change state many times. Buried side chains can stay in one state for hundreds of nanoseconds, so there is no single answer. Plot the angle against time for a few residues of interest and count the jumps.

## What should I use for `saving_frequency`?

The time between saved frames in the trajectory file, in picoseconds. It is not the integration time step. If you saved a frame every 5000 steps of 2 fs, the saving frequency is 10 ps.

## What do the lines "had statistically fewer transitions" mean?

They are information, not errors. CONDACT prints one for every pair of residues that it could not analyse because one of the two had too few transitions.

## Why is `A[X][Y]` different from `A[Y][X]`?

Because the measure is directional. `A[X][Y]` describes how X responds after Y changes; `A[Y][X]` describes how Y responds after X changes. One can be large while the other is close to zero.

## Is a large conditional activity enough to claim two residues are coupled?

No. Always run the significance test ({doc}`tutorials/03-significance`) and base conclusions on the q-values.

## Which files does CONDACT write?

All in the current working directory: `Dihedral_Angle.csv`, `transition_times.pkl`, `wait_times.pkl`, `Conditional_Activity_list.npy`, `Conditional_Activity_matrix.npy` and `Dynamic_Memory.npy`. The significance test adds `CA_all_tested_edges.csv` and `CA_significant_edges.csv`.

:::{warning}
These files are overwritten each time you run CONDACT in the same folder. Run each analysis in its own folder.
:::

## Can I use a degree of freedom other than χ1?

The current version uses χ1 for amino acids and the sugar-base dihedral for nucleotides. Other degrees of freedom would require changing the code.

## Where do I report a problem?

On the [issue tracker](https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity/issues) of the source repository.
