---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Tutorial 4. Connected domains and communication routes

Once you have a conditional activity matrix for many residues, two further analyses summarise it: which residues form a group that moves in step, and along which chain of residues a signal travels.

:::{important}
The example has only three residues, so the results here are trivial. The purpose of this tutorial is to show the commands. Both analyses become informative for a whole protein or complex.
:::

```{code-cell} ipython3
:tags: [hide-input, remove-stdout]

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances

import conditional_activity
from conditional_activity.analysis.condact import CONDACT
from conditional_activity.analysis.trace_route import CONDACT_CommPath

example_dir = os.path.join(os.path.dirname(conditional_activity.__file__), "analysis")
u = mda.Universe(os.path.join(example_dir, "pytest_LYS.prmtop"),
                 os.path.join(example_dir, "pytest_LYS.xtc"))
SELECTED_RESID = "1-3"

study = CONDACT(universe=u, selected_resid=SELECTED_RESID, saving_frequency=10)
conditional_activity_list, persistence_times, exchange_times = study.mean_conditional_activity()
CA = np.array(conditional_activity_list, dtype=float)
_, residue_names, residue_type = study.residues()
labels = [f"{name} {res.resid}" for name, res in zip(residue_names, residue_type)]
```

The cell above (click to show it) repeats {doc}`01-quickstart` and leaves the matrix in `CA` and the residue names in `labels`.

## 4.1 Dynamically connected domains

Residues that are all coupled to one another form a *dynamically connected domain*. The principal eigenvector of the matrix (the eigenvector with the largest eigenvalue) picks out the dominant one: residues with a large component belong to it, residues near zero do not.

Two preparation steps are needed. `NaN` entries are set to zero (no evidence of coupling), and the matrix is made symmetric by averaging the two directions.

```{code-cell} ipython3
CA_filled = np.nan_to_num(CA, nan=0.0)
CA_symmetric = 0.5 * (CA_filled + CA_filled.T)

eigenvalues, eigenvectors = np.linalg.eigh(CA_symmetric)
principal = np.abs(eigenvectors[:, np.argmax(eigenvalues)])

print(f"Largest eigenvalue: {eigenvalues.max():.2f}")
pd.Series(principal, index=labels, name="weight in the main domain").round(2)
```

GLN 1 and ASP 3 form the domain; ASN 2, which had too few transitions, is left out.

:::{tip}
For a large system, plot the weights against residue number or colour a structure by them. The following eigenvectors (second largest eigenvalue, and so on) can reveal further domains.
:::

## 4.2 The distance matrix

A signal passes directly only between residues that are close in space. `CONDACT_CommPath` therefore needs the average distance between every pair of residues, in the same order as the matrix.

```{code-cell} ipython3
residues_sel = u.select_atoms(f"resid {SELECTED_RESID} and not (resname ALA or resname GLY)")
dist = np.zeros((len(labels), len(labels)))
for frame in u.trajectory:
    centers = residues_sel.center_of_geometry(compound="residues")
    dist += distances.distance_array(centers, centers, box=u.dimensions)
dist /= len(u.trajectory)

pd.DataFrame(dist, index=labels, columns=labels).round(1)   # in Angstrom
```

:::{caution}
Use exactly the same selection for the distance matrix as for CONDACT, including the exclusion of alanine and glycine. Otherwise the rows of the two matrices will not refer to the same residues.
:::

## 4.3 Build the network

`CONDACT_CommPath` links Y → X when the conditional activity `A[X][Y]` is greater than `min_ca` **and** the two residues are closer than `cutoff` (in Å).

```{code-cell} ipython3
comm = CONDACT_CommPath(CA_filled, dist, labels=labels)

neighbors = comm.build_adjacency(cutoff=30.0, min_ca=0.0)
for source, targets in zip(labels, neighbors):
    for target, value in targets:
        print(f"{source} → {labels[target]}   (conditional activity {value:.2f})")
```

:::{warning}
The three example residues are far apart, so a very loose cutoff of 30 Å is used here to have something to show. For a real system use a contact distance, for example 5 to 8 Å.
:::

## 4.4 Trace routes

**List the best routes from the strongest link.** With `start_mode="global_max"` the search starts at the source of the strongest link in the network and lists routes by score (the sum of the conditional activities along the route).

```{code-cell} ipython3
paths = comm.enumerate_paths(start=None, start_mode="global_max",
                             cutoff=30.0, min_ca=0.0, max_depth=5, top_k=10)
comm.print_paths(paths, top=10)
```

**Follow the strongest link at each step** from a residue you choose, optionally stopping at a target:

```{code-cell} ipython3
paths = comm.trace_greedy(start="ASP 3", target="GLN 1",
                          cutoff=30.0, min_ca=0.0, branch_k=1)
comm.print_paths(paths)
```

| Option | Meaning |
|---|---|
| `start`, `target` | A residue label such as `"ASP 3"`, or its row number. `target` may be `None`. |
| `cutoff` | Maximum distance, in Å, for two residues to be linked. |
| `min_ca` | Minimum conditional activity for a link. |
| `top_k` | Keep only the strongest `top_k` links leaving each residue. |
| `max_depth` | Maximum number of steps in a route (`enumerate_paths`). |
| `branch_k` | `1` follows the single strongest link; larger values keep several candidate routes (`trace_greedy`). |

:::{tip}
To trace routes through statistically supported links only, set the entries of the matrix that are not significant to zero before creating `CONDACT_CommPath`, using the `significant` matrix from {doc}`03-significance`:

```python
CA_significant = np.where(significant, CA_filled, 0.0)
comm = CONDACT_CommPath(CA_significant, dist, labels=labels)
```
:::
