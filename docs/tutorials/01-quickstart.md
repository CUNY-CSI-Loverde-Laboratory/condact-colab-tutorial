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

# Tutorial 1. Quick start

This tutorial takes you from a trajectory to a conditional activity matrix in four steps.

:::{note}
You can run this tutorial in your browser, with nothing to install: [open it in Google Colab](https://colab.research.google.com/github/CUNY-CSI-Loverde-Laboratory/condact-colab-tutorial/blob/main/CONDACT_Colab_Tutorial.ipynb).
:::

## 1.1 Load the trajectory

```{code-cell} ipython3
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

import conditional_activity
from conditional_activity.analysis.condact import CONDACT

example_dir = os.path.join(os.path.dirname(conditional_activity.__file__), "analysis")
topology_file = os.path.join(example_dir, "pytest_LYS.prmtop")
trajectory_file = os.path.join(example_dir, "pytest_LYS.xtc")

SELECTED_RESID = "1-3"   # which residues to analyse, e.g. "1-3" or "5 9 12"
SAVING_FREQUENCY = 10    # time between saved frames, in picoseconds

u = mda.Universe(topology_file, trajectory_file)
TOTAL_TIME = len(u.trajectory) * SAVING_FREQUENCY

print(f"{u.atoms.n_atoms} atoms, {len(u.trajectory)} frames, {TOTAL_TIME/1000:.0f} ns")
```

:::{important}
`SAVING_FREQUENCY` is the time between frames written to the trajectory file. It is not the integration time step of the simulation.
:::

## 1.2 Look at the angles

CONDACT describes each amino acid by its first side-chain dihedral angle, χ1. Before running the analysis, look at how that angle is distributed, so that you know where one state ends and the next begins.

```{code-cell} ipython3
def chi1_atoms(residue):
    # The four atoms that define chi1, in order
    atoms = residue.atoms
    return (atoms.select_atoms("name N") + atoms.select_atoms("name CA") + atoms.select_atoms("name CB")
            + atoms.select_atoms("name CG or name OG1 or name CG1 or name OG or name SG"))

selection = u.select_atoms(f"resid {SELECTED_RESID} and not (resname ALA or resname GLY)")
groups = {f"{res.resname} {res.resid}": chi1_atoms(res) for res in selection.residues}

angles = {name: [] for name in groups}
for frame in u.trajectory:
    for name, atom_group in groups.items():
        angles[name].append(atom_group.dihedral.value() % 360)   # 0 to 360 degrees

angles = pd.DataFrame(angles)
angles.index = angles.index * SAVING_FREQUENCY / 1000
angles.index.name = "time (ns)"

fig, ax = plt.subplots(figsize=(7, 3.6))
bins = np.linspace(0, 360, 73)
for name in angles.columns:
    ax.hist(angles[name], bins=bins, density=True, histtype="step", linewidth=1.8, label=name)
for boundary in (120, 240):
    ax.axvline(boundary, color="gray", linestyle="--", linewidth=1)
ax.set_xlabel("χ1 (degrees)")
ax.set_ylabel("Probability density")
ax.set_xlim(0, 360)
ax.legend(frameon=False)
plt.tight_layout()
```

There are three peaks, separated by empty valleys near 120° and 240° (dashed lines). So there are **three states** with **boundaries at 120° and 240°**. {doc}`02-state-boundaries` looks at this choice in more detail.

## 1.3 Run CONDACT

```{code-cell} ipython3
study = CONDACT(
    universe=u,
    selected_resid=SELECTED_RESID,
    No_of_peaks_protein=3,               # number of states
    peak_boundaries_protein=[120, 240],  # boundaries between them, in degrees
    states_protein="XYZ",                # one letter per state
    saving_frequency=SAVING_FREQUENCY,
)

conditional_activity_list, persistence_times, exchange_times = study.mean_conditional_activity()

CA = np.array(conditional_activity_list, dtype=float)
_, residue_names, residue_type = study.residues()
labels = [f"{name} {res.resid}" for name, res in zip(residue_names, residue_type)]
```

:::{note}
The lines saying that a residue "had statistically fewer transitions" are information, not errors. A residue must change state at least 10 times to be analysed. Here ASN 2 changed only 3 times, so every result that involves it is `NaN`. This is why real analyses need long simulations or a simulation with residues adopting multiple transitions.
:::

`mean_conditional_activity()` returns three things: the conditional activity matrix (as a list of lists), a dictionary of persistence times, and the matrix of exchange times. It also writes its results to the current folder:

| File | Contents |
|---|---|
| `Dihedral_Angle.csv` | the angle of every selected residue in every frame |
| `transition_times.pkl` | the times at which each residue changed state |
| `wait_times.pkl` | the waiting times between transitions |
| `Conditional_Activity_matrix.npy` | the conditional activity matrix |
| `Dynamic_Memory.npy` | the diagonal of the matrix |

## 1.4 Read the results

One line per residue:

```{code-cell} ipython3
summary = pd.DataFrame({
    "number of transitions": [len(study.transition_times[l]) for l in labels],
    "persistence time (ps)": [persistence_times[l] for l in labels],
    "dynamical memory A[X][X]": np.diag(CA),
}, index=labels)
summary.loc[summary["number of transitions"] < 10, "persistence time (ps)"] = np.nan
summary.round(2)
```

The dynamical memory of GLN 1 and ASP 3 is positive: their transitions come in bursts separated by quiet periods.

The full matrix:

```{code-cell} ipython3
display(pd.DataFrame(CA, index=labels, columns=labels).round(2))

fig, ax = plt.subplots(figsize=(4.8, 4))
image = ax.imshow(CA, cmap="hot_r", vmin=0, origin="upper")
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha="right")
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels)
ax.set_xlabel("Y: the residue that changes first")
ax.set_ylabel("X: the residue that responds")
fig.colorbar(image, label="Conditional activity A[X][Y]")
plt.tight_layout()
```

Each cell answers the question: after the residue in the **column** changes state, does the residue in the **row** change sooner than usual?

- Row GLN 1, column ASP 3 is large and positive: after ASP 3 changes, GLN 1 changes much sooner than it normally would. The signal runs **ASP 3 → GLN 1**.
- Row ASP 3, column GLN 1 is close to zero: a change in GLN 1 has little effect on ASP 3.
- White cells are `NaN`.

:::{caution}
A large value is suggestive, but it can arise by chance. Test it before drawing conclusions: see {doc}`03-significance`.
:::
