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

# Tutorial 2. Choosing the state boundaries

CONDACT turns a continuous angle into a small number of states. Where you place the boundaries between states decides what counts as a transition, so it is the most important choice you make. This tutorial shows how to make it and what happens when it goes wrong.

```{code-cell} ipython3
:tags: [hide-input, remove-stdout]

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
u = mda.Universe(os.path.join(example_dir, "pytest_LYS.prmtop"),
                 os.path.join(example_dir, "pytest_LYS.xtc"))
SELECTED_RESID = "1-3"
SAVING_FREQUENCY = 10

def chi1_atoms(residue):
    atoms = residue.atoms
    return (atoms.select_atoms("name N") + atoms.select_atoms("name CA") + atoms.select_atoms("name CB")
            + atoms.select_atoms("name CG or name OG1 or name CG1 or name OG or name SG"))

selection = u.select_atoms(f"resid {SELECTED_RESID} and not (resname ALA or resname GLY)")
groups = {f"{res.resname} {res.resid}": chi1_atoms(res) for res in selection.residues}
angles = {name: [] for name in groups}
for frame in u.trajectory:
    for name, atom_group in groups.items():
        angles[name].append(atom_group.dihedral.value() % 360)
angles = pd.DataFrame(angles)
angles.index = angles.index * SAVING_FREQUENCY / 1000
angles.index.name = "time (ns)"
```

The cell above (click to show it) loads the example and measures χ1 for each residue in every frame, exactly as in {doc}`01-quickstart`.

## 2.1 What a transition looks like

Plot the angle of each residue against time. The angle hovers around one value, then jumps to another. Every crossing of a dashed line is one transition.

```{code-cell} ipython3
boundaries = [120, 240]

fig, axes = plt.subplots(len(angles.columns), 1, figsize=(9, 6), sharex=True)
for ax, name in zip(axes, angles.columns):
    ax.plot(angles.index, angles[name], ".", markersize=2)
    for b in boundaries:
        ax.axhline(b, color="gray", linestyle="--", linewidth=1)
    ax.set_ylim(0, 360)
    ax.set_yticks([0, 120, 240, 360])
    ax.set_ylabel(f"{name}\nχ1 (°)")
axes[-1].set_xlabel("Time (ns)")
plt.tight_layout()
```

GLN 1 jumps often, ASP 3 occasionally, and ASN 2 almost never. This already tells you which residues can be analysed: a residue needs at least 10 transitions.

## 2.2 Find the peaks and valleys

A histogram of the angles shows the preferred values. Each peak is a state, and each valley between two peaks is a boundary.

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7, 3.6))
bins = np.linspace(0, 360, 73)
for name in angles.columns:
    ax.hist(angles[name], bins=bins, density=True, histtype="step", linewidth=1.8, label=name)
for b in boundaries:
    ax.axvline(b, color="gray", linestyle="--", linewidth=1)
ax.set_xlabel("χ1 (degrees)")
ax.set_ylabel("Probability density")
ax.set_xlim(0, 360)
ax.legend(frameon=False)
plt.tight_layout()
```

Three peaks, near 60°, 190° and 300°, and two wide, empty valleys around 120° and 240°. This gives the settings

```python
No_of_peaks_protein=3, peak_boundaries_protein=[120, 240], states_protein="XYZ"
```

which are also CONDACT's defaults for proteins.

:::{important}
The number of boundaries must be one less than the number of states, and `states_protein` needs one letter per state. Otherwise CONDACT raises a `ValueError`.
:::

## 2.3 What happens if the boundaries are wrong

Run CONDACT with three different pairs of boundaries and compare.

```{code-cell} ipython3
:tags: [remove-stdout]

rows = []
for test in ([120, 240], [130, 230], [100, 260]):
    study = CONDACT(universe=u, selected_resid=SELECTED_RESID,
                    peak_boundaries_protein=test, saving_frequency=SAVING_FREQUENCY)
    ca, _, _ = study.mean_conditional_activity()
    ca = np.array(ca, dtype=float)
    rows.append({
        "boundaries": str(test),
        "transitions of GLN 1": len(study.transition_times["GLN 1"]),
        "transitions of ASP 3": len(study.transition_times["ASP 3"]),
        "A[GLN 1][ASP 3]": ca[0, 2],
        "A[ASP 3][GLN 1]": ca[2, 0],
    })
comparison = pd.DataFrame(rows).set_index("boundaries").round(2)
```

```{code-cell} ipython3
comparison
```

- Moving the boundaries to `[130, 230]` changes almost nothing, because both still lie in empty valleys.
- Moving them to `[100, 260]` nearly doubles the number of transitions of GLN 1 and cuts the main coupling to less than a third. The boundary at 260° slices into the edge of the peak near 300°, so ordinary wobbling inside one state is counted as a series of transitions.

:::{warning}
A boundary that cuts through a peak creates false transitions and changes the result. Always choose boundaries from the histogram of *your* system. Do not rely on the defaults without checking.
:::

:::{tip}
If different residues have valleys in clearly different places, analyse them in separate runs with their own boundaries. CONDACT applies one set of boundaries to all the amino acids in a run.
:::

## 2.4 DNA and RNA

For nucleotides CONDACT uses the dihedral angle between the sugar and the base, defined by the atoms O4′–C1′–N9–C4 for purines and O4′–C1′–N1–C2 for pyrimidines. The settings have the same meaning as for proteins:

```python
No_of_peaks_nucleic=3, peak_boundaries_nucleic=[140, 330], states_nucleic="SAS"
```

With these defaults, angles between 140° and 330° are the *anti* conformation (state `A`), and the two ranges on either side are both the *syn* conformation (state `S`). Using the same letter twice is how you tell CONDACT that two ranges are one state.

The example data contains no nucleotides, so nothing is run here. For your own system, make the same two plots for the sugar-base dihedral before trusting the defaults.

:::{note}
Nucleotides are recognised by residue name. CONDACT knows the AMBER-style names `DA`, `DG`, `DC`, `DT` for DNA and `A`, `G`, `C`, `U` for RNA, each optionally followed by `5` or `3` for the chain ends. A nucleotide with any other name is treated as an amino acid and will give `NaN`.
:::

## 2.5 Negative angles

MDAnalysis reports dihedral angles from −180° to 180°. By default (`keep_negative=False`) CONDACT adds 360° to negative values, so that all angles run from 0° to 360° and the boundaries above apply. Leave this setting alone unless you have a specific reason to change it.
