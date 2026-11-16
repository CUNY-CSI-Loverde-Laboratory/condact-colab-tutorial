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

# Tutorial 3. Statistical significance

A large conditional activity suggests that two residues are coupled, but two residues that both change often will sometimes change close together by chance. This tutorial shows how to tell the difference.

```{code-cell} ipython3
:tags: [hide-input, remove-stdout]

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import MDAnalysis as mda

import conditional_activity
from conditional_activity.analysis.condact import CONDACT
from conditional_activity.analysis.fdr import CONDACT_FDR

example_dir = os.path.join(os.path.dirname(conditional_activity.__file__), "analysis")
u = mda.Universe(os.path.join(example_dir, "pytest_LYS.prmtop"),
                 os.path.join(example_dir, "pytest_LYS.xtc"))
SAVING_FREQUENCY = 10
TOTAL_TIME = len(u.trajectory) * SAVING_FREQUENCY

study = CONDACT(universe=u, selected_resid="1-3", saving_frequency=SAVING_FREQUENCY)
conditional_activity_list, persistence_times, exchange_times = study.mean_conditional_activity()
```

The cell above (click to show it) repeats {doc}`01-quickstart`. The significance test reads two of the files that CONDACT wrote: `transition_times.pkl` and `Conditional_Activity_matrix.npy`.

## 3.1 The idea

If Y really influences X, the *timing* of their transitions must line up. The test breaks that timing on purpose:

1. Shift all the transition times of Y by the same random amount, wrapping around the end of the trajectory. Y keeps its own rhythm, but any relationship to X is destroyed.
2. Recompute the conditional activity.
3. Repeat many times. The collection of values shows what arises by chance alone.

The **p-value** is the fraction of shifted runs that gave a value at least as large as the real one.

Because every pair of residues is tested, some will look significant by luck. The **q-value** is the p-value corrected for this with the Benjamini–Hochberg false discovery rate (FDR) procedure.

:::{important}
Use the **q-value** to decide. A q-value below 0.05 means that, among all the pairs you call significant, about 5% are expected to be false positives.
:::

## 3.2 Run the test

```{code-cell} ipython3
fdr = CONDACT_FDR(
    transition_times_file="transition_times.pkl",
    CA_matrix_file="Conditional_Activity_matrix.npy",
    total_time=TOTAL_TIME,   # length of the simulation in ps
    n_perm=2000,             # number of random shifts
    q=0.05,                  # significance threshold for the q-value
    CA_cutoff=0.0,           # only test entries with conditional activity above this
    random_seed=0,           # makes the result reproducible
)
pvals, qvals, significant = fdr.run()
```

`run()` returns three matrices of the same shape as the conditional activity matrix: the p-values, the q-values, and `True`/`False` for significance. Entries that were not tested are `NaN`.

## 3.3 Read the results

The same information is written to `CA_all_tested_edges.csv` (every pair tested) and `CA_significant_edges.csv` (significant pairs only).

```{code-cell} ipython3
results = pd.read_csv("CA_all_tested_edges.csv")
results = results.rename(columns={"Residue_i": "X (responds)", "Residue_j": "Y (changes first)"})
results[["Y (changes first)", "X (responds)", "Conditional_Activity", "p_value", "q_value", "Significant_FDR"]].round(3)
```

- **ASP 3 → GLN 1** is significant: its q-value is below 0.05.
- **GLN 1 → ASP 3** is not: values this large turn up often by chance.

:::{note}
Only two pairs were tested. A pair is tested only if both residues have at least 11 recorded transitions (`min_transitions`) and the conditional activity is at least `CA_cutoff`. The diagonal (dynamical memory) is not tested.
:::

## 3.4 Keep only the significant couplings

After the FDR correction, only the conditional activity values that are statistically significant should be used for downstream analysis and for interpreting the results.

Keep the entries with a q-value of 0.05 or less and set all the others to `NaN`. This removes the pairs that failed the correction (q > 0.05) and the entries that were not tested, which have a `NaN` q-value. The diagonal (dynamical memory) is not tested, so it is removed as well.

```{code-cell} ipython3
CA = np.array(conditional_activity_list, dtype=float)
_, residue_names, residue_type = study.residues()
labels = [f"{name} {res.resid}" for name, res in zip(residue_names, residue_type)]

CA_filtered = np.where(qvals <= 0.05, CA, np.nan)
pd.DataFrame(CA_filtered, index=labels, columns=labels).round(2)
```

Only ASP 3 → GLN 1 remains. {doc}`04-domains-routes` uses this filtered matrix.

## 3.5 How many shifts?

The smallest p-value the test can report is 1 / (`n_perm` + 1), and the estimate becomes less noisy as `n_perm` grows. Here the test is repeated with five different random seeds, for a small and a large number of shifts.

```{code-cell} ipython3
:tags: [remove-stdout]

rows = []
for n_perm in (200, 5000):
    for seed in range(5):
        p, _, _ = CONDACT_FDR(total_time=TOTAL_TIME, n_perm=n_perm, random_seed=seed).run()
        rows.append({"n_perm": n_perm, "seed": seed, "p-value, ASP 3 → GLN 1": p[0, 2]})
spread = pd.DataFrame(rows).pivot(index="seed", columns="n_perm", values="p-value, ASP 3 → GLN 1")
```

```{code-cell} ipython3
spread.round(4)
```

With 200 shifts the p-value changes noticeably from one seed to the next; with 5000 it is much more stable.

:::{tip}
Use a few hundred or a few thousand shifts while exploring, and 10,000 or more for results you intend to publish. The time taken grows in proportion to `n_perm` and to the number of pairs tested.
:::

## 3.6 Options

| Option | Default | Meaning |
|---|---|---|
| `total_time` | read from `Dihedral_Angle.csv` | Length of the simulation in ps. Passing it explicitly is recommended. |
| `n_perm` | 200 | Number of random shifts per pair. |
| `q` | 0.05 | Threshold on the q-value. |
| `CA_cutoff` | 0.0 | Only entries with conditional activity at or above this value are tested. Raising it reduces the number of tests. |
| `min_transitions` | 11 | Minimum number of recorded transition times for a residue to be tested. |
| `random_seed` | 0 | Seed for the random shifts. |

:::{caution}
With a large selection the number of pairs grows quickly: 200 residues give almost 40,000 directed pairs. Start with a region of interest, or raise `CA_cutoff`, before testing a whole complex with many shifts.
:::
