# Reference

A summary of the three classes and their options. The source code is in the [`analysis` folder](https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity/tree/main/conditional_activity/analysis) of the repository.

## `CONDACT`

```python
from conditional_activity.analysis.condact import CONDACT
```

Computes transition times, persistence and exchange times, and the conditional activity matrix. Used in {doc}`tutorials/01-quickstart`.

| Option | Default | Meaning |
|---|---|---|
| `universe` | required | An MDAnalysis `Universe`. |
| `selected_resid` | required | Residues to analyse, in MDAnalysis `resid` syntax, for example `"1-50"` or `"5 9 12"`. Alanine and glycine are skipped. |
| `No_of_peaks_protein` | `3` | Number of χ1 states for amino acids. |
| `peak_boundaries_protein` | `[120, 240]` | Boundaries between those states, in degrees. One fewer than the number of states. |
| `states_protein` | `"XYZ"` | One letter per state. |
| `No_of_peaks_nucleic` | `3` | Number of sugar-base dihedral states for nucleotides. |
| `peak_boundaries_nucleic` | `[140, 330]` | Boundaries between those states, in degrees. |
| `states_nucleic` | `"SAS"` | One letter per state. A repeated letter joins two ranges into one state. |
| `saving_frequency` | required | Time between saved frames, in picoseconds. |
| `keep_negative` | `False` | If `False`, 360° is added to negative angles so that all angles run from 0° to 360°. |

| Method | Returns |
|---|---|
| `mean_conditional_activity()` | Runs the whole analysis. Returns the conditional activity matrix (list of lists), a dictionary of persistence times, and the exchange-time matrix. Writes the output files. |
| `run()` | Reads the trajectory and returns the dictionaries of transition times and waiting times. Called by `mean_conditional_activity()`. |
| `residues()` | Returns the selected atoms, the residue names, and the residues, in the order used for the matrix. |

After a run, `transition_times` and `wait_times` are available as attributes (dictionaries keyed by labels such as `"GLN 1"`).

## `CONDACT_FDR`

```python
from conditional_activity.analysis.fdr import CONDACT_FDR
```

Permutation test with Benjamini–Hochberg FDR correction. Used in {doc}`tutorials/03-significance`.

| Option | Default | Meaning |
|---|---|---|
| `transition_times_file` | `"transition_times.pkl"` | File written by `CONDACT`. |
| `CA_matrix_file` | `"Conditional_Activity_matrix.npy"` | File written by `CONDACT`. |
| `dihedral_csv` | `"Dihedral_Angle.csv"` | Only used to work out the simulation length when `total_time` is not given. |
| `total_time` | `None` | Length of the simulation in picoseconds. |
| `n_perm` | `200` | Number of random shifts per pair. |
| `q` | `0.05` | Threshold on the q-value. |
| `CA_cutoff` | `0.0` | Only entries at or above this value are tested. |
| `min_transitions` | `11` | Minimum number of recorded transition times for a residue to be tested. |
| `random_seed` | `0` | Seed for the random shifts. |

`run()` returns the matrices of p-values, q-values and significance, and writes `CA_all_tested_edges.csv` and `CA_significant_edges.csv`.

## `CONDACT_CommPath`

```python
from conditional_activity.analysis.trace_route import CONDACT_CommPath
```

Traces directed communication routes. Used in {doc}`tutorials/04-domains-routes`.

```python
comm = CONDACT_CommPath(CA, dist, labels=labels)
```

`CA` is the conditional activity matrix with `NaN` replaced by numbers, `dist` is the matrix of average distances between residues in Å, and `labels` is the list of residue labels. An entry `CA[X, Y]` is read as a link Y → X.

| Method | Purpose |
|---|---|
| `build_adjacency(cutoff, min_ca, top_k)` | For each residue, the list of residues it links to, strongest first. |
| `enumerate_paths(start, target, cutoff, min_ca, max_depth, top_k, start_mode)` | All routes without repeated residues, up to `max_depth` steps, sorted by score. |
| `trace_greedy(start, target, cutoff, min_ca, top_k, branch_k, max_steps)` | Follows the strongest link at each step (`branch_k=1`) or keeps several candidates (`branch_k>1`). |
| `print_paths(paths, top)` | Prints routes in readable form. |
