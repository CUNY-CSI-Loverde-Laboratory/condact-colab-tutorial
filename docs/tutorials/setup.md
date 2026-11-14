# Before you start

## What you need

- A working installation of CONDACT (see {doc}`../installation`), or a Google account if you prefer to run the [Colab notebook](https://colab.research.google.com/github/CUNY-CSI-Loverde-Laboratory/condact-colab-tutorial/blob/main/CONDACT_Colab_Tutorial.ipynb).
- For the tutorials: `pandas`, `matplotlib` and Jupyter (`pip install pandas matplotlib jupyterlab`).

## What you should already know

The tutorials assume that you know what an MD trajectory is and that you can run Python code in a Jupyter notebook. You do not need previous experience with MDAnalysis. If you are new to it, the [MDAnalysis user guide](https://userguide.mdanalysis.org) is a good place to start.

## The example data

Every tutorial uses the small example that is installed with CONDACT, so there is nothing to download:

| | |
|---|---|
| System | three amino-acid residues (GLN 1, ASN 2, ASP 3) taken from a larger simulation |
| Length | 30 ns, one frame every 10 ps (3000 frames) |
| Files | `pytest_LYS.prmtop` (topology) and `pytest_LYS.xtc` (trajectory) |

The tutorials locate these files for you with:

```python
import os, conditional_activity
example_dir = os.path.join(os.path.dirname(conditional_activity.__file__), "analysis")
```

:::{important}
The example is deliberately tiny so that every tutorial runs in seconds. It is enough to learn the workflow, but three residues cannot show you a real allosteric network. Trajectories of the nucleosome from the CONDACT paper are available on [Zenodo](https://zenodo.org/records/17193099).
:::

## How the tutorials are organised

| Tutorial | What you learn |
|---|---|
| {doc}`01-quickstart` | The whole workflow, from a trajectory to a conditional activity matrix |
| {doc}`02-state-boundaries` | How to choose the state boundaries, and what happens when they are wrong |
| {doc}`03-significance` | How to decide which couplings are statistically significant |
| {doc}`04-domains-routes` | Dynamically connected domains and communication routes |

Work through them in order the first time. Each one can also be run on its own.

:::{note}
The numbers and figures on these pages are produced by running the code shown, each time the site is built. If you run the same code you should get the same results.
:::
