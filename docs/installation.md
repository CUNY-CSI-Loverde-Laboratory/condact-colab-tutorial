# Installation

## Requirements

- Python 3.9 or later
- [MDAnalysis](https://www.mdanalysis.org) 2.0 or later (installed automatically)
- NumPy and Matplotlib (installed automatically with MDAnalysis)
- pandas and Jupyter are used in the tutorials but are not needed by CONDACT itself

## Install with pip

It is good practice to install into a fresh environment:

```bash
conda create -n condact python=3.11
conda activate condact
```

Then install CONDACT directly from GitHub:

```bash
pip install git+https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity.git
```

or from a local copy, if you want to read or change the code:

```bash
git clone https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity.git
cd conditional_activity
pip install .
```

## Check that it works

```bash
python -c "from conditional_activity.analysis.condact import CONDACT; print('CONDACT is installed')"
```

:::{note}
Importing CONDACT prints a line beginning `Total runtime of the program`. It is harmless and can be ignored.
:::

## No installation: Google Colab

The quick-start tutorial is also available as a notebook that runs in your browser: [open it in Google Colab](https://colab.research.google.com/github/CUNY-CSI-Loverde-Laboratory/condact-colab-tutorial/blob/main/CONDACT_Colab_Tutorial.ipynb).

## Input files

CONDACT accepts any topology and trajectory pair that MDAnalysis can read (AMBER, GROMACS, NAMD/CHARMM and others), and residues are chosen with the usual MDAnalysis `resid` selection syntax.

:::{tip}
Remove water and ions from the trajectory before analysis. CONDACT only needs the protein and nucleic acid atoms, and smaller files load much faster.
:::
