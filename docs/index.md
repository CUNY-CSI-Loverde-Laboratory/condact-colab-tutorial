# CONDACT

**CONDitional ACTivity, a kinetic framework for mapping allosteric communication**

CONDACT is an open-source Python library from the Loverde Laboratory (CUNY). It reads a molecular dynamics (MD) trajectory and asks a timing question about every pair of residues:

> When residue Y changes its conformation, does residue X tend to change soon afterwards?

Pairs for which the answer is yes are *kinetically coupled*. From these couplings CONDACT identifies residues with dynamical memory, directional communication between residues, groups of residues that move in step, and routes along which a signal can travel. It works for proteins, DNA and RNA.

:::{note}
New to CONDACT? Start with {doc}`tutorials/01-quickstart`. You can also run that tutorial in your browser, with nothing to install: [open it in Google Colab](https://colab.research.google.com/github/CUNY-CSI-Loverde-Laboratory/condact-colab-tutorial/blob/main/CONDACT_Colab_Tutorial.ipynb).
:::

Questions and bug reports are welcome on the [issue tracker](https://github.com/CUNY-CSI-Loverde-Laboratory/conditional_activity/issues).

```{toctree}
:caption: Manual
:maxdepth: 2

about
installation
citations
faq
api
```

```{toctree}
:caption: Tutorials
:maxdepth: 2

tutorials/setup
tutorials/01-quickstart
tutorials/02-state-boundaries
tutorials/03-significance
tutorials/04-domains-routes
```
