# Sphinx configuration for the CONDACT documentation and tutorials.
project = "CONDACT"
author = "Loverde Laboratory, CUNY"
copyright = "2026, Loverde Laboratory, CUNY"

extensions = [
    "myst_nb",
    "sphinx.ext.mathjax",
]

# Tutorials are text notebooks that are run every time the site is built,
# so the numbers and figures always match the current CONDACT code.
nb_execution_mode = "force"
nb_execution_timeout = 600
nb_execution_raise_on_error = True
nb_merge_streams = True

myst_enable_extensions = ["colon_fence", "dollarmath", "amsmath"]
myst_heading_anchors = 3

exclude_patterns = ["_build", "**.ipynb_checkpoints", "requirements.txt"]

html_theme = "sphinx_rtd_theme"
html_title = "CONDACT documentation"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False}
html_context = {
    "display_github": True,
    "github_user": "CUNY-CSI-Loverde-Laboratory",
    "github_repo": "condact-colab-tutorial",
    "github_version": "main",
    "conf_py_path": "/docs/",
}
