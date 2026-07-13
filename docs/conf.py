# Configuration file for the Sphinx documentation builder.
#
# For a full list of options see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from datetime import datetime

# -- Project information -----------------------------------------------------

project = "UKBBGym"
author = "Shubhankar Londhe"
copyright = f"{datetime.now():%Y}, {author}."
release = "0.1.0"
version = release

repository_url = "https://github.com/ShubhankarLondhe/ukbgym"

html_context = {
    "display_github": True,
    "github_user": "ShubhankarLondhe",
    "github_repo": "ukbgym",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.mathjax",
    "sphinxext.opengraph",
]

nitpicky = False
needs_sphinx = "4.0"

myst_heading_anchors = 6
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")

nb_execution_mode = "off"
nb_output_stderr = "remove"
nb_merge_streams = True

source_suffix = {
    ".rst": "restructuredtext",
    ".myst": "myst-nb",
    ".md": "myst-nb",
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

html_title = project
html_logo = "_static/img/logo.png"
html_favicon = "_static/img/favicon.png"

html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "path_to_docs": "docs/",
    "navigation_with_keys": False,
}

pygments_style = "default"
