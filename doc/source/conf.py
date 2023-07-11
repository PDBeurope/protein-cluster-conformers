# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import cluster_conformers

project = "protein-cluster-conformers"
copyright = "2023, Protein Data Bank in Europe"
author = "Protein Data Bank in Europe"
version = cluster_conformers.__version__
release = cluster_conformers.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx_markdown_tables",
    "sphinx.ext.coverage",
    "sphinx.ext.intersphinx",
    "myst_parser",
]


intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_parsers = {".rst": "restructuredtext", ".txt": "markdown", ".md": "markdown"}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [
    "_static",
    # "../../conf/config.yaml",
    "../../.pre-commit-config.yaml",
]

html_css_files = ["css/styles.css"]


html_theme_options = {
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
}
html_logo = "_static/logo.png"
pygments_style = "sphinx"

# region Extension configuration
github_doc_root = "https://github.com/rtfd/recommonmark/tree/master/doc/"


# def setup(app):
#     app.add_config_value(
#         "recommonmark_config",
#         {
#             "url_resolver": lambda url: github_doc_root + url,
#             "auto_toc_tree_section": "Contents",
#         },
#         True,
#     )
#     app.connect("autodoc-process-docstring", no_namedtuple_attrib_docstring)


# def no_namedtuple_attrib_docstring(app, what, name, obj, options, lines):
#     is_namedtuple_docstring = len(lines) == 1 and lines[0].startswith(
#         "Alias for field number"
#     )
#     if is_namedtuple_docstring:
#         # We don't return, so we need to purge in-place
#         del lines[:]


# endregion
