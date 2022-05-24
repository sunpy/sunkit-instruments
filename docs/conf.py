"""
Configuration file for the Sphinx documentation builder.
"""
# -- stdlib imports ------------------------------------------------------------
import os
import datetime
from packaging.version import Version

# -- Read the Docs Specific Configuration --------------------------------------
on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if on_rtd:
    os.environ["SUNPY_CONFIGDIR"] = "/home/docs/"
    os.environ["HOME"] = "/home/docs/"
    os.environ["LANG"] = "C"
    os.environ["LC_ALL"] = "C"
    os.environ["HIDE_PARFIVE_PROGESS"] = "True"

# -- Non stdlib imports --------------------------------------------------------
from sunkit_instruments import __version__  # NOQA

# -- Project information -------------------------------------------------------
project = "sunkit_instruments"
author = "The SunPy Community"
copyright = f"{datetime.datetime.now().year}, {author}"
release = __version__
sunkit_instruments_version = Version(__version__)
is_release = not (
    sunkit_instruments_version.is_prerelease or sunkit_instruments_version.is_devrelease
)
# -- General configuration -----------------------------------------------------
suppress_warnings = [
    "app.add_directive",
]
extensions = [
    "sphinx_gallery.gen_gallery",
    "matplotlib.sphinxext.plot_directive",
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx_changelog",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
]
html_extra_path = ["robots.txt"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for intersphinx extension -----------------------------------------
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://numpy.org/doc/stable/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/numpy.inv"),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/scipy.inv"),
    ),
    "matplotlib": (
        "https://matplotlib.org/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/matplotlib.inv"),
    ),
    "sunpy": (
        "https://sunpy.org/",
        (None, "https://docs.sunpy.org/en/stable/"),
    ),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "sqlalchemy": ("https://docs.sqlalchemy.org/en/latest/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
    "drms": ("https://docs.sunpy.org/projects/drms/en/stable/", None),
    "parfive": ("https://parfive.readthedocs.io/en/stable/", None),
    "reproject": ("https://reproject.readthedocs.io/en/stable/", None),
    "aiapy": ("https://aiapy.readthedocs.io/en/stable/", None),
}
default_role = "obj"

# -- Options for HTML output ---------------------------------------------------
from sunpy_sphinx_theme.conf import *  # NOQA

# -- Graphviz Output ------------------------------------------------------------
graphviz_output_format = "svg"
graphviz_dot_args = [
    "-Nfontsize=10",
    "-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Efontsize=10",
    "-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Gfontsize=10",
    "-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
]

# -- Sphinx Gallery ------------------------------------------------------------
sphinx_gallery_conf = {
    "backreferences_dir": os.path.join("generated", "modules"),
    "filename_pattern": "^((?!skip_).)*$",
    "examples_dirs": os.path.join("..", "examples"),
    "gallery_dirs": os.path.join("generated", "gallery"),
    "matplotlib_animations": True,
    # Comes from the theme.
    "default_thumb_file": png_icon,
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("sunpy"),
    "only_warn_on_example_error": True,
}
