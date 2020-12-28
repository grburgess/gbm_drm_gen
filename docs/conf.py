# -*- coding: utf-8 -*-
#
# gbm_drm_gen documentation build configuration file, created by
# sphinx-quickstart on Sun Oct  8 13:17:28 2017.
#
#



# General information about the project.
project = u'gbm_drm_gen'
copyright = u'2017-2020, J. Michael Burgess'
author = u'J. Michael Burgess'


import os
import sys
sys.path.insert(0, os.path.abspath('.'))



# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['nbsphinx',
              'sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.autodoc',
              'sphinx.ext.githubpages',
              'sphinx.ext.napoleon'

napoleon_google_docstring = True
napoleon_use_param = False


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','**.ipynb_checkpoints']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'

# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ['_static']


# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'css/custom.css',
]


html_js_files = [
    'css/custom.js',
]




html_show_sourcelink = False
html_favicon = "media/favicon.ico"

html_show_sphinx = False


# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

autosectionlabel_prefix_document = True

# avoid time-out when running the doc
nbsphinx_timeout = 30 * 60

nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 96}",
]

# autodoc_member_order = 'bysource'

# autoclass_content = 'both'


# edit_on_github_project = 'JohannesBuchner/UltraNest'
# edit_on_github_branch = 'master'
# #edit_on_github_url
# edit_on_github_src = 'docs/'  # optional. default: ''



html_theme_options = {
 #   'canonical_url': 'https://johannesbuchner.github.io/UltraNest/',
    'style_external_links': True,
    # 'vcs_pageview_mode': 'edit',
    'style_nav_header_background': '#470E80',
    #'only_logo': True,
}


html_logo ='media/logo.png'

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'gbm_drm_gendoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'gbm_drm_gen.tex', u'gbm\\_drm\\_gen Documentation',
     u'J. Michael Burgess', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'gbm_drm_gen', u'gbm_drm_gen Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'gbm_drm_gen', u'gbm_drm_gen Documentation',
     author, 'gbm_drm_gen', 'One line description of project.',
     'Miscellaneous'),
]



