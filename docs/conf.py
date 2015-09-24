import sys
import os
import shlex
from recommonmark.parser import CommonMarkParser

extensions = [
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']

source_suffix = ['.rst', '.md']

master_doc = 'index'

project = u'Kiva'
copyright = u'2015, Big Ladder Software'
author = u'Neal Kruis'

version = '0.2.1'
release = '0.2.1'
exclude_patterns = ['_build']
pygments_style = 'sphinx'

language = 'en'

todo_include_todos = False


html_theme = 'pydoctheme'

html_static_path = ['_static']
htmlhelp_basename = 'Kivadoc'
html_split_index = True
#html_theme_options = {'collapsiblesidebar': True}

latex_elements = {}

latex_documents = [
  (master_doc, 'Kiva.tex', u'Kiva Documentation',
   u'Neal Kruis', 'manual'),
]
man_pages = [
    (master_doc, 'kiva', u'Kiva Documentation',
     [author], 1)
]
