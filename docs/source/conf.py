import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

project = 'pynlcd'
copyright = '2024, Saswata Nandi'
author = 'Saswata Nandi'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

html_theme = 'sphinx_rtd_theme'

autodoc_mock_imports = ['gdal', 'ogr', 'osr']