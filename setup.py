from __future__ import absolute_import
from setuptools import setup, find_packages, Extension

import sys, os
from six.moves import map
sys.path.insert(0, os.path.join(os.getcwd(), 'src/'))

from sparse_pauli import __version__ as v

config = {
    'description': 'Set-based implementation of sparse Pauli operators',
    'author': 'Ben Criger',
    'url': 'https://github.com/bcriger/sparse_pauli',
    'download_url': 'https://github.com/bcriger/sparse_pauli.git',
    'author_email': 'bcriger@gmail.com',
    'version': '.'.join(map(str, v)),
    'install_requires': ['nose'],
    'package_dir': {'': 'src'},
    'packages': ['sparse_pauli'],
    'scripts': [],
    'name': 'sparse_pauli'
}

setup(**config)