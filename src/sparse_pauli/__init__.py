from __future__ import absolute_import
from six.moves import map
from functools import reduce

__version__ = (0, 0, 0)

from . import pauli as _p

from .pauli import *

__all__ = _p.__all__ + ['__version__']
