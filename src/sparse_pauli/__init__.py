__version__ = (0, 0, 0)

import pauli as _p

__modules = [_p]
map(reload, __modules)

from pauli import *

__all__ = reduce(lambda a, b: a + b, map(lambda mod: mod.__all__,
                                         __modules)) + ['__version__']
