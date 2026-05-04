"""
Top-level SPARC-atomSFE package.

Re-exports for convenience (only __init__.py, no stub files):
    from atom import AtomicDFTSolver
    from atom.solver import AtomicDFTSolver
    from atom.data import DataManager
    from atom.data.data_manager import ...
"""

import sys
import importlib

_SUBPACKAGES = ("solver", "data", "xc", "mesh", "pseudo", "utils", "scf")
__all__ = ["AtomicDFTSolver"] + list(_SUBPACKAGES)


def _register_submodules():
    """Register atom.src.xxx as atom.xxx in sys.modules so 'from atom.xxx import ...' works."""
    for name in _SUBPACKAGES:
        if f"atom.{name}" not in sys.modules:
            mod = importlib.import_module(f"atom.src.{name}")
            sys.modules[f"atom.{name}"] = mod


# Register on first import so "from atom.solver import X" works without stub files
_register_submodules()


def __getattr__(name: str):
    if name == "AtomicDFTSolver":
        from atom.src.solver import AtomicDFTSolver
        return AtomicDFTSolver
    if name in _SUBPACKAGES:
        return sys.modules[f"atom.{name}"]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


