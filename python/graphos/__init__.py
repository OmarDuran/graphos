"""graphos: a metric-free computational topology engine.

Cells stratified by dimension, signed boundary operators, chain maps, and
an operation calculus (coproduct, product, join, quotient, pushout,
excision, cutting, subdivision, collapse, amalgamation, coarsening) with
homology, manifoldness, and amalgamation-safety queries. The formal
specification lives in THEORY.md; geometry never enters — geometric
decisions arrive as identifications, markers, and labels.
"""

from ._core import *  # noqa: F401,F403
from ._core import __doc__ as _core_doc  # noqa: F401
