"""First-principles SU(3) adjoint recoupling coefficients.

The public API deliberately exposes the representation/coupling objects as
well as the final swap blocks.  This makes it possible to test the index-level
identity behind Eq. (8), instead of treating the generated numbers as an
opaque table.
"""

from .representations import Irrep, IrrepFactory, dimension
from .coupling import AdjointCouplings, Coupling, ProductDecomposition
from .dimension_only import DimensionOnlySwapTableBuilder
from .direct_specht import (
    DirectSpechtSwapOracle,
    DirectSpechtSwapTableBuilder,
    SpechtDiagnostics,
    skew_dimension,
)
from .fundamental import adjoint_split_isometry
from .fundamental_recoupling import FundamentalSplitSwapTableBuilder
from .recoupling import Path, SwapBlock, SwapTableBuilder
from .recursive_reduction import RecursiveReductionSwapTableBuilder
from .table import read_table, write_table
from .validation import assert_braid_relations, braid_residuals

__all__ = [
    "AdjointCouplings",
    "Coupling",
    "DimensionOnlySwapTableBuilder",
    "DirectSpechtSwapOracle",
    "DirectSpechtSwapTableBuilder",
    "FundamentalSplitSwapTableBuilder",
    "Irrep",
    "IrrepFactory",
    "Path",
    "ProductDecomposition",
    "RecursiveReductionSwapTableBuilder",
    "SwapBlock",
    "SwapTableBuilder",
    "SpechtDiagnostics",
    "adjoint_split_isometry",
    "dimension",
    "assert_braid_relations",
    "braid_residuals",
    "read_table",
    "skew_dimension",
    "write_table",
]
