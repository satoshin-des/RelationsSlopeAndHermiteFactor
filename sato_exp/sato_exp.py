from __future__ import annotations

import numpy as np
import importlib.resources as resources

from . import core
from . import reduction
from . import svp

class Lattice:
    """
    A class for representing a lattice in n-dimensional space.
    """
    def __init__(self, basis) -> None:
        """Initialize the Lattice class.

        Args:
            basis (array like): The basis vectors of the lattice.
        """
        self.basis = basis
        self.n, self.m = basis.shape

    def __repr__(self) -> str:
        """Return a string representation of the Lattice object.

        Returns:
            str: A string representation of the Lattice object.
        """
        return f"Lattice(basis={self.basis})"

    def __str__(self) -> str:
        """Return a string representation of the Lattice object.

        Returns:
            str: A string representation of the Lattice object.
        """
        return f"{self.n}-dimensional lattice with basis:\n{self.basis}"
    
    def enum_sv(self, pruning: bool = False, alg: str = "gs") -> np.ndarray[int]:
        """Enumerates the shortest vector in the lattice basis using the SVP algorithm.
        
        ## Reference
        - N. Gama and P. Q. Nguyen and O. Regev. Lattice enumeration using extreme pruning. 2010

        Args:
            pruning (bool, optional): Whether to use pruning. Defaults to False.
            alg (str, optional): The algorithm to use ('gs' for Gram-Schmidt, 'qr' for QR decomposition). Defaults to "gs".

        Returns:
            np.ndarray[int]: The shortest vector found in the lattice.
        """
        return svp.enum_sv(self.basis, pruning, alg)

    def bkz(self, delta: float = 0.99, beta: int = 20, max_loops: int = -1, pruning: bool = False) -> Lattice:
        """Perform BKZ reduction on the lattice basis with given delta and beta parameters.
        
        ## Reference
        - C.-P. Schnorr and M. Euchner. Lattice basis reduction: Improved practical algorithms and solving subset sum problems. 1994

        Args:
            delta (float, optional): The delta parameter for BKZ reduction. Defaults to 0.99.
            beta (int, optional): The block size parameter for BKZ reduction. Defaults to 20.
            max_loops (int, optional): The maximum number of tours through the basis. Defaults to -1 (no limit).
            pruning (bool, optional): Whether to use pruning in the SVP solver. Defaults to False.
            output_sl_log (bool, optional): Whether to output the GSA-slope log. Defaults to False.
            output_rhf_log (bool, optional): Whether to output the RHF log. Defaults to False.
            output_err (bool, optional): Whether to output the error. Defaults to False.

        Returns:
            Lattice: The reduced basis.
        """
        reduced_basis = reduction.bkz(self.basis, delta, beta, max_loops, pruning)
        return Lattice(reduced_basis)

def svp_challenge(dim: int, seed: int) -> Lattice:
    """Return the basis of the SVP challenge lattice.

    Returns:
        Lattice: The basis of the SVP challenge lattice.
    """
    with resources.files("sato_exp.svp_challenge_list").joinpath(f"svp_challenge_{dim}_{seed}.txt").open('r') as f:
        basis = np.loadtxt(f, dtype=np.int64)[: dim, :dim]
    return Lattice(basis)
