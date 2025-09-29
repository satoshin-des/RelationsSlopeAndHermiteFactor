from .loader import lib
import numpy as np
import pandas as pd
import ctypes
import os

def bkz(basis: np.ndarray[int], delta: float = 0.99, beta: int = 20, max_loops: int = -1, pruning: bool = False) -> np.ndarray[int]:
    """Performs BKZ reduction on a basis with given delta and beta parameters.

    Args:
        basis (np.ndarray[int]): The input basis vectors.
        delta (float, optional): The delta parameter for BKZ reduction. Defaults to 0.99.
        beta (int, optional): The block size parameter for BKZ reduction. Defaults to 20.
        max_loops (int, optional): The maximum number of tours through the basis. Defaults to -1 (no limit).
        pruning (bool, optional): Whether to use pruning in the SVP solver. Defaults to False.
        output_sl_log (bool, optional): Whether to output the GSA-slope log. Defaults to False.
        output_rhf_log (bool, optional): Whether to output the RHF log. Defaults to False.
        output_err (bool, optional): Whether to output the error log. Defaults to False.

    Returns:
        np.ndarray[int]: The BKZ reduced basis.
    """
    if delta <= 0.25 or delta >= 1:
        raise ValueError("Delta must be in the range (0.25, 1).")
    if beta < 2:
        raise ValueError("Beta must be at least 2.")
    if max_loops < -1:
        raise ValueError("max_loops must be -1 (no limit) or a non-negative integer.")
    
    n, m = basis.shape

    lib.BKZ.argtypes = (
        ctypes.POINTER(ctypes.POINTER(ctypes.c_long)),  # basis
        ctypes.c_double,
        ctypes.c_long,
        ctypes.c_long,
        ctypes.c_bool,
        ctypes.c_long,
        ctypes.c_long
    )
    lib.BKZ.restype = None

    basis_ptr = (ctypes.POINTER(ctypes.c_long) * n)()
    for i in range(n):
        basis_ptr[i] = (ctypes.c_long * m)()
        for j in range(m):
            basis_ptr[i][j] = ctypes.c_long(basis[i, j])

    lib.BKZ(basis_ptr, ctypes.c_double(delta), ctypes.c_long(beta), ctypes.c_long(max_loops), ctypes.c_bool(pruning), n, m)

    reduced_basis = np.zeros((n, m), dtype=np.int64)
    for i in range(n):
        for j in range(m):
            reduced_basis[i, j] = basis_ptr[i][j]

    return reduced_basis
