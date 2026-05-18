from dataclasses import dataclass, field
from typing import Callable, Tuple, Optional
import pandas as pd
import numpy as np
from rpy2.robjects import FloatVector
import torch


def is_positive_definite(matrix: torch.Tensor) -> bool:
    """
    Check if a matrix is positive definite.
    
    Args:
        matrix (torch.Tensor): The matrix to check.
    
    Returns:
        bool: True if the matrix is positive definite, False otherwise.
    """
    try:
        _ = torch.linalg.cholesky(matrix)
        return True
    except RuntimeError:
        return False