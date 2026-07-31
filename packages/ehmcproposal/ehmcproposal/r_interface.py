import math
import torch

from rpy2.robjects import FloatVector
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.rinterface_lib import callbacks

def apply_r_log_pdf(
    target_log_pdf, 
    x_tensor: torch.Tensor,
) -> torch.Tensor:
    """
    Apply an R log PDF function to a PyTorch tensor.
    
    Args:
        target_log_pdf: An R function that computes the log PDF, which takes a FloatVector as input and returns a scalar.
        x_tensor (torch.Tensor): A PyTorch tensor of shape (n_samples, n_features) for which to compute the log PDF.
    
    Returns:
        torch.Tensor: A tensor of log probabilities for each sample.
    """
    dtype = torch.float32
    finfo = torch.finfo(dtype)

    log_probs = torch.empty(
        x_tensor.shape[0],
        dtype=dtype,
        device=x_tensor.device,
    )

    old_callback = callbacks.consolewrite_warnerror
    callbacks.consolewrite_warnerror = lambda _: None

    try:
        for i, row in enumerate(x_tensor):
            r_vec = FloatVector(row.detach().cpu().tolist())

            try:
                value = float(target_log_pdf(r_vec)[0])

                if math.isfinite(value):
                    value = max(min(value, finfo.max), finfo.min)
                else:
                    value = float("-inf")

            except RRuntimeError:
                value = float("-inf")

            log_probs[i] = value

    finally:
        callbacks.consolewrite_warnerror = old_callback

    return log_probs