import torch
from rpy2.robjects import FloatVector


def apply_r_log_pdf(target_log_pdf, x_tensor: torch.Tensor) -> torch.Tensor:
    """
    Apply an R log PDF function to a PyTorch tensor.
    
    Args:
        target_log_pdf: An R function that computes the log PDF, which takes a FloatVector as input and returns a scalar.
        x_tensor (torch.Tensor): A PyTorch tensor of shape (n_samples, n_features) for which to compute the log PDF.
    
    Returns:
        torch.Tensor: A tensor of log probabilities for each sample.
    """
    log_probs = []
    
    for row in x_tensor:
        row_np = row.detach().cpu().numpy()
        r_vec = FloatVector(row_np.tolist())
        log_prob = target_log_pdf(r_vec)[0]  # Extract R scalar
        log_probs.append(log_prob)
        
    return torch.tensor(log_probs, dtype=torch.float32, device=x_tensor.device)