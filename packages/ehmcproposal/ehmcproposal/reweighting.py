import torch


def effective_sample_size(weights: torch.Tensor) -> float:
    """
    Compute the effective sample size (ESS) for a set of normalized weights.

    Args:
        weights (torch.Tensor): The normalized importance weights.

    Returns:
        float: The effective sample size.
    """
    return 1.0 / torch.sum(weights ** 2).item()
   



 
def normalize_log_weights(log_w: torch.Tensor) -> torch.Tensor:
    """
    Normalize log weights to obtain normalized importance weights.

    Args:
        log_w (torch.Tensor): Logarithm of the importance weights.

    Returns:
        torch.Tensor: Normalized importance weights.
    """
    return torch.exp(log_w - torch.logsumexp(log_w, dim=0))
   



 
def f_ess(log_weights_shift: torch.Tensor, 
          beta: float, 
          ess_t: float) -> float:
    """
    Target function for the dichotomy search to find the beta that achieves the target effective sample size (ESS).
    
    Args:
        log_weights_shift (torch.Tensor): The shifted log weights.
        beta (float): A scaling factor for the weights.
        ess_t (float): The target effective sample size.
        
    Returns:
        float: The difference between the computed ESS and the target ESS.
    """
    scaled = beta * log_weights_shift
    w = torch.exp(scaled)
    ess = (w.sum() ** 2) / (w ** 2).sum()
    return ess.item() - ess_t
   



 
def reweight_sample(
    log_weights: torch.Tensor, 
    ess_t: float, 
    beta_init: float=0.0, 
    precision: float=1e-6
) -> tuple[torch.Tensor, float, float, float]:
    """
    Reweights a sample based on the provided weights.

    Args:
        log_weights (torch.Tensor): Logarithm of the importance weights.
        ess_t (float): Target effective sample size.
        beta_init (float, optional): Initial beta value for scaling. Default is 1.0.
        precision (float, optional): Precision for binary search. Default is 1e-6.

    Returns:
        tuple:
            - scaled_weights (torch.Tensor): Normalized importance weights for the tempered target distribution.
            - scaled_ess (float): Effective sample size after reweighting.
            - beta_output (float): Final beta used after reweighting.
            - current_ess (float): Initial ESS before any reweighting.
    """
    # --- Compute the weights
    log_w_shift = log_weights - torch.logsumexp(log_weights, dim=0)
    weights = torch.exp(log_w_shift)
    weights = weights / weights.sum()
    # --- Compute current ESS
    initial_ess = effective_sample_size(weights)
    beta = beta_init
    if initial_ess < ess_t:
        a, b = beta, 1.0
        f_a = f_ess(log_w_shift, a, ess_t)
        f_b = f_ess(log_w_shift, b, ess_t)
        if (f_a * f_b > 0):
            if f_a > 0:
                beta = 1.0
            else:
                print(f'Current ESS is: {initial_ess:.2f}, which is less than the target ESS: {ess_t:.2f}.')
                print(f'Function values at the endpoints are f(a)={f_a:.2f} and f(b)={f_b:.2f}.')
                Warning = "The function f_ess does not have a root in the interval [a, b]."
                # raise ValueError("The function f_ess does not have a root in the interval [a, b].")
        else:
            while (b - a) > precision:
                beta = (a + b) / 2.0
                f_beta = f_ess(log_w_shift, beta, ess_t)
                if (f_a * f_beta < 0):
                    b = beta
                    f_b = f_beta
                else:
                    a = beta
                    f_a = f_beta
            beta = (a + b) / 2.0
        # --- Compute the final weights with the found beta
        scaled_weights = torch.exp(beta * log_w_shift)
        scaled_weights = scaled_weights / scaled_weights.sum()
        scaled_ess = effective_sample_size(scaled_weights)
    else:
        scaled_ess = initial_ess
        scaled_weights = weights
        beta = 1.0
        
    return scaled_weights, scaled_ess, beta, initial_ess