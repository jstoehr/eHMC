# -------------------- 
# --- BANANA MODEL
# --------------------
import torch
from pathlib import Path

def get_data(
    data_dir: str | Path = None, 
    device = torch.device("cpu"), 
    dtype = torch.float32
):
    return {
        "sigma1": torch.tensor(10.0, device = device, dtype = dtype),
        "lambda": torch.tensor(0.03, device = device, dtype = dtype),
        "lag": torch.tensor(100.0, device = device, dtype = dtype),
        "sigma2": torch.tensor(1.0, device = device, dtype = dtype),
    }



def get_pars_name(data = None) -> list[str]:
    return ["theta1", "theta2"]



def rprop_init(
    n: int, 
    data = None, 
    seed: int | None = None, 
    device = torch.device("cpu"), 
    dtype = torch.float32
) -> tuple[torch.Tensor, torch.Tensor]:
    """Sampling from the initial proposal distribution."""
    if seed is not None:
        torch.manual_seed(seed)

    if data is None:
        data = get_data(device = device, dtype = dtype)

    sigma1 = data["sigma1"]
    lam = data["lambda"]
    lag = data["lag"]

    # --- Gaussian distribution parameters
    mu = torch.stack([
        torch.tensor(0.0, device = sigma1.device, dtype = sigma1.dtype),
        -lam * lag,
    ])

    cov = torch.diag(torch.full(
        (2,),
        sigma1 * sigma1,
        device = sigma1.device,
        dtype = sigma1.dtype,
    ))

    # --- Sample and compute log probabilities from the multivariate normal distribution
    dist = torch.distributions.MultivariateNormal(
        loc = mu,
        covariance_matrix = cov,
    )

    sample = dist.sample((n,))
    lp__ = dist.log_prob(sample)

    return sample, lp__