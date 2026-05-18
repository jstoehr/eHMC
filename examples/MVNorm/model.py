# -------------------- 
# --- MVNORM MODEL
# --------------------
import math
import torch
from pathlib import Path

def get_data(
    data_dir: str | Path = None, 
    device = torch.device("cpu"), 
    dtype = torch.float32
    ):
    d = 100
    rho = 0.99

    A = torch.zeros((d, d), device=device, dtype=dtype)

    for i in range(d):
        for j in range(i, d):
            A[i, j] = rho ** (j - i)

    A = A + A.T - torch.eye(d, device=device, dtype=dtype)
    A_inv = torch.linalg.inv(A)

    return {
        "d": d,
        "A": A,
        "A_inv": A_inv,
    }



def get_pars_name(data = None) -> list[str]:
    if data is None:
        data = get_data()

    return [f"x[{i}]" for i in range(1, data["d"] + 1)]



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
        data = get_data(device=device, dtype=dtype)

    d = data["d"]
    A = data["A"].to(device=device, dtype=dtype)

    # --- Mixture distribution parameters
    mu = torch.zeros(d, device=A.device, dtype=A.dtype)

    eigvals, eigvecs = torch.linalg.eigh(A)

    dt_1 = torch.where(
        eigvals > 1.0,
        torch.log(eigvals),
        eigvals,
    )

    dt_2 = torch.where(
        eigvals > 1.0,
        0.6 * eigvals,
        1.1 * eigvals,
    )

    sigma_1 = eigvecs @ torch.diag(dt_1) @ eigvecs.T
    sigma_2 = eigvecs @ torch.diag(dt_2) @ eigvecs.T

    alpha = 0.1
    
    # --- Sample from the mixture distribution
    n_1 = int(math.floor(alpha * n))
    n_2 = n - n_1

    dist_1 = torch.distributions.MultivariateNormal(mu, covariance_matrix=sigma_1)
    dist_2 = torch.distributions.MultivariateNormal(mu, covariance_matrix=sigma_2)

    sample_1 = dist_1.sample((n_1,))
    sample_2 = dist_2.sample((n_2,))

    sample = torch.cat([sample_1, sample_2], dim=0)

    # --- Compute log probabilities for the mixture distribution
    log_p1 = dist_1.log_prob(sample)
    log_p2 = dist_2.log_prob(sample)

    lp__ = torch.logsumexp(
        torch.stack([
            math.log(alpha) + log_p1,
            math.log(1.0 - alpha) + log_p2,
        ]),
        dim=0,
    )

    return sample, lp__