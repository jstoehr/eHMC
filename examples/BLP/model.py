# ----------------- 
# --- BLP MODEL
# -----------------
import torch
from pathlib import Path

def _read_numeric_table(path: str | Path):
    rows = []

    with open(path, "r", encoding = "utf-8") as file:
        for line in file:
            stripped = line.strip()

            if not stripped:
                continue

            rows.append([float(x) for x in stripped.split()])

    return rows



def get_data(
    data_dir: str | Path = None, 
    device = torch.device("cpu"), 
    dtype = torch.float32
    ):
    if data_dir is None:
        data_dir = Path(__file__).resolve().parent
    else:
        data_dir = Path(data_dir)

    data_file = data_dir / "german_data_numeric.csv"

    if not data_file.exists():
        raise FileNotFoundError(f"Missing data file: {data_file}")

    df = torch.tensor(
        _read_numeric_table(data_file),
        device = device,
        dtype = dtype,
    )

    n_obs = df.shape[0]
    k = df.shape[1]

    X = torch.cat(
        [
            torch.ones((n_obs, 1), device = device, dtype = dtype),
            df[:, : (k - 1)],
        ],
        dim = 1,
    )

    y = 2.0 - df[:, k - 1]

    return {
        "N": n_obs,
        "K": k,
        "X": X,
        "y": y,
        "sig": torch.tensor(1.0, device = device, dtype = dtype),
    }



def get_pars_name(data = None) -> list[str]:
    if data is None:
        data = get_data()

    return [f"theta[{i}]" for i in range(1, data["K"] + 1)]



def rprop_init_deprecated(
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

    k = data["K"]
    sig = data["sig"]

    # --- Gaussian distribution parameters
    mu = torch.zeros(k, device = sig.device, dtype = sig.dtype)
    cov = (sig * sig) * torch.eye(k, device = sig.device, dtype = sig.dtype)
    
    # --- Sample and compute log probabilities from the multivariate normal distribution
    dist = torch.distributions.MultivariateNormal(
        loc = mu,
        covariance_matrix = cov,
    )

    sample = dist.sample((n,))
    lp__ = dist.log_prob(sample)

    return sample, lp__



def rprop_init(
    n: int,
    data=None,
    seed: int | None = None,
    device=torch.device("cpu"),
    dtype=torch.float32,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Sampling from the initial proposal distribution."""
    if seed is not None:
        torch.manual_seed(seed)

    if data is None:
        data = get_data(device=device, dtype=dtype)

    X = data["X"]
    y = data["y"]
    sig = data["sig"]
    k = data["K"]

    # --- Prior precision
    prior_prec = torch.eye(k, device=X.device, dtype=X.dtype) / (sig * sig)

    # --- Approximate mode by ridge least squares
    XtX = X.T @ X
    Xty = X.T @ y
    precision = XtX + prior_prec

    cov_prior = torch.linalg.inv(precision)
    mu_prior = cov_prior @ Xty

    # --- Inflate covariance to avoid being too narrow
    inflation = 4.0
    cov = inflation * cov_prior

    # --- Sample and compute log probabilities from the multivariate normal distribution
    dist = torch.distributions.MultivariateNormal(
        loc=mu_prior,
        covariance_matrix=cov,
    )

    sample = dist.sample((n,))
    lp__ = dist.log_prob(sample)

    return sample, lp__