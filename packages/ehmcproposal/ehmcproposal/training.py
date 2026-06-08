import torch
import torch.utils.data as data
from tqdm.auto import tqdm
import zuko
from zuko.distributions import MultivariateNormal
from zuko.lazy import UnconditionalDistribution

from .config import FlowConfig


def train_proposal(
    y: torch.Tensor,
    cov: torch.Tensor | None = None,
    config: FlowConfig | None = None,
    flow = None,
    show_progress: bool = True,
):
    """
    Train a Masked Autoregressive Flow (MAF) model to estimate the density of y 
    
    Args:
        y (torch.Tensor): The input data tensor of shape (n_samples, n_features).
        cov (Optional[torch.Tensor]): Optional covariance matrix for the base distribution.
        config (FlowConfig): Configuration parameters for the flow model.
        architecture (str): The architecture of the flow model, either "MAF" or "RNVP".
        flow (Optional[zuko.flows.MAF or zuko.flows.RNVP]): An existing flow model to be trained. If None, a new model is created.
        
    Returns:
        flow (zuko.flows.MAF or zuko.flows.RNVP): The trained flow model.
    """
    config = FlowConfig() if config is None else config

    y = y.detach()
    device = y.device
    n_samples = y.size(0)

    # --------------------------------------------------
    # --- Step 1: splitting data
    # --------------------------------------------------

    # --- Train set
    n_train = int((1.0 - config.test_fraction) * n_samples)

    if n_train == 0:
        raise ValueError("The training set is empty. Decrease config.test_fraction.")

    if n_train == n_samples:
        raise ValueError("The test set is empty. Increase config.test_fraction.")

    permuted_indices = torch.randperm(
        n_samples,
        device = device,
    )

    y_train = y[permuted_indices[:n_train], :]

    # --- Test set
    y_test = y[permuted_indices[n_train:], :]
    features = y.size(1)

    # --- Loading data
    trainloader = data.DataLoader(
        y_train,
        batch_size = config.batch_size,
        shuffle = True,
    )

    # --------------------------------------------------
    # --- Step 2: setup flow and optimizer
    # --------------------------------------------------

    # --- Setting the flow model
    if flow is None:
        if config.architecture == "MAF":
            flow = zuko.flows.MAF(
                features = features,
                transforms = config.transforms,
                hidden_features = config.hidden_features,
                randperm = config.randperm,
                activation = config.activation,
            )

        elif config.architecture == "RealNVP":
            flow = zuko.flows.RealNVP(
                features = features,
                transforms = config.transforms,
                hidden_features = config.hidden_features,
                activation = config.activation,
            )

        else:
            raise ValueError(
                f"Unknown flow architecture: {config.architecture}"
            )

    flow = flow.to(device)

    if cov is not None:
        # Setting the base distribution with a given covariance
        flow.base = UnconditionalDistribution(
            MultivariateNormal,
            torch.zeros(
                features,
                device = device,
                dtype = y.dtype,
            ),
            cov.to(device),
            buffer = True,
        )
    
    # --- Setting optimizer
    optimizer = torch.optim.Adam(
        flow.parameters(),
        lr = config.learning_rate,
    )

    # --------------------------------------------------
    # --- Step 3: training loop with early stopping
    # --------------------------------------------------

    best_loss = float("inf")
    best_state_dict = None
    patience_counter = 0

    iterator = range(config.max_epochs)

    if show_progress:
        iterator = tqdm(
            iterator,
            desc = f"Train {config.architecture}",
            unit = "epoch",
            leave = False,
        )

    flow.train()

    for epoch in iterator:
        # losses = []

        flow.train()

        for x in trainloader:
            loss = -flow().log_prob(x).mean()

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # losses.append(loss.detach())

        # losses = torch.stack(losses)
        # train_loss_mean = losses.mean().item()
        # train_loss_std = losses.std().item()

        flow.eval()

        with torch.no_grad():
            test_loss = -flow().log_prob(y_test).mean().item()

        if test_loss < best_loss:
            best_loss = test_loss
            best_state_dict = {
                key: value.detach().clone()
                for key, value in flow.state_dict().items()
            }
            patience_counter = 0

        else:
            patience_counter += 1

        if show_progress:
            iterator.set_postfix(
                validation_loss = f"{test_loss:.4f}",
                best_loss = f"{best_loss:.4f}",
                patience = patience_counter,
            )

        if patience_counter >= config.patience:
            break

    if best_state_dict is not None:
        flow.load_state_dict(best_state_dict)

    flow.eval()

    return flow