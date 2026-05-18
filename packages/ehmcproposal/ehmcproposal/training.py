from typing import Optional

import torch
import torch.utils.data as data
import zuko
from zuko.distributions import MultivariateNormal
from zuko.lazy import UnconditionalDistribution

from .config import FlowConfig


def train_proposal(
    y: torch.Tensor,
    cov: Optional[torch.Tensor] = None,
    config: FlowConfig = FlowConfig(),
    flow = None
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
    
    y = y.detach()
    n_samples = y.size(0)
    
    # --------------------------------------------------
    # --- Step 1: splitting data
    # --------------------------------------------------

    # --- Train set
    n_train = int((1 - config.test_fraction) * n_samples)
    permuted_indices = torch.randperm(n_samples)
    y_train = y[permuted_indices[:n_train],:]
   
    # --- Test set
    y_test = y[permuted_indices[n_train:],:]
    features = y.size(1)
        
    # --- Loading data
    trainloader = data.DataLoader(y_train, batch_size = config.batch_size, shuffle = True)

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
                activation = config.activation
            )
            if y_train.is_cuda:
                flow = flow.cuda()
            elif y_train.is_mps:
                flow = flow.mps()
        else:
            raise ValueError(f"Unknown flow architecture: {config.architecture}")
    
    if cov is not None:
        # Setting the base distribution with a given covariance
        flow.base = UnconditionalDistribution(
            MultivariateNormal,
            torch.zeros(features, device = y.device),
            cov,
            buffer = True
        )
    
    # --- Setting optimizer
    optimizer = torch.optim.Adam(flow.parameters(), lr = config.learning_rate)

    # --------------------------------------------------
    # --- Step 3: training loop with early stopping
    # --------------------------------------------------
    
    best_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(config.max_epochs):
        # losses = []
        for x in trainloader:
            loss = -flow().log_prob(x).mean()
            loss.backward()
            
            optimizer.step()
            optimizer.zero_grad()
            
            # losses.append(loss.detach())
        # losses = torch.stack(losses)
        # train_loss_mean = losses.mean().item()
        # train_loss_std = losses.std().item()
        # print(f'Epoch {epoch}: Training loss {train_loss_mean} ± {train_loss_std}')
            
        with torch.no_grad():
            test_loss = -flow().log_prob(y_test).mean()
            # print(f'Epoch {epoch + 1}: Validation loss {test_loss} (patience: {patience_counter + 1})')

        if test_loss < best_loss:
            best_loss = test_loss
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= config.patience:
            break

    return flow