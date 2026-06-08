from typing import Callable, Any
from dataclasses import dataclass

import torch


@dataclass
class FlowConfig:
    """
    Configuration parameters for the flow model used in learning the PDF proxy.
    
    Attributes:
        test_fraction (float): The fraction of data to be used as the test set.
        architecture (str): The architecture of the flow model, either "MAF" or "RNVP".
        transforms (int): The number of flow transforms to use in the model.
        hidden_features (Tuple[int, ...]): The number of hidden features in each layer of the flow model.
        randperm (bool): Whether to use random permutations in the MAF architecture.
        max_epochs (int): The maximum number of training epochs for the flow model.
        patience (int): The number of epochs with no improvement after which training will be stopped.
        learning_rate (float): The learning rate for the optimizer.
        batch_size (int): The batch size for training the flow model.
        activation (Callable): The activation function to use in the flow model.
        device (torch.device): The device to use for training the flow model (e.g., "cpu", "cuda", "mps").
    """
    test_fraction: float = 0.2
    architecture: str = "MAF"
    transforms: int = 5
    hidden_features: tuple[int, ...] = (64, 64)
    randperm: bool = True
    max_epochs: int = 512
    patience: int = 30
    learning_rate: float = 1e-3
    batch_size: int = 256
    activation: Callable = torch.nn.ReLU # torch.nn.Tanh
    device: torch.device = torch.device("cpu")
   


@dataclass
class ProposalConfig:
    """
    Configuration parameters for the PMC proposal adaptation process.
    
    Attributes:
        n_warmup (int): The number of warmup samples to generate before the main sampling phase.
        n_train (int): The number of training samples to use for learning the flow model.
        ess_train (float): The target effective sample size (ESS) for the training phase, expressed as a fraction of n_train.
        ess_target (float | None): The target effective sample size (ESS) for the PMC process. If None, it will be set to ess_train * n_train.
        seed (int): The random seed for reproducibility.
        max_tries (int): The maximum number of tries for the PMC process to achieve the target ESS.
        resampling (int): The resampling strategy to use (0, for no resampling, 1 for multinomial resampling).
    """
    n_warmup: int = 1000
    n_train: int = 1000
    ess_train: float = 0.8
    ess_target: float | None = None
    seed: int = 1
    max_tries: int = 10
    resampling: int = 1
    data: dict[str, Any] | None = None
    device: torch.device = torch.device("cpu")
    dtype: torch.dtype = torch.float32