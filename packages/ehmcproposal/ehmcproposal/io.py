import pandas as pd
import torch

from pathlib import Path

def save_checkpoint(
    beta_history: list,
    ess_history: list,
    scaled_ess_history: list,
    elapsed_time: float,
    path: str | Path,
    y = None,
    log_w_is = None,
):
    """
    Save a checkpoint to a file, including histories and optionally the sample and log weights.
    
    Args:
        beta_history (list): History of beta values.
        ess_history (list): History of effective sample sizes.
        scaled_ess_history (list): History of scaled effective sample sizes.
        elapsed_time (float): Total elapsed time for the process.
        path (str or Path): The file path to save the checkpoint.
        y (torch.Tensor, optional): The sample tensor to save. Default is None.
        log_w_is (torch.Tensor, optional): The log weights tensor to save. Default is None.
    """
    path = Path(path)
    
    checkpoint = {
        "beta_history": beta_history,
        "ess_history": ess_history,
        "scaled_ess_history": scaled_ess_history,
        "elapsed_time": elapsed_time
    }
    if y is not None:
        checkpoint['y'] = y.detach().cpu()
        
    if log_w_is is not None:
        checkpoint['log_w_is'] = log_w_is.detach().cpu()
        
    torch.save(checkpoint, path)
    




def load_checkpoint(path: str | Path, device: torch.device):
    """
    Load a checkpoint from a file and move tensors to the specified device.
    
    Args:
        path (str or Path): The file path to load the checkpoint from.
        device (torch.device): The device to move the loaded tensors to.
    
    Returns:
        dict: A dictionary containing the loaded checkpoint data, with tensors moved to the specified device.
    """
    path = Path(path)
    
    checkpoint = torch.load(path, map_location = device, weights_only = True)

    if "y" in checkpoint and checkpoint["y"] is not None:
        checkpoint["y"] = checkpoint["y"].to(device)

    if "log_w_is" in checkpoint and checkpoint["log_w_is"] is not None:
        checkpoint["log_w_is"] = checkpoint["log_w_is"].to(device)

    return checkpoint





def save_history_to_parquet(history: list, filename: str, column: str):
    """
    Save a history list to a Parquet file with a specified column name.
    
    Args:
        history (list): The history data to save.
        filename (str): The name of the output Parquet file.
        column (str): The name of the column in the DataFrame.
    """
    df = pd.DataFrame({column: history})
    df.to_parquet(filename, engine = "pyarrow", index = False)
    print(f"Saved {column} history to {filename}")
    




def save_warmup_sample(
    y: torch.Tensor, 
    weights: torch.Tensor, 
    filename: str,
    col_names: list[str] = None
):
    """
    Save the warmup sample and corresponding weights to a Parquet file.
    
    Args:
        y (torch.Tensor): The sample tensor of shape (n_samples, n_features).
        weights (torch.Tensor): The corresponding weights tensor of shape (n_samples,).
        filename (str): The name of the output Parquet file.
        col_names (list[str], optional): The names of the columns in the DataFrame. Default is None.
    """
    df = pd.DataFrame(
        y.detach().cpu().numpy(),
        columns = col_names if col_names is not None else [f"param_{i + 1}" for i in range(y.shape[1])],
    )
    df["w"] = weights.detach().cpu().numpy()
    df.to_parquet(filename, index = False)