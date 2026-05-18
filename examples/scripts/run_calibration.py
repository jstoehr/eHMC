import os
import argparse
import ast
import shutil
import torch
import rpy2.robjects as robjects
import importlib.util
from pathlib import Path

from ehmcproposal.config import FlowConfig, ProposalConfig
from ehmcproposal.proposal import pmc_proposal
from ehmcproposal.io import save_warmup_sample, save_history_to_parquet

def load_python_model(model_file: Path):
    """
    Dynamically load a Python module from a given file path.
        
    Args:
        model_file (Path): The file path to the Python module to be loaded.
    """
    spec = importlib.util.spec_from_file_location("model", model_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    return module



def load_r_target_log_pdf(model_name: str, model_dir: Path):
    """
    Load the target log PDF function from R for a given model.
    
    Args:
        model_name (str): The name of the model for which to load the target log PDF function.
        model_dir (Path): The directory where the model's R scripts and data are located.
    """
    model_dir = Path(model_dir)

    robjects.r("library(ehmcexamples)")
    robjects.r["source"](str(model_dir / "model.R"))

    robjects.globalenv["model_name"] = model_name
    robjects.globalenv["model_dir"] = str(model_dir)

    robjects.r("""
        .ehmc_data <- get_data(data_dir = model_dir)
        .ehmc_target_log_pdf <- ehmcexamples::get_target_log_pdf(
            model = model_name,
            data = .ehmc_data
        )
    """)

    return robjects.globalenv[".ehmc_target_log_pdf"]



def checkpoint_from_lower_target(
    proposal_dir: Path,
    checkpoint_path: Path,
    flow_path: Path,
    seed: int,
    n_train: int,
    ess_train: float,
    ess_target: float
):
    """
    Initialize the checkpoint for a given ESS target by copying from a checkpoint with a lower ESS target if it exists.
    
    Args:
        proposal_dir (Path): The directory where proposal checkpoints are stored.
        checkpoint_path (Path): The file path for the checkpoint to be initialized.
        flow_path (Path): The file path for the flow model to be initialized.
        seed (int): The random seed used in training.
        n_train (int): The number of training samples used in training.
        ess_train (float): The ESS achieved during training.
        ess_target (float): The target ESS for which to initialize the checkpoint.
    """
    if checkpoint_path.exists():
        return

    pattern = (
        f"seed_{seed}_train_{n_train}_{int(100 * ess_train)}"
        f"_ess_*_checkpoint.pt"
    )

    candidates = []

    for path in proposal_dir.glob(pattern):
        name = path.name

        try:
            old_ess = int(name.split("_ess_")[1].split("_checkpoint.pt")[0])
        except Exception:
            continue

        if old_ess < ess_target:
            candidates.append((old_ess, path))

    if not candidates:
        return

    # Use the largest lower ESS target.
    old_ess, source_path = max(candidates, key=lambda x: x[0])

    shutil.copy2(source_path, checkpoint_path)
    
    source_flow = proposal_dir / source_path.name.replace(
        "_checkpoint.pt",
        "_flow.pt",
    )

    if source_flow.exists() and flow_path is not None:
        shutil.copy2(source_flow, flow_path)

    print(
        f"Initialized checkpoint for ESS target {int(ess_target)} "
        f"from lower target {old_ess}: {source_path}"
    )



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--repo_root", type = str, default = ".", help = "Root directory of the repository (default: current directory)")
    parser.add_argument("--model", type = str, required = True, help = "Model name (e.g., 'MVNorm', 'BLP', etc.)")
    parser.add_argument("--seed", type = int, default = 123, help = "Random seed for reproducibility (default: 123)")
    parser.add_argument("--n_train", type = int, default = 500, help = "Number of training samples (default: 500)")
    parser.add_argument("--ess_train", type = float, default = 0.8, help = "Target ESS for training expressed as a fraction of n_train (default means 0.8 * n_train)")
    parser.add_argument("--ess_target", type = float, default = None, help = "Target ESS (default: None)")
    parser.add_argument("--n_warmup", type = int, default = 1000, help = "Number of warmup samples in eHMC (default: 1000)")
    parser.add_argument("--resampling", type = int, default = 1, help = "Resampling step from proposal (default: 1)")
    parser.add_argument("--max_tries", type = int, default = 20, help = "Maximum number of tries to improve the proposal (default: 20)")
    parser.add_argument("--resume", type = int, default = 1, help = "Whether to resume from the last checkpoint if it exists.")
    parser.add_argument("--device", type = str, default = "cpu", help = "Device to use for computation (default: 'cpu')")
    parser.add_argument("--test_fraction", type = float, default = 0.2, help = "The proportion of the sample to use for testing.")
    parser.add_argument("--transforms", type = int, default = 5, help = "The number of transformation blocks in the MAF.")
    parser.add_argument("--hidden_features", type = str, default = "(32, 32)", help = "A tuple defining the hidden layer sizes in the MAF (e.g., '(64, 64)').")
    parser.add_argument("--randperm", type = int, default = 1, help = "Whether to use random permutation.")
    parser.add_argument("--max_epochs", type = int, default = 2048, help = "The maximum number of epochs to train for.")
    parser.add_argument("--patience", type = int, default = 20, help = "Lag for the early stopping rule.")
    parser.add_argument("--learning_rate", type = float, default = 1e-3, help = "Learning rate of Adam optimizer.")
    parser.add_argument("--batch_size", type = int, default = 256, help = "The batch size to use during training.")

    args = parser.parse_args()
    
    # --- Set the default ESS target
    if args.ess_target is None:
        args.ess_target = 0.8 * args.n_train
    # --- Convert hidden_features to a tuple
    hidden_features = ast.literal_eval(args.hidden_features)
    
    # --------------------------------------------------
    # --- Set up directories and file paths
    # --------------------------------------------------
    
    repo_root = Path(args.repo_root).resolve()
    
    model_dir = repo_root / "examples" / args.model
    
    proposal_dir = model_dir / "proposal"
    proposal_dir.mkdir(parents = True, exist_ok = True)
    
    if args.resampling ==  0:
        out_warmup_dir = model_dir / "warmup_wo_resampling"
    elif args.resampling ==  1:
        out_warmup_dir = model_dir / "warmup_with_resampling"
    else:
        raise ValueError("resampling must be 0 or 1.")
    
    out_warmup_dir.mkdir(parents = True, exist_ok = True)
    
    basename = f"seed_{args.seed}_train_{args.n_train}_{int(100 * args.ess_train)}_ess_{int(args.ess_target)}"
    
    checkpoint_path = proposal_dir / f"{basename}_checkpoint.pt"
    flow_path = proposal_dir / f"{basename}_flow.pt"
    ess_path = proposal_dir / f"{basename}_ess_history.parquet"
    scaled_ess_path = proposal_dir / f"{basename}_scaled_ess_history.parquet"
    beta_path = proposal_dir / f"{basename}_beta_history.parquet"
    
    warmup_path = out_warmup_dir / f"{basename}_calib_ehmc_w_{args.n_warmup}_r_{args.resampling}.parquet"
    
    if os.path.exists(warmup_path):
        return  # Exit if the calibration file already exists
    
    # --------------------------------------------------
    # --- Load the model and data
    # --------------------------------------------------
    
    model_py = load_python_model(model_dir / "model.py")
    target_log_pdf = load_r_target_log_pdf(args.model, model_dir)
    
    data = model_py.get_data(
        data_dir = model_dir,
        device = torch.device(args.device),
        dtype = torch.float32,
    )
    
    # --------------------------------------------------
    # --- Configurations
    # --------------------------------------------------
    
    device = torch.device(args.device)
    fconfig = FlowConfig(
        test_fraction = args.test_fraction,
        transforms = args.transforms,
        hidden_features = hidden_features,
        randperm = bool(args.randperm),
        max_epochs = args.max_epochs,
        patience = args.patience,
        learning_rate = args.learning_rate,
        batch_size = args.batch_size,
        device = device
    )
    
    pconfig = ProposalConfig(
        n_warmup = args.n_warmup,
        n_train = args.n_train,
        ess_train = args.ess_train,
        ess_target = args.ess_target,
        seed = args.seed,
        max_tries = args.max_tries,
        resampling = args.resampling,
        data = data,
        device = device,
        dtype = torch.float32
    )
    
    # --------------------------------------------------
    # --- Run the PMC proposal adaptation
    # --------------------------------------------------
    
    checkpoint_from_lower_target(
        proposal_dir = proposal_dir,
        checkpoint_path = checkpoint_path,
        flow_path = flow_path,
        seed = args.seed,
        n_train = args.n_train,
        ess_train = args.ess_train,
        ess_target = args.ess_target
    )
    
    y, w, ess_history, scaled_ess_history, beta_history, elapsed_time = pmc_proposal(
        rprop_init = model_py.rprop_init,
        target_log_pdf = target_log_pdf,
        fconfig = fconfig,
        pconfig = pconfig,
        checkpoint_path = checkpoint_path,
        flow_path = flow_path,
        resume = bool(args.resume),
    )
    
    print(f"Elapsed time: {elapsed_time:.2f}s")
    
    # --------------------------------------------------
    # --- Save outputs to parquet files for usage in R
    # --------------------------------------------------
    
    save_history_to_parquet(ess_history, ess_path, "ess")
    save_history_to_parquet(scaled_ess_history, scaled_ess_path, "scaled_ess")
    save_history_to_parquet(beta_history, beta_path, "beta")
    
    save_warmup_sample(y, w, warmup_path, col_names = model_py.get_pars_name(data))
    print(f"Warmup sample saved to: {warmup_path}")
    

if __name__ ==  "__main__":

    main()