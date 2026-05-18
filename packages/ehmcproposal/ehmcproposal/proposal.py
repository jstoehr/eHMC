import torch
import time
from pathlib import Path
from .config import FlowConfig, ProposalConfig
from .reweighting import normalize_log_weights, effective_sample_size,  reweight_sample
from .linalg import is_positive_definite
from .training import train_proposal
from .r_interface import apply_r_log_pdf
from .io import save_checkpoint, load_checkpoint

def sample_initial_proposal(
    rprop_init, 
    target_log_pdf, 
    n: int, 
    pconfig: ProposalConfig = ProposalConfig()
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, float]:
    """
    Sample from the initial proposal and compute IS log-weights.
    
    Args:
        rprop_init: A function that samples from the initial proposal distribution and computes log probabilities.
        target_log_pdf: An R function that computes the log PDF of the target distribution.
        n (int): The number of samples to draw from the initial proposal.
        pconfig (ProposalConfig): Configuration parameters for the PMC process, used for setting the random seed.
    
    Returns:
        A tuple containing:
            - y (torch.Tensor): The samples drawn from the initial proposal distribution.
            - log_w_is (torch.Tensor): The log importance weights for the samples.
            - w_is (torch.Tensor): The normalized importance weights for the samples.
            - ess (float): The effective sample size (ESS) computed from the normalized importance weights.
    """

    y, lp__ = rprop_init(
        n = n, 
        data = pconfig.data, 
        seed = pconfig.seed, 
        device = pconfig.device, 
        dtype = pconfig.dtype
    )

    target_lp__ = apply_r_log_pdf(target_log_pdf, y)

    log_w_is = target_lp__ - lp__
    w_is = normalize_log_weights(log_w_is)
    ess = effective_sample_size(w_is)

    return y, log_w_is, w_is, ess





def resample_from_initial_proposal(
    rprop_init,
    target_log_pdf,
    pconfig: ProposalConfig = ProposalConfig(),
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Generate final eHMC warmup sample from the initial proposal.
    
    Args:
        rprop_init: A function that samples from the initial proposal distribution and computes log probabilities.
        target_log_pdf: An R function that computes the log PDF of the target distribution.
        pconfig (ProposalConfig): Configuration parameters for the PMC process.
    
    Returns:
        A tuple containing:
            - y (torch.Tensor): The final sample tensor after resampling from the initial proposal.
            - w (torch.Tensor): The normalized importance weights for the final sample.
    """

    if pconfig.resampling == 0:
        n = pconfig.n_warmup
    elif pconfig.resampling == 1:
        n = pconfig.n_train
    else:
        raise ValueError("pconfig.resampling must be 0 or 1.")

    y, log_w, w, _ = sample_initial_proposal(
        rprop_init = rprop_init,
        target_log_pdf = target_log_pdf,
        n = n,
        pconfig = pconfig
    )

    if pconfig.resampling == 0:
        return y, w
    
    indices = torch.multinomial(w, pconfig.n_warmup, replacement = True)

    w_resampled = torch.full(
        (pconfig.n_warmup,),
        1.0 / pconfig.n_warmup,
        dtype = w.dtype,
        device = w.device,
    )

    return y[indices, :], w_resampled





def pmc_proposal(
    rprop_init,
    target_log_pdf,
    fconfig: FlowConfig = FlowConfig(),
    pconfig: ProposalConfig = ProposalConfig(),
    checkpoint_path: str | Path = None,
    flow_path: str | Path = None,
    resume: bool = False
) -> tuple[torch.Tensor, torch.Tensor, list, list, list, float]:
    """
    Perform Population Monte Carlo (PMC) to adapt a proposal distribution for HMC calibration.
    
    Args:
        rprop_init: A function that samples from the initial proposal distribution and computes log probabilities.
        target_log_pdf: An R function that computes the log PDF of the target distribution.
        fconfig (FlowConfig): Configuration parameters for the flow model.
        pconfig (ProposalConfig): Configuration parameters for the PMC process.
        checkpoint_path (str or Path, optional): The file path to save checkpoints during the PMC process. Default is None (no checkpoints saved).
        flow_path (str or Path, optional): The file path to save the final trained flow model. Default is None (model not saved).
    
    Returns:
    tuple:
        - y (torch.Tensor): The final sample tensor after PMC adaptation.
        - w (torch.Tensor): The normalized importance weights for the final sample.
        - ess_history (list): History of effective sample sizes during the PMC process.
        - scaled_ess_history (list): History of scaled effective sample sizes during the PMC process.
        - beta_history (list): History of beta values during the PMC process.
        - elapsed_time (float): Total elapsed time for the PMC process.
    """
    
    # --- Set the seed for reproducibility (if provided)
    if pconfig.seed is not None:
        print(f"Setting the Python seed to {pconfig.seed} for reproducibility.")
        torch.manual_seed(pconfig.seed)
        
    # --- Initialize containeers
    checkpoint_path = Path(checkpoint_path) if checkpoint_path is not None else None
    flow_path = Path(flow_path) if flow_path is not None else None
    n_tries = 0
    beta = 0.0
    elapsed_time = 0.0
    flow = None
    lag = 0
    has_resumed = False
    
    
    # --------------------------------------------------
    # --- Step 0: check for checkpoint and resume if available
    # --------------------------------------------------
        
    if resume and checkpoint_path is not None and checkpoint_path.exists():
        print(f"Resuming from checkpoint: {checkpoint_path}")
        
        checkpoint = load_checkpoint(checkpoint_path, device=pconfig.device)
        
        y = checkpoint["y"]
        log_w_is = checkpoint["log_w_is"]
        
        ess_past = checkpoint["ess_history"]
        scaled_ess_past = checkpoint["scaled_ess_history"]
        beta_past = checkpoint["beta_history"]
        
        ess_history = ess_past + [None] * pconfig.max_tries
        scaled_ess_history = scaled_ess_past + [None] * pconfig.max_tries
        beta_history = beta_past + [None] * pconfig.max_tries
        
        elapsed_time = checkpoint["elapsed_time"]
        
        lag = len(ess_past) - 1
        
        if flow_path is not None and flow_path.exists():
            print(f"Loading flow model from: {flow_path}")
            flow = torch.load(flow_path, map_location = pconfig.device, weights_only = False)
            
        has_resumed = True
    
    
    if not has_resumed:
        # --------------------------------------------------
        # --- Step 1: evaluate initial proposal
        # --------------------------------------------------

        y, log_w_is, w_is, ess_init = sample_initial_proposal(
            rprop_init = rprop_init,
            target_log_pdf = target_log_pdf,
            n = pconfig.n_train,
            pconfig = pconfig
        )
    
        ess_history = [None] * (pconfig.max_tries + 1)
        scaled_ess_history = [None] * pconfig.max_tries
        beta_history = [None] * (pconfig.max_tries + 1)
        
        ess_history[0] = ess_init
        beta_history[0] = beta
    
        # --------------------------------------------------
        # --- Step 2: use initial proposal if good enough
        # --------------------------------------------------
    
        if ess_init >= pconfig.ess_target:
            print(
                f"Initial proposal accepted: ESS = {ess_init:.2f}, "
                f"target = {pconfig.ess_target:.2f}."
            )

            y_warmup, w_warmup = resample_from_initial_proposal(
                rprop_init = rprop_init,
                target_log_pdf = target_log_pdf,
                pconfig = pconfig,
            )

            return (
                y_warmup,
                w_warmup,
                ess_history[:1],
                scaled_ess_history[:0],
                beta_history[:1],
                elapsed_time
            )
        
    # --------------------------------------------------
    # --- Step 3: train flow if initial proposal is not enough
    # --------------------------------------------------
    
    # --- Reweight the sample to achieve the target ESS
    scaled_weights, scaled_ess, beta, current_ess = reweight_sample(
        log_w_is, 
        pconfig.ess_train * pconfig.n_train, 
        beta_init = beta
    )
    
    print(f"Initial ESS: {current_ess:.2f}")
    
    while current_ess < pconfig.ess_target and n_tries < pconfig.max_tries:
        start = time.perf_counter()
        n_tries += 1
        # --- Compute the emprical covariance
        cov_emp = torch.cov(y.T, aweights=scaled_weights)
        if not is_positive_definite(cov_emp):
            print("Covariance matrix is not positive definite. Using default value.")
            cov_emp = None
        # --- Resample the data
        indices = torch.multinomial(scaled_weights, pconfig.n_train, replacement=True)
        y_resampled = y[indices,:]
        # --- Train the MAF model with the resampled data
        flow = train_proposal(
            y_resampled,
            cov = cov_emp,
            config = fconfig,
            flow = flow if n_tries % 20 != 0 else None
        )
        with torch.no_grad():
            # --- Sample from the MAF model
            y = flow().sample((pconfig.n_train,))  
            # --- Compute the importance weights for the MAF sample
            proposal_lp__ = flow().log_prob(y)
        
        target_lp__ = apply_r_log_pdf(target_log_pdf, y)
        log_w_is = target_lp__ - proposal_lp__
        # --- Reweight the sample
        scaled_weights, scaled_ess, beta, current_ess = reweight_sample(
            log_w_is, 
            pconfig.ess_train * pconfig.n_train, 
            beta_init = 0.0
        )
        end = time.perf_counter()
        elapsed_time += end - start
        # --- Update the history
        ess_history[n_tries + lag] = current_ess
        scaled_ess_history[n_tries - 1 + lag] = scaled_ess
        beta_history[n_tries + lag] = beta
        print(
            f"Try {n_tries}: ESS = {current_ess:.2f}, "
            f"scaled ESS = {scaled_ess:.2f}, "
            f"Target ESS = {pconfig.ess_target:.2f}, beta = {beta:.4f}, "
            f"elapsed = {elapsed_time:.2f}s"
        )
        if checkpoint_path is not None and n_tries % 3 == 0:
            save_checkpoint(
                beta_history = beta_history[: n_tries + 1 + lag],
                ess_history = ess_history[: n_tries + 1 + lag],
                scaled_ess_history = scaled_ess_history[: n_tries + lag],
                elapsed_time = elapsed_time,
                path = checkpoint_path,
                y = y,
                log_w_is = log_w_is,
            )
    
    # --------------------------------------------------
    # --- Step 4: save flow/checkpoint
    # --------------------------------------------------
    if flow_path is not None:
        torch.save(flow, flow_path)

    if checkpoint_path is not None:
        save_checkpoint(
            beta_history = beta_history[: n_tries + 1 + lag],
            ess_history = ess_history[: n_tries + 1 + lag],
            scaled_ess_history = scaled_ess_history[: n_tries + lag],
            elapsed_time = elapsed_time,
            path = checkpoint_path,
            y = y,
            log_w_is = log_w_is
        )
        
    if current_ess < pconfig.ess_target:
        raise RuntimeError(
            f"PMC proposal failed: final ESS = {current_ess:.2f}, "
            f"target ESS = {pconfig.ess_target:.2f}."
        )

    # --------------------------------------------------
    # --- Step 5: final warmup sample from trained flow
    # --------------------------------------------------
    
    print("Sampling from the flow for eHMC calibration.")
    if pconfig.resampling == 0:
        n_warmup = pconfig.n_warmup
    elif pconfig.resampling == 1:
        n_warmup = pconfig.n_train
    else:
        raise ValueError("pconfig.resampling must be 0 or 1.")
            
    with torch.no_grad():
        y = flow().sample((n_warmup,))
        proposal_lp__ = flow().log_prob(y)
                        
    target_lp__ = apply_r_log_pdf(target_log_pdf, y)                
    log_w = target_lp__ - proposal_lp__
    w = normalize_log_weights(log_w)

    if pconfig.resampling == 0:
        return(
            y, 
            w, 
            ess_history[: n_tries + 1 + lag], 
            scaled_ess_history[: n_tries + lag], 
            beta_history[: n_tries + 1 + lag], 
            elapsed_time
        )

    indices = torch.multinomial(w, pconfig.n_warmup, replacement=True)
    w_resampled = torch.full(
        (pconfig.n_warmup,),
        1.0 / pconfig.n_warmup,
        dtype = w.dtype,
        device = w.device,
    )
    
    return (
        y[indices, :], 
        w_resampled, 
        ess_history[: n_tries + 1 + lag], 
        scaled_ess_history[: n_tries + lag], 
        beta_history[: n_tries + 1 + lag], 
        elapsed_time
    )