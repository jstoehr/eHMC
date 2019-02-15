#metric: minESS, meanESS, medianESS, maxESS       -- ess of mcmcse
#        minESS, meanESS, medianESS, maxESS       -- effectiveSize of coda
#        ESS of lp (mcmcse), ESS of lp (coda)
#        ESJD_x, ESJD_lp
#        Computation Cost, Tuned epsilon, Exptected L, empirical accept probability

# --- Useful functions ---
transform.output <- function(L) {
  ans <- L$output[, ,1]
  for(i in 2:dim(L$output)[3]){
    ans <- rbind(ans, L$output[, , i])
  }
  return(data.frame(sampler = L$sampler, ans))
}

# --- ESS functions ---
ESS <- function(L, index, index.lp, algo = "NUTS", method = "mcmcse") {
  # --- Extract result corresponding to the algo ---
  if(algo == "NUTS"){
    temp.L <- L$NUTS_Summary
  } else if(algo == "eHMC"){
    temp.L <- L$eHMC_Summary
  } else if(algo == "eHMCq"){
    temp.L <- L$eHMCq_Summary
  } else {
    temp.L <- L$eHMCu_Summary
  }
  # --- Computationnal cost ---
  cost <- temp.L$Comp
  # --- ESS method ---
  if(method == "mcmcse") {
    temp <- temp.L$statistics$ess_by_mcmcse[index]
    lp <- temp.L$statistics$ess_by_mcmcse[index.lp]/cost
  } else {
    temp <- temp.L$statistics$ess_by_coda[index]
    lp <- temp.L$statistics$ess_by_coda[index.lp]/cost
  }
  
  return(data.frame(
    min.ess = min(temp)/cost,
    mean.ess = mean(temp)/cost,
    median.ess = median(temp)/cost,
    max.ess = max(temp)/cost,
    ess.lp = lp
  ))
}

ESS.vec <- function (mod, n.exp, range.p, ind, ind.lp, algo = "NUTS", method = "mcmcse"){
  ans <- array(0, c(n.exp, 6, length(range.p)), dimnames = list(
    NULL, c("p", "min.ess", "mean.ess", "median.ess", "max.ess", "ess.lp"), NULL
  ))
  ans[ , 1, ] <- rep(range.p/100, each = n.exp)
  for(k in 1:length(range.p)) {
    temp <- get(load(file = paste(mod, toString(range.p[k]), "/Result.RData", sep = "")))
    for(i in 1:length(temp)) {
      val <- temp[[i]]
      ans[i, -1, k] <- as.numeric(ESS(val, ind, ind.lp, algo, method))
    }
  }
  return(list(sampler = algo, ess.method = method, output = ans))
}


# --- ESJD functions ---
ESJD <- function(L, index, index.lp, algo = "NUTS") {
  # --- Extract result corresponding to the algo ---
  if(algo == "NUTS"){
    temp.L <- L$NUTS_Summary
  } else if(algo == "eHMC"){
    temp.L <- L$eHMC_Summary
  } else if(algo == "eHMCq"){
    temp.L <- L$eHMCq_Summary
  } else {
    temp.L <- L$eHMCu_Summary
  }
  # --- Computationnal cost ---
  cost <- temp.L$Comp
  
  return(data.frame(
    esjd = temp.L$statistics$ESJD_para/cost,
    esjd.lp = temp.L$statistics$ESJD_lp/cost
  ))
}

ESJD.vec <- function (mod, n.exp, range.p, ind, ind.lp, algo = "NUTS"){
  ans <- array(0, c(n.exp, 3, length(range.p)), dimnames = list(
    NULL, c("p", "esjd", "esjd.lp"), NULL
  ))
  ans[ , 1, ] <- rep(range.p/100, each = n.exp)
  for(k in 1:length(range.p)) {
    temp <- get(load(file = paste(mod, toString(range.p[k]), "/Result.RData", sep = "")))
    for(i in 1:length(temp)) {
      val <- temp[[i]]
      ans[i, -1, k] <- as.numeric(ESJD(val, ind, ind.lp, algo))
    }
  }
  return(list(sampler = algo, output = ans))
}


# --- Summaries ---
resume <- function(L, index, index.lp, algo = "NUTS") {
  # --- Extract result corresponding to the algo ---
  if(algo == "NUTS"){
    temp.L <- L$NUTS_Summary
    dist.L <- temp.L$L_trace
  } else if(algo == "eHMC"){
    temp.L <- L$eHMC_Summary
    dist.L <- temp.L$LearnL_Summary$L_emp
  } else if(algo == "eHMCq"){
    temp.L <- L$eHMCq_Summary
    dist.L <- temp.L$LearnL_Summary$L_emp
  } else {
    temp.L <- L$eHMCu_Summary
    dist.L <- temp.L$LearnL_Summary$L_emp
  }
  
  return(data.frame(
    emp.p = temp.L$accept,
    eps = temp.L$epsilon,
    cost = temp.L$Comp,
    t(as.matrix(summary(dist.L)))
  ))
}

resume.vec <- function (mod, n.exp, range.p, ind, ind.lp, algo = "NUTS"){
  ans <- array(0, c(n.exp, 10, length(range.p)), dimnames = list(
    NULL, c("p", "emp.p", "eps", "cost", "min(L)", "1st.Q(L)", "median(L)", "mean(L)", "3rd.Q(L)", "max(L)"), NULL
  ))
  ans[ , 1, ] <- rep(range.p/100, each = n.exp)
  for(k in 1:length(range.p)) {
    temp <- get(load(file = paste(mod, toString(range.p[k]), "/Result.RData", sep = "")))
    for(i in 1:length(temp)) {
      val <- temp[[i]]
      ans[i, -1, k] <- as.numeric(resume(val, ind, ind.lp, algo))
    }
  }
  return(list(sampler = algo, output = ans))
}
