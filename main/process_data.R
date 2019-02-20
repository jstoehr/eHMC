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
ESS <- function(L, index, index.lp, algo = "NUTS_Summary", method = "ess_by_mcmcse") {
  # --- Extract result corresponding to the algo ---
  temp.L <- L[[algo]]
  # --- Computationnal cost ---
  cost <- temp.L$Comp
  # --- ESS method ---
  temp <- temp.L[["statistics"]][[method]][index]
  lp <- temp.L[["statistics"]][[method]][index.lp]
  
  return(data.frame(
    min.ess = min(temp)/cost,
    mean.ess = mean(temp)/cost,
    median.ess = median(temp)/cost,
    max.ess = max(temp)/cost,
    ess.lp = lp/cost
  ))
}

ESS.vec <- function (mod, n.exp, range.p, ind, ind.lp, 
                     algo = "NUTS_Summary", method = "ess_by_mcmcse"){
  ans <- array(0, c(n.exp, 6, length(range.p)), dimnames = list(
    NULL, c("p", "min.ess", "mean.ess", "median.ess", "max.ess", "ess.lp"), NULL
  ))
  ans[ , 1, ] <- rep(range.p/100, each = n.exp)
  
  if(algo == "NUTS_Summary"){
    algo.name <- "NUTS"
  } else {
    algo.name <- "eHMC"
  }
  
  for(k in 1:length(range.p)) {
    temp <- get(load(file = paste(mod, toString(range.p[k]), "/Result.RData", sep = "")))
    for(i in 1:n.exp) {
      ans[i, -1, k] <- as.numeric(ESS(temp[[i]], ind, ind.lp, algo, method))
    }
  }
  return(list(sampler = algo.name, ess.method = method, output = ans))
}


# --- ESJD functions ---
ESJD <- function(L, algo = "NUTS_Summary", param = "ESJD_para") {
  # --- Extract result corresponding to the algo ---
  temp.L <- L[[algo]]
  # --- Computationnal cost ---
  cost <- temp.L$Comp
  
  return(data.frame(
    esjd = temp.L[["statistics"]][[param]]/cost,
    esjd.lp = temp.L$statistics$ESJD_lp/cost
  ))
}

ESJD.vec <- function (mod, n.exp, range.p, algo = "NUTS_Summary", param = "ESJD_para"){
  ans <- array(0, c(n.exp, 3, length(range.p)), dimnames = list(
    NULL, c("p", "esjd", "esjd.lp"), NULL
  ))
  ans[ , 1, ] <- rep(range.p/100, each = n.exp)
  
  if(algo == "NUTS_Summary"){
    algo.name <- "NUTS"
  } else {
    algo.name <- "eHMC"
  }
  
  for(k in 1:length(range.p)) {
    temp <- get(load(file = paste(mod, toString(range.p[k]), "/Result.RData", sep = "")))
    for(i in 1:length(temp)) {
      ans[i, -1, k] <- as.numeric(ESJD(temp[[i]], algo, param))
    }
  }
  return(list(sampler = algo.name, output = ans))
}


# --- Summaries ---
resume <- function(L, index, index.lp, algo = "NUTS_Summary", L.name = "L_trace") {
  # --- Extract result corresponding to the algo ---
  temp.L <- L[[algo]]
  if(algo == "NUTS_Summary"){
    dist.L <- temp.L[[L.name]]
  } else {
    dist.L <- L[[L.name]][["L_emp"]]
  }
  
  return(data.frame(
    emp.p = temp.L$accept,
    eps = temp.L$epsilon,
    cost = temp.L$Comp,
    t(as.matrix(summary(dist.L)))
  ))
}

resume.vec <- function (mod, n.exp, range.p, ind, ind.lp, algo = "NUTS_Summary"){
  ans <- array(0, c(n.exp, 10, length(range.p)), dimnames = list(
    NULL, c("p", "emp.p", "eps", "cost", "min(L)", "1st.Q(L)", "median(L)", "mean(L)", "3rd.Q(L)", "max(L)"), NULL
  ))
  ans[ , 1, ] <- rep(range.p/100, each = n.exp)
  
  if(algo == "NUTS_Summary"){
    algo.name <- "NUTS"
    L.name <- "L_trace"
  } else {
    algo.name <- "eHMC"
    L.name <- "LearnL_Summary"
  }
  
  for(k in 1:length(range.p)) {
    temp <- get(load(file = paste(mod, toString(range.p[k]), "/Result.RData", sep = "")))
    for(i in 1:length(temp)) {
      ans[i, -1, k] <- as.numeric(resume(temp[[i]], ind, ind.lp, algo, L.name))
    }
  }
  return(list(sampler = algo.name, output = ans))
}
