#metric: minESS, meanESS, medianESS, maxESS       -- ess of mcmcse
#        minESS, meanESS, medianESS, maxESS       -- effectiveSize of coda
#        ESS of lp (mcmcse), ESS of lp (coda)
#        ESJD_x, ESJD_lp
#        Computation Cost, Tuned epsilon, Exptected L, empirical accept probability

# --- Libraries ---
library(ggplot2)

# --- Path ---
setwd("~/git/eHMC")
source("main/process_data.R")
mod <- "Model/SV/Diagonal/"
f.name <- paste(mod, "SV_diagonal", sep="")

# --- Number of replicated experiment ---
n.exp <- 40

# --- Targeted Acceptance Probability ---
range.p <- seq(60, 95, by=5)

# --- Indices related to parameter of interest ---
ind.x <- 1:1000
ind.phi <- 1004
ind.kappa <- 1005
ind.sigma <- 1006
ind.lp <- 1007

# --- Data Processing: ESS ---
ESS.NUTS.x <- ESS.vec(mod, n.exp, range.p, ind.x, ind.lp)
ESS.NUTS.phi <- ESS.vec(mod, n.exp, range.p, ind.phi, ind.lp)
ESS.NUTS.kappa <- ESS.vec(mod, n.exp, range.p, ind.kappa, ind.lp)
ESS.NUTS.sigma <- ESS.vec(mod, n.exp, range.p, ind.sigma, ind.lp)

ESS.eHMC.x <- ESS.vec(mod, n.exp, range.p, ind.x, ind.lp, algo = "eHMC")
ESS.eHMC.phi <- ESS.vec(mod, n.exp, range.p, ind.phi, ind.lp, algo = "eHMC")
ESS.eHMC.kappa <- ESS.vec(mod, n.exp, range.p, ind.kappa, ind.lp, algo = "eHMC")
ESS.eHMC.sigma <- ESS.vec(mod, n.exp, range.p, ind.sigma, ind.lp, algo = "eHMC")

result.ESS <- rbind(data.frame(var = "x: min(ESS) / gradient", transform.output(ESS.NUTS.x)),
                    data.frame(var = "phi: ESS / gradient", transform.output(ESS.NUTS.phi)),
                    data.frame(var = "kappa: ESS / gradient", transform.output(ESS.NUTS.kappa)),
                    data.frame(var = "sigma: ESS / gradient", transform.output(ESS.NUTS.sigma)),
                    data.frame(var = "x: min(ESS) / gradient", transform.output(ESS.eHMC.x)),
                    data.frame(var = "phi: ESS / gradient", transform.output(ESS.eHMC.phi)),
                    data.frame(var = "kappa: ESS / gradient", transform.output(ESS.eHMC.kappa)),
                    data.frame(var = "sigma: ESS / gradient", transform.output(ESS.eHMC.sigma)))
result.ESS[, 1] <- as.factor(result.ESS[, 1])
result.ESS[, 2] <- as.factor(result.ESS[, 2])

summary.ESS <- rbind(data.frame(var = "x: min(ESS) / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS.x$output[,2,], 2, median)),
                     data.frame(var = "phi: ESS / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS.phi$output[,2,], 2, median)),
                     data.frame(var = "kappa: ESS / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS.kappa$output[,2,], 2, median)),
                     data.frame(var = "sigma: ESS / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS.sigma$output[,2,], 2, median)),
                     data.frame(var = "x: min(ESS) / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC.x$output[,2,], 2, median)),
                     data.frame(var = "phi: ESS / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC.phi$output[,2,], 2, median)),
                     data.frame(var = "kappa: ESS / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC.kappa$output[,2,], 2, median)),
                     data.frame(var = "sigma: ESS / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC.sigma$output[,2,], 2, median)))


ggplot(data = result.ESS, mapping = aes(x = p, y = min.ess)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2)) +
  geom_point(alpha = 1, size = 0.4) +
  facet_grid(var ~ sampler, 
             scales = "free_y", labeller = label_parsed) +
  geom_line(data = summary.ESS, mapping = aes(x = p, y = val), color = "red", size = 1) +
  xlab("Targeted acceptance probability:" ~p[0])

ggsave(paste(f.name, "ESS.pdf", sep = "_"), device = "pdf", width = 14, height = 20, units = "cm", dpi = 600)


# --- Data Processing: ESS ---
ESJD.NUTS <- ESJD.vec(mod, n.exp, range.p, ind, ind.lp)
ESJD.eHMC <- ESJD.vec(mod, n.exp, range.p, ind, ind.lp, algo = "eHMC")

result.ESJD <- rbind(transform.output(ESJD.NUTS), transform.output(ESJD.eHMC))
summary.ESJD <- rbind(data.frame(criterion = "ESJD / gradient", sampler = "NUTS", p = range.p/100,
                                 val = apply(ESJD.NUTS$output[,2,], 2, median)),
                      data.frame(criterion = "ESJD / gradient", sampler = "eHMC", p = range.p/100,
                                 val = apply(ESJD.eHMC$output[,2,], 2, median)))

# --- Global data.frame ---
df.ans <- rbind(data.frame(criterion = "min(ESS) / gradient", 
                           result.ESS[, 1:2], val = result.ESS[, 3]),
                data.frame(criterion = "ESJD / gradient", 
                           result.ESJD[, 1:2], val = result.ESJD[, 3]))

df.sum <- rbind(summary.ESS, summary.ESJD)


ggplot(data = df.ans, mapping = aes(x = p, y = val)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank()) +
  geom_point(alpha = 1, size = 0.4) +
  facet_grid(rows = vars(factor(criterion)), cols = vars(factor(sampler)), scales = "free_y") +
  geom_line(data = df.sum, mapping = aes(x = p, y = val), color = "red", size = 1) +
  xlab("Targeted acceptance probability:" ~p[0])

ggsave(f.name, device = "pdf", width = 14, height = 10, units = "cm", dpi = 600)


