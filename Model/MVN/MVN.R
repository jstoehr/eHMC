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
mod <- "Model/MVN/Diagonal/"
f.name <- paste(mod, "MVN_diagonal.pdf", sep="")

# --- Number of replicated experiment ---
n.exp <- 40

# --- Targeted Acceptance Probability ---
range.p <- seq(60, 95, by=5)

# --- Indices related to parameter of interest ---
ind <- 1:100
ind.lp <- 101

# --- Data Processing: ESS ---
ESS.NUTS <- ESS.vec(mod, n.exp, range.p, ind, ind.lp)
ESS.eHMC <- ESS.vec(mod, n.exp, range.p, ind, ind.lp, algo = "eHMC")

result.ESS <- rbind(transform.output(ESS.NUTS), transform.output(ESS.eHMC))
summary.ESS <- rbind(data.frame(criterion = "min(ESS) / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS$output[,2,], 2, median)),
                     data.frame(criterion = "min(ESS) / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC$output[,2,], 2, median)))


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
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2)) +
  geom_point(alpha = 1, size = 0.4) +
  facet_grid(factor(criterion) ~ factor(sampler), scales = "free_y", labeller = label_parsed) +
  geom_line(data = df.sum, mapping = aes(x = p, y = val), color = "red", size = 1) +
  xlab("Targeted acceptance probability:" ~p[0])

ggsave(f.name, device = "pdf", width = 14, height = 10, units = "cm", dpi = 600)


