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
mod <- "Model/IRT//Diagonal/"
f.name <- paste(mod, "IRT_diagonal", sep="")

# --- Number of replicated experiment ---
n.exp <- 40

# --- Targeted Acceptance Probability ---
range.p <- seq(60, 95, by=5)

# --- Indices related to parameter of interest ---
ind.eta <- 2:101
ind.a <- 148:167
ind.b <- 125:144
ind.hyper <- c(123,145,146,147)
ind.lp <- 168
# c(1,102,124)  index for log(sigma_theta), log(sigma_a), log(sigma_b)
# 103:122 index for phi -- log(a)?!

# --- Data Processing: ESS ---
ESS.NUTS.eta <- ESS.vec(mod, n.exp, range.p, ind.eta, ind.lp)
ESS.NUTS.a <- ESS.vec(mod, n.exp, range.p, ind.a, ind.lp)
ESS.NUTS.b <- ESS.vec(mod, n.exp, range.p, ind.b, ind.lp)
ESS.NUTS.hyper <- ESS.vec(mod, n.exp, range.p, ind.hyper, ind.lp)

ESS.eHMC.eta <- ESS.vec(mod, n.exp, range.p, ind.eta, ind.lp, algo = "eHMC")
ESS.eHMC.a <- ESS.vec(mod, n.exp, range.p, ind.a, ind.lp, algo = "eHMC")
ESS.eHMC.b <- ESS.vec(mod, n.exp, range.p, ind.b, ind.lp, algo = "eHMC")
ESS.eHMC.hyper <- ESS.vec(mod, n.exp, range.p, ind.hyper, ind.lp, algo = "eHMC")

result.ESS <- rbind(data.frame(var = "eta: min(ESS) / gradient", transform.output(ESS.NUTS.eta)),
                    data.frame(var = "a: min(ESS) / gradient", transform.output(ESS.NUTS.a)),
                    data.frame(var = "b: min(ESS) / gradient", transform.output(ESS.NUTS.b)),
                    data.frame(var = "eta: min(ESS) / gradient", transform.output(ESS.eHMC.eta)),
                    data.frame(var = "a: min(ESS) / gradient", transform.output(ESS.eHMC.a)),
                    data.frame(var = "b: min(ESS) / gradient", transform.output(ESS.eHMC.b)))
result.ESS[, 1] <- as.factor(result.ESS[, 1])
result.ESS[, 2] <- as.factor(result.ESS[, 2])

summary.ESS <- rbind(data.frame(var = "eta: min(ESS) / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS.eta$output[,2,], 2, median)),
                     data.frame(var = "a: min(ESS) / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS.a$output[,2,], 2, median)),
                     data.frame(var = "b: min(ESS) / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS.b$output[,2,], 2, median)),
                     data.frame(var = "eta: min(ESS) / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC.eta$output[,2,], 2, median)),
                     data.frame(var = "a: min(ESS) / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC.a$output[,2,], 2, median)),
                     data.frame(var = "b: min(ESS) / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC.b$output[,2,], 2, median)))

ggplot(data = result.ESS, mapping = aes(x = p, y = min.ess)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2)) +
  geom_point(alpha = 1, size = 0.4) +
  facet_grid(var ~ sampler, 
             scales = "free_y", labeller = label_parsed) +
  geom_line(data = summary.ESS, mapping = aes(x = p, y = val), color = "red", size = 1) +
  xlab("Targeted acceptance probability:" ~p[0])


x <- apply(ESS.NUTS.eta$output[,2,], 1, max)
y <- apply(ESS.eHMC.eta$output[,2,], 1, max)
summary(x)
summary(y)
boxplot(x,y)
mean(x)
sd(x)
mean(y)
sd(y)
mean(y)/mean(x)

plot(apply(ESS.NUTS.b$output[, 2,], 2, mean))
lines(apply(ESS.eHMC.b$output[, 2,], 2, mean))


