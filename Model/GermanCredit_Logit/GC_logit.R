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
mod <- "Model/GermanCredit_Logit/Diagonal/"
f.name <- paste(mod, "GC_logit_diagonal.pdf", sep="")

# --- Number of replicated experiment ---
n.exp <- 40

# --- Targeted Acceptance Probability ---
range.p <- seq(60, 95, by=5)

# --- Indices related to parameter of interest ---
ind <- 1:25
ind.lp <- 26

# --- Data Processing: ESS ---
ESS.NUTS <- ESS.vec(mod, n.exp, range.p, ind, ind.lp)
ESS.eHMC <- ESS.vec(mod, n.exp, range.p, ind, ind.lp, algo = "eHMC_Summary")

result.ESS <- rbind(transform.output(ESS.NUTS), transform.output(ESS.eHMC))
summary.ESS <- rbind(data.frame(criterion = "min(ESS) / gradient", sampler = "NUTS", p = range.p/100,
                                val = apply(ESS.NUTS$output[,2,], 2, median)),
                     data.frame(criterion = "min(ESS) / gradient", sampler = "eHMC", p = range.p/100,
                                val = apply(ESS.eHMC$output[,2,], 2, median)))


# --- Data Processing: ESS ---
ESJD.NUTS <- ESJD.vec(mod, n.exp, range.p)
ESJD.eHMC <- ESJD.vec(mod, n.exp, range.p, algo = "eHMC_Summary")

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

x <- apply(ESS.NUTS$output[,2,], 1, max)
y <- apply(ESS.eHMC$output[,2,], 1, max)
mean(x)
sd(x)
mean(y)
sd(y)
mean(y)/mean(x)

x <- resume.vec(mod, n.exp, range.p, ind, ind.lp, algo = "NUTS_Summary")
y <- resume.vec(mod, n.exp, range.p, ind, ind.lp, algo = "eHMC_Summary")
x$output[,,5]
y$output[,,5]

# --- Old version ---
# p1 <- ggplot(data = result.ESS[, 1:3], mapping = aes(x = p, y = min.ess)) + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
#   geom_point(alpha = 1, size = 0.4) +
#   geom_line(data = summary.ESS, mapping = aes(x = p, y = val), color = "red", size = 1) +
#   facet_grid(cols = vars(factor(sampler))) + 
#   #xlab("Targeted acceptance probability:" ~p[0]) + 
#   ylab("min(ESS)/gradient") + 
#   #ggtitle("Minimum ESS of" ~theta) +
#   #theme(plot.title = element_text(hjust = 0.5)) +
#   #scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=3)) + 
#   ylim(c(0,0.18))
# 
# p2 <- ggplot(data = result.ESJD[, 1:3], mapping = aes(x = p, y = esjd)) + 
#   geom_point(alpha = 1, size = 0.4) +
#   geom_line(data = summary.ESJD, mapping = aes(x = p, y = val), color = "red", size = 1) +
#   facet_grid(cols = vars(factor(sampler))) + 
#   xlab("Targeted acceptance probability:" ~p[0]) + 
#   ylab("ESJD/gradient") + 
#   #ggtitle("ESJD of" ~theta) +
#   #theme(plot.title = element_text(hjust = 0.5)) +
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=3)) + 
#   ylim(c(0,0.18))
# 
# grid.arrange(p1, p2, ncol = 1, heights = c(2, 2))