#metric: minESS, meanESS, medianESS, maxESS       -- ess of mcmcse
#        minESS, meanESS, medianESS, maxESS       -- effectiveSize of coda
#        ESS of lp (mcmcse), ESS of lp (coda)
#        ESJD_x, ESJD_lp
#        Computation Cost, Tuned epsilon, Exptected L, empirical accept probability
library(mcmcse)
setwd("~/git/eHMC/prHMC/Model/MVN/SimulationResult")

NUTS_ESS <- numeric(40)
eHMC_ESS <- numeric(40)
prHMC_ESS <- numeric(40)

temp <- get(load(file="Result.RData"))

for(i in 1:40){
  val <- temp[[i]]
  NUTS_ESS[i] <- min(val$NUTS_Summary$statistics$ess_by_mcmcse)/val$NUTS_Summary$Comp
  eHMC_ESS[i] <- min(val$eHMC_Summary$statistics$ess_by_mcmcse)/val$eHMC_Summary$Comp
  #prHMC_ESS[i] <- min(val$eHMCq_Summary_Summary$statistics$ess_by_mcmcse)/val$NUTS_Summary$Comp
}


################################ Data Processing #################################
NUTS_KS <- matrix(0,nrow=40,ncol=100)
eHMC_KS <- matrix(0,nrow=40,ncol=100)
eHMCq_KS <- matrix(0,nrow=40,ncol=100)
eHMCu_KS <- matrix(0,nrow=40,ncol=100)
MCHMC_KS <- matrix(0,nrow=40,ncol=100)

for(i in 1:40)
{
  rt <- temp[[i]]
  NUTS_KS[i,] <- rt$NUTS_KS
  eHMC_KS[i,] <- rt$eHMC_KS
  eHMCq_KS[i,] <- rt$eHMCq_KS
  eHMCu_KS[i,] <- rt$eHMCu_KS
  MCHMC_KS[i,] <- rt$CVHMC_KS
}

boxplot(apply(NUTS_KS, 1, max), apply(eHMC_KS, 1, max), apply(MCHMC_KS, 1, max))

ESSTransformed <- NULL
ESSTransformed <- rbind(ESSTransformed, cbind(0, apply(NUTS_KS, 1, max)))
ESSTransformed <- rbind(ESSTransformed, cbind(1, apply(eHMC_KS, 1, max)))
ESSTransformed <- rbind(ESSTransformed, cbind(2, apply(eHMCq_KS, 1, max)))
ESSTransformed <- rbind(ESSTransformed, cbind(3, apply(eHMCu_KS, 1, max)))
ESSTransformed <- rbind(ESSTransformed, cbind(4, apply(MCHMC_KS, 1, max)))

ESSTransformed <- data.frame(ESSTransformed)
colnames(ESSTransformed) <-  c("Sampler", "y")
IndexSampler <- which(ESSTransformed[,"Sampler"] == 0)
ESSTransformed[IndexSampler, "Sampler"] <- "NUTS"
IndexSampler <- which(ESSTransformed[,"Sampler"] == 1)
ESSTransformed[IndexSampler, "Sampler"] <- "eHMC"
IndexSampler <- which(ESSTransformed[,"Sampler"] == 2)
ESSTransformed[IndexSampler, "Sampler"] <- "eHMCq"
IndexSampler <- which(ESSTransformed[,"Sampler"] == 3)
ESSTransformed[IndexSampler, "Sampler"] <- "eHMCu"
IndexSampler <- which(ESSTransformed[,"Sampler"] == 4)
ESSTransformed[IndexSampler, "Sampler"] <- "prHMC"

neworder <- c("NUTS","eHMC","eHMCq","eHMCu","prHMC")
library(plyr)  ## or dplyr (transform -> mutate)
library(ggplot2)
ESSTransformed <- arrange(transform(ESSTransformed,
                                    Sampler=factor(Sampler,levels=neworder)),Sampler)

ggplot(data=ESSTransformed, mapping = aes(x=Sampler, y=y)) +
  geom_boxplot(mapping=aes(color=Sampler, fill=Sampler)) +
  ylab("Kolmogorov-Smirnov distance")+
  ggtitle("MVN: Max KS of" ~theta) +
  theme(plot.title = element_text(hjust = 0.5))

#ggplot(data=ESSTransformed,mapping=aes(x=x,y=y)) + 
#  geom_point(alpha=1,size=0.4) +
#  geom_line(data=ESSSummary, mapping=aes(x=x, y=y), color="red", size=1) +
#  facet_grid(cols = vars(factor(Sampler))) + xlab("Desired accept probability:" ~delta) + 
#  ylab("ESS per gradient") + 
#  ggtitle("Logit: Min ESS of" ~theta) +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=3)) + 
#  ylim(c(0,0.15))


#neworder <- c("NUTS","eHMC","eHMCq","eHMCu","MCHMC")
Ind1 <- which(ESSTransformed[,"Sampler"] == "NUTS")
Ind2 <- which(ESSTransformed[,"Sampler"] == "eHMC")
Ind3 <- which(ESSTransformed[,"Sampler"] == "prHMC")
neworder <- c("NUTS","eHMC","prHMC")
ESSTransformed <- ESSTransformed[c(Ind1, Ind2, Ind3),]
ESSTransformed <- arrange(transform(ESSTransformed,
                                    Sampler=factor(Sampler,levels=neworder)),Sampler)
ggplot(data = ESSTransformed, mapping = aes(x = Sampler, y = y)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none") +
  geom_boxplot(mapping = aes(color = Sampler, fill = Sampler))
#ggtitle("Logit: Max KS of" ~theta) +
# theme(plot.title = element_text(hjust = 0.5))

ggsave("MVN_diagonal_KS.pdf", device = "pdf", width = 14, height = 7, units = "cm", dpi = 600)



