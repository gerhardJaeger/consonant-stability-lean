library(phytools)
library(plotrix)
library(dplyr)
library(tidyr)
library(phangorn)
library(ggplot2)
library(OUwie)
library(castor)


tree <- read.nexus("data/tree_pruned.nex")

# Generate Bounded Brownian Motion data

dir.create("data/simulated/", recursive = TRUE, showWarnings = FALSE)

for (s in c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)) {
  x <- fastBM(tree, sig2=s, a=4, bounds = c(1,10), internal=TRUE)
  
  x<-x[tree$tip.label]
  
  x<-as.data.frame(x)
  
  x <- dplyr::rename(x, series_markedness_fullness = x)
  
  #hist(x$series_markedness_fullness)
  
  write.csv(x, file=paste0("data/simulated/BBM_sig_",s,".csv"), row.names=TRUE)
}

# Generate Bounded OU Brownian Motion data

for (sd in c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)) {
  for (t_half in c(0.5, 1, 2, 5, 10, 15, 20, 30, 40, 50, 70, 100)) {
    dec <- log(2) / t_half
    OUwie_dat <- data.frame(species = tree$tip.label, regime = "1", x = 0)

    OUwie_dat$x <- simulate_ou_model(tree, stationary_mean = 4, stationary_std = sd, decay_rate = dec)$tip_states

    #hist(OUwie_dat$x)
    
    x <- OUwie_dat %>% dplyr::select(species, x) %>% dplyr::rename(series_markedness_fullness = x, glottocode = species)
    
    write.csv(x, file=paste0("data/simulated/OU_sd_",sd,"_thalf_",t_half,".csv"), row.names=FALSE)
  }
}


