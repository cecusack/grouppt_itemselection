#########################
### synthetic example ###
###    supp file 2    ###
#########################

#### libs ####
library(psych)
library(truncnorm)
library(clairecleans)
library(reshape2)
library(tidyverse)
library(mlVAR)
library(qgraph)
library(plyr)
library(ggbeeswarm)

#### synthetic dat ####
dat <- read.csv("supp2dat.csv")
items <-  seq(4,8,by =1) # let's say there are 5 nodes we're considering. Feeling tired, anxious, sad, concentration problems, binge eating
# all possible combinations of 3 items
all_combination <- combn(x = items, m = 3) # all_combination is an array containing the 10 unique 3-item symptom combinations
# note: because there are only 10 combinations, I will not randomly sample from here but instead estimate all 10.
# in the actual study there were over 18 million 8-item combinations, so we randomly selected 1,000 unique combinations

#### random-item network ####
boot <- dim(all_combination)[2] # set empty array to store all cent values
all_cents <- array(NA,dim = c(4,3,boot)) # create an empty 3D array

# set progress bar
pb <- txtProgressBar(min=0, max=boot, style=3,char="=")
init <- numeric(boot)
end <- numeric(boot)
pcc_plots <- list()
#### estimate random-item nets ####
for(j in 1:boot){
  
  print(j) # so we know which iteration we're on
  init[j] <- Sys.time()
  network_variance <- sd(melt(dat[,all_combination[,j]],id.vars=NULL)$value)
  
  # set empty matrix with 4 rows (item names, instrength, outstrength, and strength) and 8 columns (centrality estimates for 8 nodes)
  network_1 <- matrix(0, nrow=4, ncol=3,
                      dimnames = list(c("nodes", "InStrength", "OutStrength", "Strength"),
                                      colnames(dat)[all_combination[,j]]))
  
  network_1[1,] <- colnames(dat)[all_combination[,j]] # make first row items selected for permutation
  
  var_group <- mlVAR(dat, vars= colnames(dat[,all_combination[,j]]), idvar="ID", beepvar="beepvar", dayvar="dayvar", lags = 1,temporal = "correlated", estimator = "lmer", contemporaneous = "correlated", verbose = TRUE, scale = TRUE, scaleWithin = FALSE, AR = FALSE, MplusSave = TRUE, MplusName = "mlVAR", iterations = "(2000)", chains = nCores)
  
  pdcplot <- plot(var_group) # temporal networks
  tempcent <- if("InStrength" %in% names(centrality_auto(pdcplot)[[1]])) { # store instrength and outstrength estimates
    centrality_auto(pdcplot)[[1]][,c("InStrength", "OutStrength")]
  } else { # if instrength and outstrenth = 0, these columns won't be in centrality plot so
    data.frame(matrix(0, nrow=3,ncol=2)) # create 8x2 matrix of 0's
  }
  
  pccplot <- plot(var_group, "contemporaneous", edge.labels=TRUE) # contemporaneous network
  strength <- centrality_auto(pccplot)[[1]][3] # store strength centrality estimates in a vector
  
  # store centrality estimates into matrix 
  network_1[2,] <- unlist(tempcent[1]) # place instrength estimates in 2nd row
  network_1[3,] <- unlist(tempcent[2]) # place oustrength estimates in 3rd row
  network_1[4,] <- unlist(strength) # place strength centrality estimates in 4th row
  
  rm(pdcplot, pccplot) # remove network objects
  
  
  # store directed network results
  all_cents[,,j] <- network_1 # store all persons network_1 in the all_cents array
  
  end[j] <- Sys.time()
  setTxtProgressBar(pb,j)
  time <- round(lubridate::seconds_to_period(sum(end-init)),0)
  est <- boot*(mean(end[end!=0]-init[init!=0]))-time
  remaining <- round(lubridate::seconds_to_period(est),0)
  cat(paste(" // Execution time:", time,
            " // Estimated time remaining:", remaining), "")
}

close(pb)

#### tidy up ####
all_cents <- plyr::alply(all_cents, 3, .dims = TRUE) # split on third dimension. keep dimension names
# set column names
all_cents <- lapply(all_cents, function(x){ # this function makes row 1 (node items) column names across the list
  colnames(x) <- x[1,]
  x <- x[-1,] # then removes row 1
  return(x)
})

# set rownames
all_cents <- lapply(all_cents, function(x){ row.names(x)<-c("InStrength", "OutStrength", "Strength"); x})


flatten <- function(x){  
  islist <- sapply(x, class) %in% "list"
  r <- c(x[!islist], unlist(x[islist],recursive = F))
  if(!sum(islist))return(r)
  flatten(r)
}

nodenames_iter <- unlist(Map(colnames,flatten(all_cents)))
table(nodenames_iter) # 220 - 276

# class(all_cents[[1]]) # matrix array
all_cents <- lapply(all_cents, function(x) {if(any(class(x)=="matrix")) as.data.frame(x) else x}) # make df
# make numeric instead of character so that I can rank
all_cents<-lapply(all_cents, function(x) {
  x[] <- lapply(x, as.numeric)
  x
})

for(i in seq_along(all_cents)){
  all_cents[[i]] <- as.data.frame(t(all_cents[[i]])) # transpose
  all_cents[[i]]$rank <- NA # add rank col
  all_cents[[i]]$rank <- ifelse(!is.na(all_cents[[i]]$Strength), rank(-all_cents[[i]]$Strength, na.last = NA), NA) # order rank descending. it's fine with neg and pos
  all_cents[[i]]$names <- rownames(all_cents[[i]]) # add col for rownames
}

all_cents <- lapply(all_cents, function(x){ row.names(x)<-NULL; x}) # no longer need rownames
all_cents <- do.call(rbind, all_cents) # put in a df
all_cents <- all_cents %>% relocate(names) # names first

all_cents$iter <- rownames(all_cents) # add a column for combination number
all_cents$iter <- sub("^(\\d+).*", "\\1", all_cents$iter) 

#### a priori network ####
##### highest means #####
item_sel(dat[4:ncol(dat)], 3) # highest means are anxious, diff_conc, sad
