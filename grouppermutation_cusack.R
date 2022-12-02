#-------------------------------------------------------
# Central symptoms stability analysis at the group-level 
# 8-node networks (permuted) for PT pilot study

# Claire E. Cusack - 06/10/2022
# code adapted from Sean Kelley (sensitivity_general_update.R)
#-------------------------------------------------------

# rm(list=ls()) # clear env
#### libs and dat ####
library(reshape2)
library(tidyverse)
library(mlVAR)
library(qgraph)
library(plyr)
library(ggbeeswarm)
library(here)
here()

#### network combinations ####
# column numbers for target items
# items <-  seq(9,42,by =1) # items are in cols 9 thru 42
# all possible combinations of 8 items
# all_combination <- combn(x = items, m = 8)
# dim(all_combination) # whew 8 by 18,156,204

# write.csv(all_combination,file = 'node8_all33.csv',row.names = F)

# all_10000 <- all_combination[,sample(ncol(all_combination),10000)]
# dim(all_10000) # 8 x 10000

# write.csv(all_10000,file = 'node8_10000.csv',row.names = F)

# select a subset of 1,000 8-node networks 
# all_1000 <- all_combination[,sample(ncol(all_combination),1000)]
# dim(all_1000) # 8 x 1000

# write.csv(all_1000,file = 'node8_1000.csv',row.names = F)

# setwd("/Users/clairecusack/Dropbox/PT group network paper (PT open series)/grouppt_dev/IdioImpute")
# dir.create("/Users/clairecusack/Dropbox/PT group network paper (PT open series)/grouppt_dev/GroupCent_results")

# esm data
dat <- read.csv("pt_imputed_hasinitialqdat.csv")
dim(dat) # 6060 x 42
length(unique(dat$ID)) # 71
names(dat)

# read combination randomly sampled 1000
all_combination <- read.csv("node8_1000.csv")

# read baseline
baseline <- read.csv("initquestion_idsclean.csv")
# baseline <- read.csv("baselineforsk.csv") # for sk no identifiers. just study id and ede-q
baseline <- baseline %>% dplyr::rename(ID=unique_id)
# load("allcent_boot10.RData") # the test set 10 combinations if you didn't want to run lines 41-93

#### group level mlVAR ####
boot <- dim(all_combination)[2] # set boot to equal number of permutations
# set empty array to store all cent values
all_cents <- array(NA,dim = c(4,8,boot)) # create an empty 3D array
# set progress bar
pb <- txtProgressBar(min=0, max=boot, style=3,char="=")
init <- numeric(boot)
end <- numeric(boot)

for(j in 1:boot){
  
  print(j) # so we know which iteration we're on
  init[j] <- Sys.time()
  network_variance <- sd(melt(dat[,all_combination[,j]],id.vars=NULL)$value)
  # dat[,all_combination[,j]] <- as.data.frame(lapply(dat[,all_combination[,j]], sample))
  # set empty matrix with 4 rows (item names, instrength, outstrength, and strength) and 8 columns (centrality estimates for 8 nodes)
  network_1 <- matrix(0, nrow=4, ncol=8,
                      dimnames = list(c("nodes", "InStrength", "OutStrength", "Strength"),
                                      colnames(dat)[all_combination[,j]]))
  
  network_1[1,] <- colnames(dat)[all_combination[,j]] # make first row items selected for permutation
  
  # var_group <- mlVAR(dat[,c(which(colnames(dat)=="ID"), all_combination[,j])], vars = colnames(dat[,all_combination[,j]]), idvar = "ID", lags = 1)
  
  var_group <- mlVAR(dat[,c(which(colnames(dat)=="ID"), which(colnames(dat)=="dayvar"):which(colnames(dat)=="beepvar"), all_combination[,j])],  vars= colnames(dat[,all_combination[,j]]), idvar="ID", beepvar="beepvar", dayvar="dayvar", lags = 1,temporal = "correlated", estimator = "lmer", contemporaneous = "correlated", verbose = TRUE, scale = TRUE, scaleWithin = FALSE, AR = FALSE, MplusSave = TRUE, MplusName = "mlVAR", iterations = "(2000)", chains = nCores)
  
  pdcplot <- plot(var_group) # temporal networks
  tempcent <- if("InStrength" %in% names(centrality_auto(pdcplot)[[1]])) { # store instrength and outstrength estimates
    centrality_auto(pdcplot)[[1]][,c("InStrength", "OutStrength")]
  } else { # if instrength and outstrenth = 0, these columns won't be in centrality plot so
    data.frame(matrix(0, nrow=8,ncol=2)) # create 8x2 matrix of 0's
  }
  
  pccplot <- plot(var_group, "contemporaneous") # contemporaneous network
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

#### save ####
resultsdir <- here("Group_cent")
dir.create(resultsdir)
name_cent= paste0(resultsdir,'/all_cents')
filetype= '.rds'
filename1= paste(name_cent, Sys.Date(), filetype, sep='')
saveRDS(all_cents, filename1)

#--------------------------------------------------------------------------------
# Distribution of rank
#--------------------------------------------------------------------------------
# which nodes were included and how many times 
flatten <- function(x){  
  islist <- sapply(x, class) %in% "list"
  r <- c(x[!islist], unlist(x[islist],recursive = F))
  if(!sum(islist))return(r)
  flatten(r)
}

nodenames_iter <- unlist(Map(colnames,flatten(all_cents)))
table(nodenames_iter)

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
all_cents$iter <- sub("^(\\d+).*", "\\1", all_cents$iter) # keep only the combination number

name_centdf= paste0(resultsdir,'/all_centsDF')
filename2= paste(name_centdf, Sys.Date(), filetype, sep='')
saveRDS(all_cents, filename2)

all_cents$names<-as.factor(all_cents$names) # char to factor
all_cents %>% group_by(names) %>% filter(n()<2) # which items were only included once

#### plot distribution ####
# all_cents <- readRDS("~/Dropbox/PT group network paper (PT open series)/grouppt_dev/Group_cent/all_centsDF2022-06-19.rds")
View(all_cents)
apriori_cent <- read_csv("ImputedGroupNet/contemporaneous centrality imp.csv")
names(apriori_cent) <- c("names", "betweenness", "closeness", "strength", "ei")
apriori_cent$rank <- rank(-apriori_cent$strength, na.last = NA)

cent_dist <- 
  all_cents %>% 
  group_by(names) %>% filter(n()>1) %>% ungroup() %>% # remove items only included once
  ggplot(aes(x=reorder(names,rank, mean),y=rank, color=names,fill=names)) +
  geom_violin(alpha=0.7) +
  # geom_quasirandom(groupOnX = TRUE, width = .2, aes(color = names)) + # may turn this off when we have 1000 or adjust alpha / flip them so violin is set high and points lower. tbd
  
  geom_point(data = apriori_cent, aes(x = names, y = rank, shape=8), size=1.3, color="black")+
  scale_shape_identity()+
  theme_classic()+
  labs(x="", y="Strength Centrality Rank Order") +
  ggtitle("Distribution of Strength Centrality Rank Order Across 8-Node Networks") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position ="none",
        panel.grid.major = element_line(color="lightgrey"),
        plot.title = element_text(hjust = 0.5)) # ,
#         panel.grid.major.y = element_blank())

cent_dist # see it
name_violin= paste0(resultsdir,'/centdist_violin_nojitter.png')
ggsave(filename = name_violin, device = "png", plot = cent_dist, width = 12, height = 5, units = "in")

#### plot strength centrality not just ranks ####
#### summarize ####
freq <- as.data.frame(table(all_cents$names))
names(freq) <- c("names", "total_n")
# View(all_cents[which(all_cents$names=="shame"),])
df <- all_cents %>% group_by(names, rank) %>% tally()
df <- left_join(df, freq) # add total n
rm(freq)
df$perc <- round((df$n/df$total_n)*100,2)
# df <- left_join(df, all_cents[-c(2:3)])
# save
name_sum=paste0(resultsdir,'/summaryranks')
filetype= '.csv'
filename= paste(name_sum, Sys.Date(), filetype, sep='')
write.csv(df,file=filename, row.names = FALSE)


cis <- all_cents %>%
  group_by(names) %>%
  summarise_at(vars(Strength), list(m = mean, sd=sd)) 
cis <- left_join(cis, freq)
cis <- cis %>% mutate(upper = m + 1.96*(sd/sqrt(total_n)), lower = m - 1.96*(sd/sqrt(total_n)))
all_cents <- left_join(all_cents, cis)

cent_dist2 <- all_cents %>% 
  group_by(names) %>% filter(n()>1) %>% ungroup() %>% # remove items only included once
  ggplot(aes(x=reorder(names,Strength, mean, decreasing = TRUE),y=Strength, color=names, fill=names)) +
  geom_violin(alpha=0.7) +
  # geom_point(stat="summary", fun.y="mean", color="black", size=.5)+
  # geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), color="black", width=0) +
  # stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
  #              colour = "black")+
  geom_point(data = apriori_cent, aes(x = names, y = strength, shape=8), size=1.3, color="black") +
  scale_shape_identity()+
  geom_errorbar(mapping = aes(x = names, ymin = lower, ymax = upper), width = 0, color="black", size=.3)+
  geom_point(mapping = aes(names, m), color="black", size=.3) +
  # geom_quasirandom(groupOnX = TRUE, width = .2, aes(color = names)) + # may turn this off when we have 1000 or adjust alpha / flip them so violin is set high and points lower. tbd
  theme_classic()+
  labs(x="Node Names", y="Strength Centrality Estimates") +
  ggtitle("Distribution of Strength Centrality Estimates Across 8-Node Networks") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position ="none",
        panel.grid.major = element_line(color="lightgrey"),
        plot.title = element_text(hjust = 0.5))
name_violin= paste0(resultsdir,'/centdist_violin_valuesebs.png')
ggsave(filename = name_violin, device = "png", plot = cent_dist2, width = 12, height = 5, units = "in")

library(patchwork)
patchworkplot <- cent_dist+cent_dist2+plot_layout(ncol = 1)
name_patch= paste0(resultsdir,'/patchwork.png')
ggsave(name_patch, patchworkplot, width = 15, height = 7)

all_cents[which(all_cents$names=="shame"),]
all_cents[which(all_cents$names=="shame" & all_cents$rank==6),]
# where shame == 6 iter 403, 410, 741
all_cents[which(all_cents$iter=="403"|all_cents$iter=="410"|all_cents$iter=="741"),]

includesshame=subset(df, iter %in% unlist(unique(df[which(df$names=="shame"),]$iter)))
includesshame=includesshame %>% arrange(iter)
table(includesshame$names)

View(df[which(df$names=="drivethin"),] )
# 1st 10.96% 24/219
# 2nd 12.33% 27/219
# 3rd 19.18% 42/219
# 4th 14.61% 32/219
View(df[which(df$names=="alwaysworry"),] )
# 1st 39.22% 91/232
# 2nd 25.43% 59/232
# 3rd 16.38% 38/232
View(df[which(df$names=="fowg"),] )
# 1st 8.23% 19/231
# 2nd 12.12% 28/231
# 3rd 12.12% 28/231
# 4th 10.82% 25/231
# 5th 19.91% 46/231
# 6th 21.65% 50/231
# 7th 10.39% 24/231
# 8th 4.76 11/231

#### z scores ####
# how likely is the observed centrality estimate from the sample distribution?
unique(all_cents$names)
length(unique(all_cents$names)) # 34
apriori_cent
# subset to include items in the a priori network, calculate their mean strength centrality estimates and sds
ztab <- subset(all_cents, names %in% apriori_cent$names) %>% group_by(names) %>% summarise(mean_strength= mean(Strength), sd_strength = sd(Strength), mean_rank=mean(rank), sd_rank=sd(rank), n_count = max(total_n))
ztab <- left_join(ztab, apriori_cent[-c(2:3,5)])
ztab <- ztab %>% rename(observed_cent=strength) %>% rename(observed_rank=rank)
# calculate z scores
ztab$cent_zscore <- (ztab$observed_cent-ztab$mean_strength)/ztab$sd_strength
ztab$rank_zscore <- (ztab$observed_rank-ztab$mean_rank)/ztab$sd_rank

# if observed is lower look at lower tail; if observed is higher than mean, look at upper tail
for(i in 1:nrow(ztab)){
  if(ztab$observed_cent[i] < ztab$mean_strength[i]) {
    ztab$prob_cent[i] <- pnorm(q= ztab$cent_zscore[i], lower.tail = TRUE)
  } else {
    ztab$prob_cent[i] <- pnorm(q= ztab$cent_zscore[i], lower.tail = FALSE)
  }
  if(ztab$observed_rank[i] < ztab$mean_rank[i]){
    ztab$prob_rank[i] <- pnorm(ztab$rank_zscore[i], lower.tail = TRUE)
  } else {
    ztab$prob_rank[i] <- pnorm(ztab$rank_zscore[i], lower.tail = FALSE)
  }
}

ztab$prob_cent <- format.pval(pv = ztab$prob_cent, digits = 2,eps = 0.001, nsmall = 3)
ztab$prob_rank <- format.pval(pv = ztab$prob_rank, digits = 2,eps = 0.001, nsmall = 3)
ztab <- ztab %>% relocate(n_count, .after = names)

write.csv(ztab, paste0(resultsdir, "/ztest_df.csv"), row.names=FALSE)

#### AIM 2 results section ####
length(unique(all_cents$names)) # 34 unique items included in 1000 random networks
test <- as.data.frame(table(all_cents$names))
min(test$Freq) # 205
max(test$Freq) # 259
mean(test$Freq) # 235.2941
sd(test$Freq) # 13.43566
mean(as.data.frame(table(all_cents$names))[,2]) # mean time items were included 235.2941
sd(as.data.frame(table(all_cents$names))[,2]) # sd time items included in networks 13.43566
# mean(ztab$n_count) # 231.375
# sd(ztab$n_count) # 8.227784
min(ztab$n_count) # 219
max(ztab$n_count) # 247

View(ztab[which(ztab$prob_cent < .05 & which(ztab$prob_rank < .05)),])
View(ztab[which(ztab$prob_cent >= .05 & which(ztab$prob_rank >= .05)),])
# get a sense of distribution 
describeBy(all_cents$Strength, all_cents$names) # skew = -0.65 - 0.64
# save for labels
skew_labels <-do.call("rbind", describeBy(all_cents$Strength, all_cents$names))
skew_labels <- cbind(rownames(skew_labels), skew_labels) # make row name (var name) the first column
names(skew_labels)[1] <- "names" # rename this column
skew_labels <- skew_labels[,-2] # get rid of vars column
skew_labels$skew <- round(skew_labels$skew, 2) # round 2 decimal points
skew_labels <- data.frame(names = unlist(skew_labels$names), label = unlist(skew_labels$skew))
# plot distribution
cent_randomdist <-all_cents %>% 
  ggplot(aes(x=Strength)) +
  geom_histogram() +
  facet_wrap(~names) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  geom_text(x = .65, y = 45, aes(label = label), size =2.3, color="red", data = skew_labels) +
  annotate("text", x = .30, y = 45, label= "Skew = ", size =2.3,fontface="italic", color="red")
randist= paste0(resultsdir,'/centdist_randomnets.png')
ggsave(filename = randist, device = "png", plot = cent_randomdist, units = "in")

# qqplot 
all_cents %>% 
  ggplot(aes(sample=Strength)) +
  stat_qq() + stat_qq_line()+
  facet_wrap(~names) +
  theme_classic() 
