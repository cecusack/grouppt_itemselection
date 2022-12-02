load("~/Dropbox/PT group network paper (PT open series)/grouppt_dev/R15 EMA post hoc/R15_posthoc.RData")
rm(all_cents, btwn, btwn_cent, data_lst, fit_r15, freq, pcc, pdc, 
   pdc_cent, pdc_imp, test, zscore_r15, badids, common_path, i, nodes, top8)

#### subset keep treatment targets that match PT open series ####
library(reshape2)
library(tidyverse)
library(mlVAR)
library(qgraph)
library(plyr)
library(ggbeeswarm)
library(clairecleans)
library(here)
here()

r15_subset = dat %>% 
  select(ID, dayvar, beepvar, shame, guilty, overvalwtsh, worryoverwhelm, 
         fearloc, mistakes, emo_overwhelm, ineffective, bodydiss, rumination, 
         intrus_thought, reject, mistakes, iuc, relax, eat_anx, food_intrus, 
         desirethin, saa, fowg, avoid_food, food_intrus, feelfat, hunger_anx, 
         physsens, bodycheck, foodrules, exercise, highstand, sens_body, avoid_emo, 
         urge_restrict, binge, eat_public)


#### identify means ####
clairecleans::item_sel(r15_subset[-c(1:3)], 8)
# bodydiss, feel fat, fowg, urge_restrict, mistakes, desirethin, overvalwtsh, highstand

#### combinations ####
# column numbers for target items
items <-  seq(4,35,by =1) # items are in cols 9 thru 42
# all possible combinations of 8 items
# all_combinationr15 <- combn(x = items, m = 8)
dim(all_combinationr15) # whew 8 by 10,518,300

# write.csv(all_combinationr15,file = 'node8_all32_r15.csv',row.names = F)

# all_10000_r15 <- all_combinationr15[,sample(ncol(all_combinationr15),10000)]
# dim(all_10000) # 8 x 10000

# write.csv(all_10000_r15,file = 'node8_10000_r15.csv',row.names = F)

# select a subset of 1,000 8-node networks 
# all_1000_r15 <- all_combinationr15[,sample(ncol(all_combinationr15),1000)]
# dim(all_1000) # 8 x 1000

# write.csv(all_1000_r15,file = 'node8_1000_r15.csv',row.names = F)
# all_combination <- all_1000_r15
# rm(all_10000_r15, all_combinationr15, all_1000_r15)
all_combination <- read.csv("node8_1000_r15.csv")

#### set up ####
dim(r15_subset) # 4575 x 35
round(sum(is.na(r15_subset[,-c(1:3)]))/prod(dim(r15_subset[,-c(1:3)]))*100,2) # 31.3%
length(unique(r15_subset$ID)) # 41
boot <- dim(all_combination)[2]
# set empty array to store all cent values
all_cents <- array(NA,dim = c(4,8,boot)) # create an empty 3D array
# set progress bar
pb <- txtProgressBar(min=0, max=boot, style=3,char="=")
init <- numeric(boot)
end <- numeric(boot)

#### estimate 1000 nets ####
for(j in 1:boot){
  
  print(j) # so we know which iteration we're on
  init[j] <- Sys.time()
  network_variance <- sd(melt(r15_subset[,all_combination[,j]],id.vars=NULL)$value)
  # r15_subset[,all_combination[,j]] <- as.data.frame(lapply(r15_subset[,all_combination[,j]], sample))
  # set empty matrix with 4 rows (item names, instrength, outstrength, and strength) and 8 columns (centrality estimates for 8 nodes)
  network_1 <- matrix(0, nrow=4, ncol=8,
                      dimnames = list(c("nodes", "InStrength", "OutStrength", "Strength"),
                                      colnames(r15_subset)[all_combination[,j]]))
  
  network_1[1,] <- colnames(r15_subset)[all_combination[,j]] # make first row items selected for permutation
  
  # var_group <- mlVAR(r15_subset[,c(which(colnames(r15_subset)=="ID"), all_combination[,j])], vars = colnames(r15_subset[,all_combination[,j]]), idvar = "ID", lags = 1)
  
  var_group <- mlVAR(r15_subset[,c(which(colnames(r15_subset)=="ID"), which(colnames(r15_subset)=="dayvar"):which(colnames(r15_subset)=="beepvar"), all_combination[,j])],  vars= colnames(r15_subset[,all_combination[,j]]), idvar="ID", beepvar="beepvar", dayvar="dayvar", lags = 1,temporal = "correlated", estimator = "lmer", contemporaneous = "correlated", verbose = TRUE, scale = TRUE, scaleWithin = FALSE, AR = FALSE, MplusSave = TRUE, MplusName = "mlVAR", iterations = "(2000)", chains = nCores)
  
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
resultsdir <- here("/Users/clairecusack/Dropbox/PT group network paper (PT open series)/grouppt_dev/R15 EMA post hoc/RandomNets")
# dir.create(resultsdir)
name_cent= paste0(resultsdir,'/all_cents_r15')
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
all_cents$iter <- sub("^(\\d+).*", "\\1", all_cents$iter) # keep only the combination number

name_centdf= paste0(resultsdir,'/all_centsDF')
filename2= paste(name_centdf, Sys.Date(), filetype, sep='')
saveRDS(all_cents, filename2)

all_cents$names<-as.factor(all_cents$names) # char to factor
all_cents %>% group_by(names) %>% filter(n()<2) # which items were only included once 0

#### plot distribution ####
# cent from study 1
df <- readRDS("~/Dropbox/PT group network paper (PT open series)/grouppt_dev/R15 EMA post hoc/RandomNets/all_centsDF2022-11-07.rds")
all_cents_study1 <- readRDS("~/Dropbox/PT group network paper (PT open series)/grouppt_dev/Group_cent/all_centsDF2022-06-19.rds")
View(all_cents_study1)
study1_strength <- all_cents_study1 %>% 
  group_by(names) %>% 
  summarise(mean=mean(Strength), rank=mean(rank))

setdiff(study1_strength$names, all_cents$names)
#  [1] "alwaysworry" "avoidfood"  [3] "cogrest"     "drivethin"  
# [5] "eatanx"      "eatpublic" [7] "guilt"       "hungeranx"  
# [9] "interocept"  "obsess"     [11] "overval"     "pastmeal"   
# [13] "pays_sens"   "ruminate"   [15] "selfcrit"   
all_cents$names <- mgsub::mgsub(all_cents$names, c( "worryoverwhelm", "avoid_food", "urge_restrict", "desirethin",
                           "eat_anx", "eat_public", "guilty", "hunger_anx",
                           "sens_body", "intrus_thought", "overvalwtsh",
                           "physsens", "rumination"),
             c("alwaysworry", "avoidfood", "cogrest", "drivethin",
               "eatanx", "eatpublic", "guilt", "hungeranx", 
               "interocept", "obsess", "overval", "pays_sens", "ruminate"))
setdiff(all_cents$names, study1_strength$names)
setdiff(study1_strength$names, all_cents$names) # "past_meal" "selfcrit"  

cent_dist <- 
  all_cents %>% 
  group_by(names) %>% filter(n()>1) %>% ungroup() %>% # remove items only included once
  ggplot(aes(x=reorder(names,rank, mean),y=rank, color=names,fill=names)) +
  geom_violin(alpha=0.7) +
  # geom_quasirandom(groupOnX = TRUE, width = .2, aes(color = names)) + # may turn this off when we have 1000 or adjust alpha / flip them so violin is set high and points lower. tbd
  
  geom_point(data = study1_strength %>% filter(names %in% all_cents$names), 
             aes(x = names, y = rank, shape=8), size=1.3, color="black")+
  scale_shape_identity()+
  theme_classic()+
  labs(x="", y="Strength Centrality Rank Order") +
  ggtitle("Distribution of Strength Centrality Rank Order Across Study 1 and Study 2") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position ="none",
        panel.grid.major = element_line(color="lightgrey"),
        plot.title = element_text(hjust = 0.5)) # ,
#         panel.grid.major.y = element_blank())

cent_dist # see it
name_violin= paste0(resultsdir,'/centrank_violin_nojitter_r15.png')
ggsave(filename = name_violin, device = "png", plot = cent_dist, width = 12, height = 5, units = "in")

#### plot strength centrality not just ranks ####
#### summarize ####
freq <- as.data.frame(table(all_cents$names))
names(freq) <- c("names", "total_n")
# View(all_cents[which(all_cents$names=="shame"),])
df <- all_cents %>% group_by(names, rank) %>% tally()
df <- left_join(df, freq) # add total n
# rm(freq)
df$perc <- round((df$n/df$total_n)*100,2)
# df <- left_join(df, all_cents[-c(2:3)])
# save
name_sum=paste0(resultsdir,'/summaryranks_r15')
filetype= '.csv'
filename= paste(name_sum, Sys.Date(), filetype, sep='')
write.csv(df,file=filename, row.names = FALSE)


cis <- all_cents %>%
  group_by(names) %>%
  summarise_at(vars(Strength), list(m = mean, sd=sd)) 
total_n <- df %>% group_by(names) %>% select(names, total_n) %>% unique
cis <- left_join(cis, total_n)
cis <- left_join(cis, freq)
cis <- cis %>% mutate(upper = m + 1.96*(sd/sqrt(total_n)), lower = m - 1.96*(sd/sqrt(total_n)))
all_cents <- left_join(all_cents, cis)
setdiff(all_cents_study1$names, all_cents$names) # "pastmeal" "selfcrit"
setdiff(all_cents$names,all_cents_study1$names) # 0

cent_dist2 <- all_cents %>% 
 #  group_by(names) %>% filter(n()>1) %>% ungroup() %>% # remove items only included once
  ggplot(aes(x=reorder(names,Strength, mean, decreasing = TRUE),y=Strength, color=names, fill=names)) +
  geom_violin(alpha=0.7) +
  # geom_point(stat="summary", fun.y="mean", color="black", size=.5)+
  # geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), color="black", width=0) +
  # stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
  #              colour = "black")+
  geom_point(data = study1_strength %>% filter(names %in% all_cents$names), 
             aes(x = names, y = mean, shape=8), size=1.3, color="black") +
  scale_shape_identity()+
  geom_errorbar(mapping = aes(x = names, ymin = lower, ymax = upper), width = 0, color="black", size=.3)+
  geom_point(mapping = aes(names, m), color="black", size=.3) +
  # geom_quasirandom(groupOnX = TRUE, width = .2, aes(color = names)) + # may turn this off when we have 1000 or adjust alpha / flip them so violin is set high and points lower. tbd
  theme_classic()+
  labs(x="Node Names", y="Strength Centrality Estimates") +
  ggtitle("Distribution of Strength Centrality Estimates Across Study 1 and Study 2") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position ="none",
        panel.grid.major = element_line(color="lightgrey"),
        plot.title = element_text(hjust = 0.5))
name_violin= paste0(resultsdir,'/centdist_violin_valuesebs_r15.png')
ggsave(filename = name_violin, device = "png", plot = cent_dist2, width = 12, height = 5, units = "in")

library(patchwork)
patchworkplot <- cent_dist+cent_dist2+plot_layout(ncol = 1)
name_patch= paste0(resultsdir,'/patchwork_r15.png')
ggsave(name_patch, patchworkplot, width = 15, height = 7)

# cat(paste(shQuote(unlist(unique(dat$ID)), type="cmd"), collapse=", "))


#### z scores ####
# how likely is the observed centrality estimate from the sample distribution?
unique(all_cents$names)
length(unique(all_cents$names)) # 32
# subset to include items in the a priori network, calculate their mean strength centrality estimates and sds
ztab <- subset(all_cents, names %in% study1_strength$names) %>% group_by(names) %>% summarise(mean_strength= mean(Strength), sd_strength = sd(Strength), mean_rank=mean(rank), sd_rank=sd(rank), n_count = max(total_n))
ztab <- left_join(ztab, study1_strength)
ztab <- ztab %>% rename(study1_meancent=mean) %>% rename(study1_meanrank=rank)
# calculate z scores
ztab$cent_zscore <- (ztab$mean_strength-ztab$study1_meancent)/ztab$sd_strength
ztab$rank_zscore <- (ztab$mean_rank-ztab$study1_meanrank)/ztab$sd_rank

# if observed is lower look at lower tail; if observed is higher than mean, look at upper tail
for(i in 1:nrow(ztab)){
  if(ztab$mean_strength[i] < ztab$study1_meancent[i]) {
    ztab$prob_cent[i] <- pnorm(q= ztab$cent_zscore[i], lower.tail = TRUE)
  } else {
    ztab$prob_cent[i] <- pnorm(q= ztab$cent_zscore[i], lower.tail = FALSE)
  }
  if(ztab$mean_strength[i] < ztab$study1_meanrank[i]){
    ztab$prob_rank[i] <- pnorm(ztab$rank_zscore[i], lower.tail = TRUE)
  } else {
    ztab$prob_rank[i] <- pnorm(ztab$rank_zscore[i], lower.tail = FALSE)
  }
}

ztab$prob_cent <- format.pval(pv = ztab$prob_cent, digits = 2,eps = 0.001, nsmall = 3)
ztab$prob_rank <- format.pval(pv = ztab$prob_rank, digits = 2,eps = 0.001, nsmall = 3)
ztab <- ztab %>% relocate(n_count, .after = names)

write.csv(ztab, paste0(resultsdir, "/ztestr15_df.csv"), row.names=FALSE)

#### study 2 results section ####
length(unique(all_cents$names)) # 32 unique items included in 1000 random networks
test <- as.data.frame(table(all_cents$names))
min(test$Freq) # 220
max(test$Freq) # 277
mean(test$Freq) # 250
sd(test$Freq) # 14.55357
mean(as.data.frame(table(all_cents$names))[,2]) # mean time items were included 250.00
sd(as.data.frame(table(all_cents$names))[,2]) # sd time items included in networks 14.55357
# mean(ztab$n_count) # 250
# sd(ztab$n_count) # 14.55357
min(ztab$n_count) # 220
max(ztab$n_count) # 277

View(ztab[which(ztab$prob_cent < .05),])
# cog restraint and obsess were outside expected range
# both higher in r15
View(ztab[which(ztab$prob_rank < .05),])
# rank order expected
View(ztab[which(ztab$prob_cent < .05 & which(ztab$prob_rank < .05)),])
# 0 items had centrality estimates were outside the range expected from study 1
View(ztab[which(ztab$prob_cent >= .05 & which(ztab$prob_rank >= .05)),])
# 30 items had both cent and rank expected based on study 1
choose <- c("obsess", "shame", "bodydiss", "mistakes", "guilt", "fearloc",
            "alwaysworry",  "eatanx")

View(ztab %>% filter(names %in% choose) %>% 
       select(names, n_count, mean_strength, sd_strength, study1_meancent, cent_zscore, prob_cent))

View(ztab %>% filter(names %in% choose) %>% 
       select(names, n_count, mean_rank, sd_rank, study1_meanrank, rank_zscore, prob_rank))




