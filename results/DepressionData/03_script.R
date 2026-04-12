################## Analyses reported in ######################
################## Elmer & Stadtfeld  #######################
#### Depressive symptoms are associated with social #########
#### isolation in face-to-face interaction networks #########
########### submitted to Scientific Reports ################


rm(list = ls()) # empty the environment

# load the data ###
load("Data.RData")
# the age variable was replaced with random data to prevent participants from being identified. 
str(dat)

## correlation table (Table 1 of the manuscript)
corstarsl(dat,which.spearman = 3)


#### MRQAP analyses: testing depression-isolation, depression-homphily, and depression-friendship hypotheses
## define the multi-group MRQAP function ##

QAP.MG <- function(dvs, ivs,  iv.names = iv.names, mode = "yQAP" ,samples = 1000, diag = F, directed = T, 
                   global.deltas = T, ecdf.plot = F, return.perms = F){
  pb <- txtProgressBar(min = 0, max = samples, style = 3) # set progress bar
  
  
  if(!diag){ # get rid of diagonal values
    for(DV in 1:length(dvs)){
      diag(dvs[[DV]]) <- NA
    }
    
    for(IV in 1:length(ivs)){  
      IV.dat <- ivs[[IV]]
      for(subIV in 1:length(ivs[[IV]])){
        diag(ivs[[IV]][[subIV]]) <- NA
      }
    }
  }
  
  
  
  # get the observed estimates
  
  if(directed){
    observedLm <- lm(unlist(dvs) ~ Reduce(rbind,lapply(1:length(ivs),
                                                       function(x) sapply(ivs[[x]], function(y) y)
    )))
    
  }else{
    
    observedLm <- lm(unlist(lapply(dvs, function(x) x[lower.tri(x, diag = diag)])) ~
                       Reduce(rbind,lapply(1:length(ivs),function(x) sapply(ivs[[x]], function(y) y[lower.tri(y, diag = diag)]))))
    
  }
  
  
  if(mode == "linearregression"){
    names(observedLm$coefficients) <- iv.names
    return(observedLm) 
  }
  
  # prepare output file
  observedEstimates <- observedLm$coefficients
  output <- data.frame(Estimates = observedEstimates)
  rownames(output) <- iv.names
  
  r.squared <- c(summary(observedLm)$r.squared, summary(observedLm)$adj.r.squared)
  names(r.squared) <- c("r.squared","adj.r.squared")
  
  
  
  if(mode == "yQAP"){
    sampledEstimates <- data.frame()
    for(sampleNr in 1:samples){
      if(directed){
        sampledEstimates <- rbind(sampledEstimates, lm(unlist(lapply(dvs, function(x) sample(x))) ~ Reduce(rbind,lapply(1:length(ivs),
                                                                                                                        function(x) sapply(ivs[[x]], function(y) y)
        )))$coefficients)
      }else{
        sampledEstimates <- rbind(sampledEstimates, lm(unlist(lapply(dvs, function(x) sample(x[lower.tri(x, diag = diag)]))) ~ 
                                                         Reduce(rbind,lapply(1:length(ivs), function(x) sapply(ivs[[x]], function(y) y[lower.tri(y, diag = diag)])
                                                         )))$coefficients)
      }
      
      setTxtProgressBar(pb, sampleNr)
    }
    
    ecdf.plots <- list()
    for(est in 1:length(observedEstimates)){

      output[est,"p(1sided)"] <- ecdf(sampledEstimates[,est])(observedEstimates[est])
      output[est,"abs(p)"] <- output[est,"p(1sided)"] 
      output[est,"abs(p)"][output[est,"p(1sided)"] > 0.5] <- 1-output[est,"abs(p)"][output[est,"p(1sided)"] > 0.5]
      output[est,"adj.d"] <- abs(mean(sampledEstimates[,est]) - observedEstimates[est])/sd(sampledEstimates[,est])
      output[est,"Exp.V"] <- mean(sampledEstimates[,est], na.rm = T)
      output[est,"Exp.V.sd"] <- sd(sampledEstimates[,est], na.rm = T)
      output[est,"2.5th P"] <- quantile(sampledEstimates[,est],probs=c(.025))
      output[est,"97.5th P"] <- quantile(sampledEstimates[,est],probs=c(.975))
     
      
      if(ecdf.plot){
        
        require(ggplot2)

        ecdf.plots[[est]] <- plot(ecdf(sampledEstimates[,est])(observedEstimates[est]))
        
      }
    }
    
  }
  
  
  output <- round(output,5)
  
  stars <- function(p.value){
    ifelse(p.value > 0.90 | p.value < 0.10,
           ifelse(p.value > 0.95 | p.value < 0.05, 
                  ifelse(p.value > 0.99 | p.value < 0.01, 
                         ifelse(p.value > 0.999 | p.value < 0.001, "***", "**"), "*"), "x"), "")
  }
  output[,"significance"] <- sapply(output$`p(1sided)`, stars)
  if(return.perms) return(list(sampledEstimates, observedEstimates))
  return(list(mode = c(mode, samples), plots = ecdf.plots,output = output, r.squared = r.squared))
  
}

### estimate the Multi-group QAP model reported in Table 2

dvs[[1]] # interaction matrix of sample 1 (in seconds per hour)
dvs[[2]] # interaction matrix of sample 2 (in seconds per hour)

# independent variables of sample 1 (ivs[[1]]), sample 2 can be accesed with ivs[[2]]
ivs[[1]][[1]] # dummy matrix for sample 2 
ivs[[1]][[2]] # dummy matrix at least one female 
ivs[[1]][[3]] # dummy matrix both female
ivs[[1]][[4]] # age mean (centered)
ivs[[1]][[5]] # age similarity # random data, for anonymity reasons
ivs[[1]][[6]] # one student organization # random data, for anonymity reasons
ivs[[1]][[7]] # same student status (being in a student organization)
ivs[[1]][[8]] # being friends 
ivs[[1]][[9]] # depression mean -> testing the depression-isolation hypothesis
ivs[[1]][[10]] # depression similarity -> testing the depression-homphily hypothesis
ivs[[1]][[11]] # depression mean * depression similarity  
ivs[[1]][[12]] # depression mean * being friends -> testing the depression-friendship hypothesis


model.fit <- QAP.MG(dvs = dvs, # 2 list of one matrix each containing the dependent networks (sample 1, sample 2)
                    ivs = ivs, # 2 list of multiple matrices each containing the independent variables (sample 1, sample 2)
                   iv.names =  iv.names, # names of the independent variables
                   directed = F, # we have an unidrected dependent nectwork
                   samples = 5000, # number of permuted matrices
                   mode = "yQAP")

model.fit # results
# results differ slightly from those reported in the manuscript due to the removal of the age variable (for anonymity reasons)

### testing the dyadic-isolation hypothesis ####

psych::pairs.panels(dat[,c("depression.L1","ratio")])
cor.test(dat$depression.L1,dat$ratio)
estimate <- cor.test(dat$depression.L1,dat$ratio)$estimate 

## a y-permutation test:
corr.coeffs <- list()
for(perm in 1:5000){
  corr.coeffs[[perm]] <- cor.test(y = sample(dat$ratio), x = dat$depression.L1)$estimate
}

# compare
output.cor <- data.frame(Estimates = estimate)

output.cor[,"p(1sided)"] <- 1-ecdf(t(sapply(corr.coeffs, function(x) x)))(estimate)

output.cor #result described in the manuscript

### testing the QAP.MG function ####

library(sna)
# create test data #
# inspired by the example funciton in ?netlm
ivnet1<-rgraph(20,4)
ivnet2<-rgraph(20,4)

dv1<-ivnet1[1,,]+4*ivnet1[2,,]+2*ivnet1[3,,]   # Note that the fourth graph is unrelated
dv1 <- dv1 + rnorm(400,mean = 1, sd = 1)

dv2 <- 2*ivnet2[1,,]+3*ivnet2[2,,]+3*ivnet2[3,,]
dv2 <- dv2 + rnorm(400,mean = 1, sd = 1)
dvs <- list(dv1, dv2)

iv1 <- list(ivnet1[1,,],ivnet1[2,,],ivnet1[3,,], ivnet1[4,,])
iv2 <- list(ivnet2[1,,],ivnet2[2,,],ivnet2[3,,], ivnet2[4,,])
ivs <- list(iv1, iv2)

# both groups together
QAP.MG(dvs, ivs, iv.names = c("intercept",paste0("IV",1:4)))

# group 1 separately
QAP.MG(list(dv1), list(iv1), iv.names = c("intercept",paste0("IV",1:4)), samples = 3000)

# comparison with the netlm function from the sna-package
netlm(dv1, iv1, nullhyp = "qapy", reps = 3000)


