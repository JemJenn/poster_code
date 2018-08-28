# August 2018
# Jeremy Bingham
# Script developing on top of concepts from Weusten et al 2011 paper "refinements on a residual risk model for donated blood products"

if(Sys.info()['login']=='jeremyb'){
  setwd("C:\\Users\\jeremyb\\Documents\\A_Master\\Code")
}else if(Sys.info()['login']=='jeremy') {
  setwd("C:/Users/JumpCo Vostro3700/infection-dating-tool/manuscripts/figures")
}else{
  setwd(".") #what does this do?
}

source("part_Onefunctions.R")
source("parameters.R")

file.create("C:\\Users\\jeremyb\\Documents\\A_Master\\Code\\runlog_poster.txt")
runlog <- file("runlog.txt")
writeLines(capture.output(sessionInfo()),runlog)
close(runlog)


#Question list
# don't we just want total days at risk? for that, given an individual, doesn't it just depend on viral growth rate, not timing?
# buuuut maybe the point is that if we only consider the mean curve we miss out on the nuance that the distribution adds...

set.seed = 921002

#Master equation
# overallRisk <- rdays/tbetween * donors_converting/donations_total

# tbetween is a number... we just use the average. since tbetween (FOR SEROCONVERTING DONORS) matters just for the total time
# and we compare it to the total risk days. so for each converting donor, how much time (on average) do they have in which
# they COULD have been infected, and for what `fraction' (risk days) of that time would they pose a risk of infection?
#  so we begin by calculating risk days

# Viral loads! We begin by calculating a distribution of doubling times, using an inverse cumulative normal distribution to
# approximate the average that would be obtained by sampling many times from a normal distribution with the mean and 
# standard deviation garnered from Weusten et al.

detail <- 10
timeaxis <- seq(0,50,1/detail)

n <- 15
nDisp <- 100


dblTs <- generate_positions_cumulative_normal(n=n, mean_center_position = dblT_mean, sd_size = dblT_sd)
C0s <- seq(.1,.1,length.out=n)#generate_positions_cumulative_normal(n=n, mean_center_position = 0, sd_size = 10) #arbitrary parameters. This shouldn't matter in any case...

eclipses <- generate_positions_cumulative_normal(n=n, mean_center_position = 20, sd_size = 2)

# eclipses <- sample(generate_positions_cumulative_normal(n=n, mean_center_position = 20, sd_size = 2),n,replace=FALSE)

# limitation: we don't know how pV (probability of one virion causing an infection, which varies a lot) is correlated (or not) with 
# the doubling time or the eclipse phase (maybe the latter matters not but the former surely does!)

#generate a representative family of viral loads, from some arbitrary base time point....

VLs <- fvirions(C0s, dblTs, eclipses, timeaxis)

family_infects <- finfects(VLs=VLs, vol=vol, pV = pV_mean, times = timeaxis)

family_nondetects <- fnondetects(VLs = VLs, X50=X50, X95=X95, Z=Z, times = timeaxis)

summary(family_nondetects)

# rdays <- AUC

#Probability of nondetection and infection occuring =
plot(timeaxis,family_infects[,1],type='l')
for (i in seq(1:n)){
  lines(timeaxis,family_infects[,i])
  lines(timeaxis,family_nondetects[,i],col='blue')
}

# pr_ndi <- integrate()


# Individual VIRAL LOAD growth. eats doubling time, initial viral load & vector of times

VL <- function(C0, dblT, eclipse, times){
  return(C0*2^((times-eclipse)/dblT))
}

# Family of VIRAL LOAD growths. eats matching vectors for doubling times and initial viral loads (equivalent to eclipse lengths)

fvirions <- function(C0s, dblTs, eclipses, times){
  if(length(C0s) != length(dblTs)){stop("Vector of initial VL must be the same length as vector of doubling times")}
  fam <- matrix(nrow = length(times), ncol = length(C0s))
  for (i in seq(1:length(C0s))){
     fam[,i] <- VL(C0 = C0s[i],dblT = dblTs[i],eclipses[i], times = times) #could also use sample(eclipses,1)
  }
  return(fam)
}

finfects <- function(VLs, vol, pV, times){
  fam <- matrix(nrow = length(times), ncol = length(VLs))
  for (i in seq(1:ncol(VLs))){
    fam[,i] <- pinfects(VLs[,i], vol, pV)
  }
  return(fam)
}

# probability of infecting. eats VIRAL LOAD, VOLUME TRANSFUSED, INDIVIDUAL VIRION'S PROBABILITY OF INFECTION.

pinfects <- function(VL,vol,pV){
  return( 1 - exp(-VL * vol * pV) )
}

# probability of nondetection in ONE TEST EVENT (ie probability of testing negative). eats VIRAL LOAD (of sample) and TEST'S X50 & X95

pnondetects <- function(VL, X50, X95, Z){
  return(1 - pnorm(Z*(log(VL/X50))/(log(X95/X50)),mean=0,sd=1))
}


fnondetects <- function(VLs, X50, X95, Z, times){
  fam <- matrix(nrow = length(times), ncol = length(VLs))
  for (i in seq(1:ncol(VLs))){
    fam[,i] <- pnondetects(VL = VLs[,i],X50 = X50, X95 = X95, Z = Z)
  }
  return(fam)
}
