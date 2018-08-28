if(Sys.info()['login']=='jeremyb'){
  setwd("C:\\Users\\jeremyb\\Documents\\A_Master\\Code")
}else if(Sys.info()['login']=='jeremy') {
  setwd("C:/Users/JumpCo Vostro3700/infection-dating-tool/manuscripts/figures")
}else{
  setwd(".") #what does this do?
}

library("DescTools")
library("shape")

#Tool for exploring the effect of incorporating more realistic nuances of test `sensitivies' on prospective 

# Initiated March 2018
# Jeremy Bingham



# New section, August 2018

# Individual VIRAL LOAD growth. eats doubling time, initial viral load & vector of times

VL <- function(C0, dblT, times){
  return(C0*2^(times/dblT))
}

# Family of VIRAL LOAD growths. eats matching vectors for doubling times and initial viral loads (equivalent to eclipse lengths)

fvirions <- function(C0s, dblTs, times){
  if(length(C0s != dblTs)){stop("Vector of initial VL must be the same length as vector of doubling times")}
  fam <- matrix(nrow = length(times), ncol = length(C0s))
  for (i in range(length(C0s))){
    fam[,i] <- VL(C0s[i],dblTs[i],times)
  }
  return(fam)
}

# probability of infecting. eats VIRAL LOAD, VOLUME TRANSFUSED, INDIVIDUAL VIRION'S PROBABILITY OF INFECTION.

pinfects <- function(VL,vol,pV){
  return( 1 - exp(-VL * Vol * pV) )
  }

# probability of nondetection in ONE TEST EVENT (ie probability of testing negative). eats VIRAL LOAD (of sample) and TEST'S X50 & X95

pnondetects <- function(VL, X50, X95, Z){
  return(1 - pnorm(Z*(log(VL/X50))/(log(X50/X95)),mean=0,sd=1))
}

# combined likelihood (with RETESTS)
# SAVE FOR LATER
# 
# pnondetTot <- function(VL, X50, X95, Spool, RNAperV, retests, idnat=FALSE){
#   ifelse(idnat,yes = return(  ),no = return(  ) )
# }

# rdays <- function(










#Definitions
# All test functions [unless specified otherwise] are sensitivity functions ie probability of testing positive as a function of time since exposure

# The next section of code defines functions that are used to:
# choose a set of values for one parameter (usually position), with spacing defined using an inverse cumulative normal distribution.
# define and generate individual curves - a sensitivity curve for each type of test, and a description of how to include the population variability to generate a set of curves
# Generate families of n curves for sensitivity   
# Generate families of n curves for likelihood of a certain test result. These functions take a 'test' object which specifies the parameter values
#   (including the set of values for the varying parameter(s) that represent the population-level variability) and the test result ('negative' or 'positive')
#   and a set of times. They return a matrix, each column of which contains the likelihood values for one of the `individual` curves; each curve specifies 
#   the likelihood of the test result given for each hypothetical date of first detectable infection.
# Note: generate_family functions don't return the range, granularity etc of thetime axis - this should be defined somewhere more central in the script.


generate_positions_for_individual_curves = function(n,scale,shape,mean_center_position,sd_size){
  shift_to_half_likelihood = scale*(-1*log(1/2))^(1/shape)
  myseq <- seq(1/(2*n),1-1/(2*n),1/n)
  positions <- qnorm(myseq,mean=mean_center_position-shift_to_half_likelihood,sd=sd_size)
  return(positions)
}

generate_positions_cumulative_normal = function(n,mean_center_position,sd_size){
  myseq <- seq(1/(2*n),1-1/(2*n),1/n)
  positions<-qnorm( myseq,mean=mean_center_position,sd=sd_size )
  return(positions)
}

#a sensitivity is a probability of correctly detecting an infection as a function of time after infection (or fddi doesn't matter here)
#a sensitivity is a probability of a positive result as a function of time since infection
#in our graphs (of likelihood functions) the time since infection is the difference between the time we observe and the test-time. This increases as we move left, which is 
# why we need to reverse the ``direction of time" when we convert from a sensitivity function to a likelihood. We should call the sensitivity function, reversed (comment this clearly so people aren't confused or angry), 
#the probability of testing negative is just 1 -probability of testing positive. no indeterminate results here




# Hello and Welcome

# These comments will guide you gently and clearly throught this code

# if you have any questions don't hesistate to email us at jeremyb@sun.ac.za

#questions I'm wondering about: legends instead of long y-axis labels?

# SECTION ONE: THE SENSITIVITY CURVES

# We begin by defining any and all sensitivity functions we will use. make sure to include an adjustable position parameter
# Each function you define here will have a function of its own in each of the following sections


shift_to_half_likelihood_weibul <- function(scale,shape){
  return(scale*(-1*log(1/2))^(1/shape))
}

sensitivity_weibul <- function(x,scale,shape,position){
  return(ifelse(x-position<0,0,1-exp(-((x-position)/scale)^shape))) #Note that the zero in the ifelse implies that pre-infection test results are never positive. 
  #adjust to incorporate imperfect specificity JEREMY CHECK THIS
}

#SECTION TWO: THE LIKELIHOOD CURVES

# we now write functions to represent each sensitivity function as a likelihood of testing positive/negative, given an infection/detectability time.

# the functional-form name at the end of each function name (eg "weibul") refers to the shape of the SENSITIVITY function for a particular test
# position is determined by the delay and the time of the test (0 for sensitivity) otherwise test_time


individual_sensitivity_weibul <- function(times,scale,shape,delay){
  return(sensitivity_weibul(x=times,scale=scale,shape=shape,position=delay-shift_to_half_likelihood_weibul(shape=shape,scale=scale)))
}

individual_negative_likelihood_weibul  <- function(times, scale, shape, delay, test_time){ #times is a vector of all the times we consider
  return(1-sensitivity_weibul(x=test_time-times,scale=scale,shape=shape,position=delay-shift_to_half_likelihood_weibul(scale=scale,shape=shape)))
}

individual_positive_likelihood_weibul <- function(times, scale, shape, delay, test_time){
  return (sensitivity_weibul(x=test_time-times,scale=scale,shape=shape,position=delay-shift_to_half_likelihood_weibul(scale=scale,shape=shape)))
}

individual_biomarker_weibul <- function(times,scale,shape,delay,height){ # so delays are delays to half-likelihood
  return(biomarker_weibul(x=times,scale=scale,shape=shape,position=delay-shift_to_half_likelihood_weibul(shape=shape,scale=scale),height=height))
}


family_biomarker_weibul<- function(times,set_of_scales,set_of_shapes,set_of_delays, set_of_heights){ #totally generic function that just 
  # creates a set of heights. The order of variable combinations, if any, must be determined via the 
  if(length(set_of_scales)!= length(set_of_shapes)| length(set_of_shapes) != length(set_of_delays) | length(set_of_delays)!= length(set_of_heights)){stop("Some sizes are not lining up!")}
  n = length(set_of_scales)
  set_of_positive_curves <- matrix(nrow = length(times),ncol=n)
  for(i in seq(1,n)){
    set_of_positive_curves[,i] <- individual_biomarker_weibul(times=times,scale=set_of_scales[i],shape=set_of_shapes[i],delay=set_of_delays[i],height = set_of_heights[i])
  } 
  return (set_of_positive_curves)
}


timeaxis <- seq(0,40,1/2)
plotdata_biomarker <- family_biomarker_weibul(timeaxis, set_of_scales = c(1,2,3,4,5,6,7) , set_of_shapes = c(2,3,4,5,6,7,8) , set_of_delays = generate_positions_cumulative_normal(7,10,5) , set_of_heights = seq(0.4,1,.1))

plot(timeaxis,plotdata_biomarker[,1],type = 'n',ylim=c(0,1))
for (i in 1:ncol(plotdata_biomarker)){
  print(i)
  lines(timeaxis,plotdata_biomarker[,i])
}

lines(timeaxis,proportions_recent(timeaxis,plotdata_biomarker,.4),lty=1,lwd=4,col=2)



# We now have code  generate the likelihood curves for our chosen tests

#SECTION THREE: THE FAMILIES

#note that the factor which varies within the family can be any of the parameters or even a combination of parameters
# currently the family-generating functions have variability driven by one input

family_sensitivity_weibul <- function(times, scale, shape, mean_delay, sd_size, n){ #could also list all parameters (with default values) and use if else to select particular sensitivity shape
  list_of_positions <- generate_positions_cumulative_normal(n=n, mean_center_position=mean_delay, sd_size=sd_size)
  set_of_sensitivity_curves <- matrix(nrow=length(times),ncol=n)
  for(i in seq(1,n)){
    set_of_sensitivity_curves[,i] <- individual_sensitivity_weibul(times,scale,shape,list_of_positions[i])
  }
  return(set_of_sensitivity_curves)
}

family_negative_likelihood_weibul <- function(times, scale, shape, mean_delay, sd_size, n, test_time){ #could also list all parameters (with default values) and use if else to select particular sensitivity shape
  list_of_positions <- generate_positions_cumulative_normal(n=n, mean_center_position=mean_delay, sd_size=sd_size)
  set_of_negative_curves <- matrix(nrow=length(times),ncol=n)
  for(i in seq(1,n)){
    set_of_negative_curves[,i] <- individual_negative_likelihood_weibul(times=times,scale=scale,shape=shape,delay = list_of_positions[i],test_time=test_time)
  }
  return(set_of_negative_curves)
}

family_positive_likelihood_weibul <- function(times, scale, shape, mean_delay, sd_size, n, test_time){
  list_of_positions = generate_positions_cumulative_normal(n=n, mean_center_position=mean_delay, sd_size=sd_size)
  set_of_positive_curves <- matrix(nrow=length(times),ncol=n)
  for(i in seq(1,n)){
    set_of_positive_curves[,i] <- individual_positive_likelihood_weibul(times=times,scale=scale,shape=shape,delay = list_of_positions[i],test_time=test_time)
  }
  return(set_of_positive_curves)
}

family_sensitivity_weibul_pos <- function(times, scale, shape, positions, n){ #could also list all parameters (with default values) and use if else to select particular sensitivity shape
  list_of_positions <- positions
  set_of_sensitivity_curves <- matrix(nrow=length(times),ncol=n)
  for(i in seq(1,n)){
    set_of_sensitivity_curves[,i] <- individual_sensitivity_weibul(times,scale,shape,list_of_positions[i])
  }
  return(set_of_sensitivity_curves)
}

family_negative_likelihood_weibul_pos <- function(times, scale, shape, positions, n, test_time){ #could also list all parameters (with default values) and use if else to select particular sensitivity shape
  list_of_positions <- positions
  set_of_negative_curves <- matrix(nrow=length(times),ncol=n)
  for(i in seq(1,n)){
    set_of_negative_curves[,i] <- individual_negative_likelihood_weibul(times=times,scale=scale,shape=shape,delay = list_of_positions[i],test_time=test_time)
  }
  return(set_of_negative_curves)
}

family_positive_likelihood_weibul_pos <- function(times, scale, shape, positions, n, test_time){
  list_of_positions = positions
  set_of_positive_curves <- matrix(nrow=length(times),ncol=n)
  for(i in seq(1,n)){
    set_of_positive_curves[,i] <- individual_positive_likelihood_weibul(times=times,scale=scale,shape=shape,delay = list_of_positions[i],test_time=test_time)
  }
  return(set_of_positive_curves)
}

#calculate the mean of a generic family of curves over a specified set of t-values
generate_mean_of_family = function(family_of_curves){ #input should be a matrix
  return(rowMeans(family_of_curves))
}

#calculate the product of two curves
generate_product_curve = function(curve1,curve2){ #this obviously doesn't need to be its own function I'm just renaming it 
  return(curve1*curve2)
}

#given a curve and find the likelihood for all times the combined test results.
likelihood_by_DDI = function(set_of_positive_curves,set_of_negative_curves,times){
  if(ncol(set_of_negative_curves) != ncol(set_of_positive_curves)) {stop ("Huh?? Why are there different numbers of curves?")}
  if(nrow(set_of_negative_curves) != nrow(set_of_positive_curves)) {stop ("Huh?? Why are there different numbers of time steps?")}
  likelihoods_per_time = rep(0,length(times))
  for (tpos in seq(1:length(times))){
    cumu_likely <- 0
    for (curvenumber in seq(1:ncol(set_of_negative_curves))){ #chronological order makes a trivial difference - there is 
      #just a product per person
      likelihood_forward <- set_of_positive_curves[tpos,curvenumber]*set_of_negative_curves[tpos,curvenumber]
      likelihood_backward <- "The same thing unless I don't actually understand this task"
      cumu_likely <- cumu_likely + likelihood_forward
    }
    likelihoods_per_time[tpos] <- cumu_likely / ncol(set_of_negative_curves)
  }
  return(likelihoods_per_time)
}

#ALSO
#nice plotting function showing likelihood for an individual time point (including dashed lines for unlikely people)

# let's try the same functiont again but generating the curves from scratch each time, or at least input the delay and test_time - is this a good idea or should I just calculate which lines should be dotted based on the actual value (> or < 1/2 at considered time). using the actual value is MUCH less efficient...
# instead, I'll calculate a cutoff curve-number. Ie plot dotted from the cuttoff-th curve
# Okay, I think that may just make it moree difficult than it needs to be....
#goto_now
simple_plot_individual_time_likelihood <- function(n,times,set_of_positive_curves,
                                                   set_of_negative_curves,
                                                   set_of_positive_curves_background,
                                                   set_of_negative_curves_background,
                                                   time,test_of_interest,
                                                   col_negative = col_negative,
                                                   col_positive = col_positive,
                                                   col_likelihood = col_likelihood,
                                                   col_mean = col_mean,
                                                   col_dotted = col_dotted,
                                                   lwd_ind,
                                                   lwd_means,
                                                   lwd_likelihood,
                                                   curve_level_cutoff_probability,
                                                   plot_full_likelihood = FALSE){
  positive_mean_naive = rowMeans(set_of_positive_curves_background)
  negative_mean_naive = rowMeans(set_of_negative_curves_background)
  likelihood_naive <- negative_mean_naive*positive_mean_naive
  #position in timeaxis
  time_position <- which(times==time)
  #calculate_likelihood
  cumu_likely=0
  for (person in seq(1:ncol(set_of_negative_curves_background))){ #chronological order makes a trivial difference - there is 
    #just a product per person
    likelihood_forward = set_of_positive_curves_background[time_position,person]*set_of_negative_curves_background[time_position,person]
    likelihood_backward = "The same thing unless I don't actually understand this task"
    cumu_likely <- cumu_likely + likelihood_forward
  }
  likelihood_DDI_at_time_given_both_test_results = 1/ncol(set_of_negative_curves_background) *cumu_likely
  #plot the curves
  plot(times,plotdata_negative[,1],type='c',xlim=c(times[1],times[length(times)]),ylim=c(-0.002,1.002),xlab='',ylab="",col=col_negative,xaxt='n',yaxt='n',xaxs='i',yaxs='i',bty='l')
  
  # Time of Hypothetical DDI
  # Likelihood of observed test results
  #points(plotdata[,1],plotdata[,3])
  
  # normalising factor
  normalise_negative_test <- 0
  normalise_positive_test <- 0
  for (curve_number in seq(1:ncol(set_of_negative_curves))){
    normalise_negative_test <- normalise_negative_test + set_of_negative_curves[time_position,curve_number]
    normalise_positive_test <- normalise_positive_test + set_of_positive_curves[time_position,curve_number]
  }
  # normalise_negative_test <- normalise_negative_test/n
  # normalise_positive_test <- normalise_positive_test/n
  # print(normalise_negative_test)
  # print(normalise_positive_test)
  
  
  
  lwd_means <- 4
  lwd_ind <- 1.37
  lwd_likelihood <- lwd_means - 3
  
  #col_negative <- "green" #rgb(27/255,158/255,119/255)
  #col_positive <- "red" #rgb(217/255,95/255,2/255)
  #col_mean <- rgb(231/255,41/255,138/255) #"blue" # rgb(231/255,41/255,138/255)
  #col_likelihood <- "purple" #rgb(117/255,112/255,179/255)
  #col_dotted <- rgb(3/7,3/7,3/7) #color for dotted lines
  
  lines(timeaxis,positive_mean_naive,col=col_positive,lwd=lwd_means)
  lines(timeaxis,negative_mean_naive,col=col_negative,lwd=lwd_means)
  
  # col_negative <- 'green'
  # col_positive <- 'red'
  # col_likelihood <- 'purple'
  
  #goto_do
  ## frame the dottedness as posteriors
  ## OR frame is as a direct probability statment
  
  scale_t1= 4
  shape_t1 = 1.73
  
  mean_delay_t1 = 12
  sd_size_t1 = 3
  #   Time of negative test (relative to arbitrary t=0)
  
  test_time_1 = 28
  
  ##                      TEST 2 (positive)
  
  scale_t2 = scale_t1
  shape_t2 = shape_t1
  
  mean_delay_t2 = mean_delay_t1
  sd_size_t2 = sd_size_t1
  
  #   Time of positive test
  test_time_2 = timeaxis[length(timeaxis)]-10
  
  
  ## Generate the individual likelihood curves for the first (negative) and second (positive) test
  #Test 1
  plotdata_negative = family_negative_likelihood_weibul(n=n, scale=scale_t1, shape=shape_t1, mean_delay=mean_delay_t1, sd_size= sd_size_t1, times = timeaxis, test_time = test_time_1)
  #for generating mean curve
  plotdata_negative_background <- family_negative_likelihood_weibul(n=n+50, scale=scale_t1, shape=shape_t1, mean_delay=mean_delay_t1, sd_size= sd_size_t1, times = timeaxis, test_time = test_time_1)
  #Test 2
  plotdata_positive = family_positive_likelihood_weibul(n=n, scale=scale_t2, shape=shape_t2, mean_delay=mean_delay_t2, sd_size= sd_size_t2, times = timeaxis, test_time = test_time_2)
  #for generating mean curve
  plotdata_positive_background <- family_positive_likelihood_weibul(n=n+50, scale=scale_t2, shape=shape_t2, mean_delay=mean_delay_t2, sd_size= sd_size_t2, times = timeaxis, test_time = test_time_2)
  
  
  
  
  #title(xlab="Time", line=1.5, cex.lab=1.2)
  title(ylab=expression("Likelihood"), line=2, cex.lab=1.05)
  #expression(atop("Histogram of "*hat(mu), Bootstrap~samples*','~Allianz))
  
  yaxis_pos <- c(0,1)
  yaxis_names <- c('0','1')
  
  axis(side=2, at=yaxis_pos, labels= yaxis_names,tck=-0.037, padj=.17)
  
  segments(x0=test_time_1,y0=0,x1=test_time_1,y1=1,lty=3)
  segments(x0=test_time_2,y0=0,x1=test_time_2,y1=1,lty=3)
  
  
  if (!plot_full_likelihood){
    if(test_of_interest=="negative"){
      for (curve in seq(1:ncol(set_of_negative_curves))){
        if(!(set_of_negative_curves[time_position,curve]==0)&&!(set_of_positive_curves[time_position,curve]==0)&&(is.nan(set_of_negative_curves[time_position,curve]/normalise_negative_test) | is.nan(set_of_positive_curves[time_position,curve]/normalise_positive_test))){print("NaN Error! at",time,curve)}
        # print(set_of_negative_curves[time_position,curve]/normalise_negative_test)
        if (!(set_of_negative_curves[time_position,curve]==0) && set_of_negative_curves[time_position,curve]/normalise_negative_test>curve_level_cutoff_probability){
          lines(times,set_of_positive_curves[,curve],col=col_positive, lwd=lwd_ind) #what color should the positive curves be 
          lines(times,set_of_negative_curves[,curve],col=col_negative, lwd=lwd_ind) 
        }else{
          lines(times,set_of_positive_curves[,curve],col=col_positive,lty=3,lwd=lwd_ind) #what color should the positive curves be 
          lines(times,set_of_negative_curves[,curve],col=col_negative,lty=3,lwd=lwd_ind)
        }
      }}
    else if(test_of_interest=="positive"){
      for (curve in seq(1:ncol(set_of_negative_curves))){
        if (!(set_of_positive_curves[time_position,curve]==0) &&set_of_positive_curves[time_position,curve]/normalise_positive_test>curve_level_cutoff_probability){
          lines(times,set_of_positive_curves[,curve],col=col_positive,lwd=lwd_ind) #what color should the positive curves be 
          lines(times,set_of_negative_curves[,curve],col=col_negative,lwd=lwd_ind) 
        }else{
          lines(times,set_of_positive_curves[,curve],col=col_positive,lty=3,lwd=lwd_ind) #what color should the positive curves be 
          lines(times,set_of_negative_curves[,curve],col=col_negative,lty=3,lwd=lwd_ind)
        }
      }
    }
    
    points(cex=1.5,pch=19,times[time_position],likelihood_by_DDI(set_of_positive_curves = set_of_positive_curves_background ,set_of_negative_curves = set_of_negative_curves_background,times=times)[time_position])#,lwd=lwd_likelihood,col=col_likelihood)
    
    
    
    # prettify_plot this is messy I know. it makes the gif-generation a bit easier
    
    
    
    xaxis_pos <- c(time,test_time_1,test_time_2)
    xaxis_names <- c(expression('t'['inf']),expression('t'['1']*'(-)'),expression('t'['2']*'(+)'))
    
    axis(side=1, at=xaxis_pos, labels= xaxis_names,padj=-.35,hadj=-.037)
    # axis(side=1, at=zero_pos, labels=zero_name,padj=-0.45,hadj=0.37)
    
    
    #positive_mean_background <- rowMeans(plotdata_positive_background)
    #negative_mean_background <- rowMeans(plotdata_negative_background)
    
    #lines(timeaxis,negative_mean_background, lwd=lwd_means, col=col_negative)      ##taking these out cause I included them in the "simple plot" function
    #lines(timeaxis,positive_mean_background, lwd=lwd_means, col=col_positive)
    
    
    segments(x0=time, y0=0, x1=time,y1=1,lty=2,lwd = 2) }
  
  else{
    
    xaxis_pos <- c(test_time_1,test_time_2)
    xaxis_names <- c(expression('t'['1']*'(-)'),expression('t'['2']*'(+)'))
    
    axis(side=1, at=xaxis_pos, labels= xaxis_names,padj=-.35,hadj=-.037)
    
    for(curve in seq(1:ncol(set_of_negative_curves))){
      lines(times,set_of_positive_curves[,curve],col=col_positive, lwd=lwd_ind) #what color should the positive curves be 
      lines(times,set_of_negative_curves[,curve],col=col_negative, lwd=lwd_ind) 
    }
    
    # WARNING: THIS IS NOT THE TRUE LIKELIHOOD AND WILL BE WRONG FOR OTHER FIGURES
    lines(times,likelihood_naive,lwd=4)
  }
}


plot_individual_time_likelihood = function(n,times,set_of_positive_curves,set_of_negative_curves,time, test_of_interest){ #times is all times, ie timeaxis, time is particular time of interest
  #calculate means
  positive_mean_naive = rowMeans(set_of_positive_curves)
  negative_mean_naive = rowMeans(set_of_negative_curves)
  likelihood_naive <- negative_mean_naive*positive_mean_naive
  #position in timeaxis
  time_position <- which(times==time)
  #calculate_likelihood
  cumu_likely=0
  for (person in seq(1:ncol(set_of_negative_curves))){ #chronological order makes a trivial difference - there is 
    #just a product per person
    likelihood_forward = set_of_positive_curves[time_position,person]*set_of_negative_curves[time_position,person]
    likelihood_backward = "The same thing unless I don't actually understand this task"
    cumu_likely <- cumu_likely + likelihood_forward
  }
  likelihood_DDI_at_time_given_both_test_results = 1/ncol(set_of_negative_curves) *cumu_likely
  #plot the curves
  plot(times,plotdata_negative[,1],type='c',xlim=c(times[1],times[length(times)]),ylim=c(-0.002,1.002),xlab="Hypothetical DDI (days)",ylab="Likelihood of observed test results",col='green')
  #points(plotdata[,1],plotdata[,3])
  if(test_of_interest=="negative"){
    for (curve in seq(1:ncol(set_of_negative_curves))){
      if (!(set_of_negative_curves[time_position,curve]==0) && set_of_negative_curves[time_position,curve]>0.5){
        lines(times,set_of_positive_curves[,curve],col='red', lwd=1.2) #what color should the positive curves be 
        lines(times,set_of_negative_curves[,curve],col='green', lwd=1.2) 
      }else{
        lines(times,set_of_positive_curves[,curve],col='red',lty=2,lwd=1.2) #what color should the positive curves be 
        lines(times,set_of_negative_curves[,curve],col='green',lty=2,lwd=1.2)
      }
    }}
  else if(test_of_interest=="positive"){
    for (curve in seq(1:ncol(set_of_negative_curves))){
      if (!(set_of_negative_curves[time_position,curve]==0) && set_of_positive_curves[time_position,curve]>0.5){
        lines(times,set_of_positive_curves[,curve],col='red',lwd=1.2) #what color should the positive curves be 
        lines(times,set_of_negative_curves[,curve],col='green',lwd=1.2) 
      }else{
        lines(times,set_of_positive_curves[,curve],col='red',lty=2,lwd=1.2) #what color should the positive curves be 
        lines(times,set_of_negative_curves[,curve],col='green',lty=2,lwd=1.2)
      }
    }
  }
  
  lines(times,positive_mean_naive,lwd=2,col='red')
  lines(times,negative_mean_naive,lwd=2,col='green')
  lines(times,likelihood_naive, lwd=4,col='grey') ### todo: transparency? or thin/thick lines? some way to clarify
  points(time,likelihood_DDI_at_time_given_both_test_results)
}


require(DescTools)

find_L_or_E_PDDI <- function(likelihood,times,cutoff){ 
  # get true cutoff not naive threshold. normalise this joint curve (integrate or can I just divide somehow?)
  # 'cutoff' probablity doesn't belong on the y-axis (not based on height based on area under normalised-likelihood pdf)
  total_area <- AUC(times,likelihood,method="spline")
  likelihood_as_probability <- likelihood / total_area
  print(AUC(times,likelihood_as_probability,method="spline"))
  #returns the indices of the EPToi and LPToi
  infection_window <- c(0,times[length(times)])
  found<-0
  for (time in  seq(1,length(times))){
    if(DescTools::AUC(times[1:time],likelihood_as_probability[1:time],method="spline")>cutoff/2 && found==0){
      infection_window[1] <- time #need to avoid this check once it's been found
      found <- 1
    }
    else if(AUC(times[1:time],likelihood_as_probability[1:time],method="spline")>1-cutoff/2 && found==1){
      infection_window[2] <- time
      found <- 2
    }
  }
  return(infection_window)
}




#
# 
# calculate_cuttoff_from_negative_curves <- function(time,timeaxis,set_of_negative_curves){
#  # given a hypothetical DDI returns the index of the first curve in the set  
#   # for which the likelihood of a negative test is smaller than 1/2
#    step <- which.min(abs(timeaxis-time))
#   nearest_curve <- which.min(abs(set_of_negative_curves[step,]-0.5))
#   if(set_of_negative_curves[step,nearest_curve]<0.5){
#     return(nearest_curve+1)}
#   else{return (nearest_curve)
#   }
# } 

# gotoscript

###     Some notes on parameters

#Scale defines interval from the last time where sample almost definitely tests negative to the first time when sample almost definitely tests positive testing positive. 
#Shape #defines the shape of the sensitivity curve as it increases from around 0 to around 1. Also affects the amount of time this takes, though not as much as scale.

