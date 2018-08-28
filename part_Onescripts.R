if(Sys.info()['login']=='jeremyb'){
  setwd("C:\\Users\\jeremyb\\Documents\\A_Master\\Code")
  }else if(Sys.info()['login']=='jeremy') {
  setwd("C:/Users/JumpCo Vostro3700/infection-dating-tool/manuscripts/figures")
  }else{
  setwd(".") #what does this do?
  }

source("part_Onefunctions.R")

file.create("C:\\Users\\jeremyb\\Documents\\A_Master\\Code\\runlog_main.txt")
runlog <- file("runlog.txt")
writeLines(capture.output(sessionInfo()),runlog)
close(runlog)
?writeLines
?date

##                    SCRIPT
# Now we are going towards risk estimation
# So we define the quick transition from non-infectious to infectious, and the probability of detection (the 'sensitivity' of the test), both with Weibul distributions

n = 5

detail <- 10
timeaxis <- seq(0,90,1/detail)

shape_infect <- 1
scale_infect <- 5
mean_delay_infect <- 7
population_sd_infect <- 10

shape_detect <- 5
scale_detect <- 5
mean_delay_detect <- mean_delay_infect + 14
population_sd_detect <- population_sd_infect

donation_time <- 77

# now we generate and plot the curves

infectious_delays <- generate_positions_cumulative_normal(n=n,mean_center_position = mean_delay_infect, sd_size = population_sd_infect)
detectable_delays <- infectious_delays + 14

infectious_curves <- family_positive_likelihood_weibul_pos(times = timeaxis, shape = shape_infect, scale = scale_infect, positions = infectious_delays, n = n, test_time = donation_time)
undetected_curves <- family_negative_likelihood_weibul_pos(times = timeaxis, shape = shape_detect, scale = scale_detect, positions = detectable_delays, n = n, test_time = donation_time)





# I want a function to take a set of infectious curves and undetected curves and do the three-layered analysis we're proposing
# or it could take the parameters required for all three calculations, and do them.
# then we can take sets of parameters gleaned from literature and stick them together...

#but now how do I calculate the background ones....... if I decide to shift the likelihood by some amount? I'll probably have to define a function or a block of 
# code which shifts the infectious delays to get corresponding detectable delays, then just run that block for a much bigger infectiousness set.. no problem



# need to think about whether there is any difference between defining that everyone waits two weeks from infectious to detectable,but the total delay is distributed, and defining that the infectiousness
# infectiousness aand detectability as distributed, using the same distribution and different means (delta is two weeks). No difference in the non-mixed case, but say we want to 
# allow the time between infectiousness and detectability to vary. In that case it is more transparent to define a bunch of infectiousness times and delay-between-infectiousness-and-detectability-times,
# then calculate detectability times accordingly. That way we can change the map from 1 - 1 to 1 - many.... Defining them conditionally seeeems to make more sense but I'm not sure.

# anyway that shouldn't take me too long, so I can start with that tonight. #DONEZOES. Now I've specified the detectability delays ACCORDING to the 
# next on the list is to check the old figures with some mixing included. (imperfect correlation). ALSO, check that zero correlation is what they say it is....



#  Figure 1
################
# : produce a family of sensitivity curves and their average
#goto_1

#Parameters

n=10

# test_details
mean_delay_t1=55
mean_delay_t2=45       #we have a second test listed so that we can easily swap between them when generating the sensitivity curves

sd_size_t1 = 5
sd_size_t2 = 1
detail = 10 # 1/Step size
timeaxis = seq(0,100,1/detail)

scale_t1= 3
shape_t1 = 2



scale_t2= 3
shape_t2 = 2


#visuals
# col_negative <- rgb(27/255,158/255,119/255)
# col_positive <- rgb(217/255,95/255,2/255)
# col_mean <- rgb(231/255,41/255,138/255)
# col_truth <- rgb(117/255,112/255,179/255)
# col_likelihood <- col_truth
# col_dotted <- rgb(3/7,3/7,3/7) #color for dotted lines

#col_negative <- "green" #rgb(27/255,158/255,119/255)
#col_positive <- "red" #rgb(217/255,95/255,2/255)
col_negative <- rgb(27/255,158/255,119/255)
col_positive <- rgb(217/255,95/255,2/255)
col_mean <- rgb(231/255,41/255,138/255)
#col_truth <- "purple"  #rgb(117/255,112/255,179/255)
col_truth <- rgb(117/255,112/255,179/255)
col_likelihood <- col_truth
col_dotted <- rgb(3/7,3/7,3/7) #color for dotted lines

#y-axis limits
ylim_chosen <- c(-.002,1.007)

#Generate Data
sensitivity_family_1 <-  family_sensitivity_weibul(n=n, scale=scale_t1,shape=shape_t1,mean_delay = mean_delay_t1 , sd_size=sd_size_t1,times=timeaxis)
sensitivity_family_background <- family_sensitivity_weibul(n=n+150,scale=scale_t1,shape=shape_t1,mean_delay=mean_delay_t1,sd_size=sd_size_t1,times=timeaxis)

sensitivity_average <- generate_mean_of_family(sensitivity_family_background) 


# Plot

# pdf(file='figure_1.pdf',width=7.3,height=4.7)

plot(timeaxis,sensitivity_family_1[,1],type='l',xaxt='n',xaxs='i',yaxs='i',xlim=c(mean_delay_t1-4*sd_size_t1,mean_delay_t1+shift_to_half_likelihood_weibul(shape=shape_t1,scale=scale_t1)+2.35*sd_size_t1),ylim=c(-.005,1.005),xlab="",ylab="",col='green',yaxt='n',bty='L' )

# plot(timeaxis,sensitivity_average,col=col_truth,lwd=5,type='l',xaxt='n',xaxs='i',yaxs='i',xlim=c(mean_delay_t1-4*sd_size_t1,mean_delay_t1+shift_to_half_likelihood_weibul(shape=shape_t1,scale=scale_t1)+2.35*sd_size_t1),ylim=c(-.005,1.005),xlab="",ylab="",yaxt='n',bty='L')

# todo:   title
# label sizes
# colors
# line widths
# scale
#  you're right (alex) about the size of the lines rendering correctly in the pdf
title(xlab="Time since infection", line=1.5, cex.lab=1.2)
title(ylab='Probability of Infection', line=2, cex.lab=1.2)

## goto_fix
# axis ticks, remove box, bring lower axis up to zero or thereabout
yaxis_pos <- c(0,1)
yaxis_names <- c('0','1')


zero_pos <- c(0)
zero_name <- c(expression('0'['']))

axis(side=2, at=yaxis_pos, labels= yaxis_names,tck=-0.037, padj=.237)
#axis(side=1, at=xaxis_pos, labels= xaxis_names,padj=-.35,hadj=-.137)
axis(side=1, at=zero_pos, labels=zero_name,padj=-0.45,hadj=0.37)




for (i in seq(1:n)){
  lines(timeaxis,sensitivity_family_1[,i],col=col_negative,lwd=1.5)
}
lines(timeaxis,sensitivity_average,col=col_truth,lwd=5)


# Average delay arrows

Arrows(x0 = mean_delay_t1-4*sd_size_t1, y0 = 0.501 ,x1= timeaxis[which.min(abs(sensitivity_average-0.5))], y1 = 0.5,code=3 ,arr.adj=1)

text(x=(mean_delay_t1-4*sd_size_t1+timeaxis[which.min(abs(sensitivity_average-0.5))])/2, y=c(0.53), pos=4, labels=expression(italic("d")))

#Standard deviation arrows

Arrows(x1=timeaxis[which.min(abs(sensitivity_average-0.5))],y0=0.5, x0= timeaxis[which.min(abs(sensitivity_average-0.5))]+sd_size_t1,y1=0.5,code=3,arr.adj=1)
text(x=  timeaxis[which.min(abs(sensitivity_average-0.5))]+sd_size_t1/2-.5,y=0.53,pos=4,labels=expression(sigma))
# Arrows(c(0,1.7),c(1.3,-1.8),c(0.8,1.1),c(1.2,-1), lwd=2

# dev.off()


#######

#Figure 2 adapted to represent blood testing scenario






# Figure 2: Likelihood of observed discordant test results, t1 negative t2 positive - different times
#goto_2

# Parameters

n=10
detail=10
timeaxis=seq(0,100,1/detail)

#                       TEST 1 (negative)

#     Describe individual (person) test sensitivity form with population mean-delay and standard deviation of delay
# delay is the variable we distribute across the population - it could in principle be anything else of course
#so each individual has the same SHAPE of sensitivity, but different delays

## Visuals

lwd_means <- 4
lwd_ind <- 1.37
lwd_likelihood <- lwd_means - 3




#High scale causes slower swap
#high shape causes quicker and steeper swap AND more symetrical swap

scale_t1= 2
shape_t1 = 1.73

mean_delay_t1 = 12
sd_size_t1 = 10
#   Time of negative test (relative to arbitrary t=0)

test_time_1 = 48

##                      TEST 2 (positive)

scale_t2 = scale_t1*4.7
shape_t2 = shape_t1

mean_delay_t2 = mean_delay_t1*2+10
sd_size_t2 = 1.3*sd_size_t1

#   Time of positive test
test_time_2 = timeaxis[length(timeaxis)]-5



## Generate Data
# Data = the individual likelihood curves for the first (negative) and second (positive) test
#Test 1
plotdata_negative = family_negative_likelihood_weibul(n=n, scale=scale_t1, shape=shape_t1, mean_delay=mean_delay_t1, sd_size= sd_size_t1, times = timeaxis, test_time = test_time_1)
#for generating mean curve
plotdata_negative_background <- family_negative_likelihood_weibul(n=n+150, scale=scale_t1, shape=shape_t1, mean_delay=mean_delay_t1, sd_size= sd_size_t1, times = timeaxis, test_time = test_time_1)
#Test 2
plotdata_positive = family_positive_likelihood_weibul(n=n, scale=scale_t2, shape=shape_t2, mean_delay=mean_delay_t2, sd_size= sd_size_t2, times = timeaxis, test_time = test_time_2)
#for generating mean curve
plotdata_positive_background <- family_positive_likelihood_weibul(n=n+150, scale=scale_t2, shape=shape_t2, mean_delay=mean_delay_t2, sd_size= sd_size_t2, times = timeaxis, test_time = test_time_2)


##        COMMUNICATE

# pdf(file = "simple_case.pdf", width = 7.3, height = 4.7)

plot(timeaxis,plotdata_negative[,1],type='l',xlim=c(timeaxis[1],timeaxis[length(timeaxis)]),ylim=c(-.002,1.002),xaxt='n',yaxt='n',xaxs='i',yaxs='i',bty='l',xlab='',ylab='',col='green') #clarify label in comment
title(xlab="Hypothetical time of infection", line=1.5, cex.lab=1.4)
title(ylab=expression('Likelihood'), line=2, cex.lab=1.4)


yaxis_pos <- c(0,0.5,1)
yaxis_names <- c('0',"",'1')

xaxis_pos <- c(test_time_2)
xaxis_names <- c(expression('T'['Donation']))

zero_pos <- c(0)
zero_name <- c(expression('0'['']))
#shift axes
axis(side=2, at=yaxis_pos, labels= yaxis_names,tck=-0.037, padj=.437)
axis(side=1, at=xaxis_pos, labels= xaxis_names,padj=-.35,hadj=0)#-.137)
#axis(side=1, at=zero_pos, labels=zero_name,padj=-0.45,hadj=0.37)

# goto_do

#points(plotdata[,1],plotdata[,3])
for (i in seq(1:n)){
  lines(timeaxis,plotdata_negative[,i],col=col_negative,lwd=lwd_ind)
}
#points(plotdata[,1],plotdata[,3])
for (i in seq(1:n)){
  lines(timeaxis,plotdata_positive[,i],col=col_positive,lwd=lwd_ind)
}

positive_mean_background <- rowMeans(plotdata_positive_background)
negative_mean_background <- rowMeans(plotdata_negative_background)
lines(timeaxis,negative_mean_background, lwd=lwd_means, col=col_negative)
lines(timeaxis,positive_mean_background, lwd=lwd_means, col=col_positive)


naive_likelihood <- positive_mean_background*negative_mean_background
lines(timeaxis,naive_likelihood,col='grey42',lwd=4,lty=5)

true_likelihood <- likelihood_by_DDI(plotdata_positive_background,plotdata_negative_background,timeaxis)
lines(timeaxis,true_likelihood,col=col_truth,lwd=4)


# segments(x0=test_time_1,y0=0,x1=test_time_1,y1=1,lty=4)
segments(x0=test_time_2,y0=0,x1=test_time_2,y1=1,lty=4)

# mean delay arrows (shape package)
# 
# Arrows(x0 = timeaxis[which.min(abs(negative_mean_background-0.5))],y0=0.5,x1=test_time_1,y1=0.5,code=3 ,arr.adj=1,arr.width = .1/2)
# text(x=(timeaxis[which.min(abs(negative_mean_background-0.5))]+test_time_1)/2-2.4,y=0.53,pos=4,labels=expression(italic('d')['1']))
# 
# Arrows(x0 = timeaxis[which.min(abs(positive_mean_background-0.5))],y0=0.5,x1=test_time_2,y1=0.5,code=3 ,arr.adj=1,arr.width = .1/2)
# text(x=(timeaxis[which.min(abs(positive_mean_background-0.5))]+test_time_2)/2-2.5,y=0.53,pos=4,labels=expression(italic('d')['2']))

# Standard deviation arrows

Arrows(x0 = timeaxis[which.min(abs(negative_mean_background-0.5))]-sd_size_t1,y0=0.5,x1= timeaxis[which.min(abs(negative_mean_background-0.5))], y1=0.5, code=3,arr.adj=1,arr.length = .1,arr.width = .1/2)
text(x=timeaxis[which.min(abs(negative_mean_background-0.5))]-sd_size_t1/2-2.4,y=0.53,pos=4,labels=expression(sigma['1']))

Arrows(x0 = timeaxis[which.min(abs(positive_mean_background-0.5))]-sd_size_t2,y0=0.5,x1= timeaxis[which.min(abs(positive_mean_background-0.5))], y1=0.5, code=3,arr.adj=1,arr.length = .1,arr.width = .1/2)
text(x=timeaxis[which.min(abs(positive_mean_background-0.5))]-sd_size_t2/2-2.5,y=0.53,pos=4,labels=expression(sigma['2']))

# Delta arrow
# Arrows(x0 = test_time_1,y0=0.75,x1=test_time_2,y1=0.75, code=3 ,arr.adj=1,arr.width = .1/2)
# text(x = test_time_1 + (test_time_2 - test_time_1)/2 - 5, y=0.8,pos=4,labels=expression(delta))


#product of means
#true likelihood





# dev.off()
