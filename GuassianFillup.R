# This software is Copyright © 2017 The Regents of the University of California. All Rights Reserved.
#  
# Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies.
#  
# Permission to make commercial use of this software may be obtained by contacting:
# Office of Innovation and Commercialization
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# invent@ucsd.edu

# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#  
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS.


#!/usr/bin/R
# Contact: Shicheng Guo
# Version 1.3
# Update: Dec/29/2016
# Function: Estimated missing single-tailed data from Guassian distribution (one-side)
# Input: input variable, assume to be Guassian distribution, other distribution should be replace with kenal

# usage: 
GuassianFillup(input)

# please load the 
GuassianFillup<-function(input){
input<- as.vector(input)
data<-density(input)
peak<-data$x[which.max(data$y)]
bound5<-data$x[which.max(data$y)]-sd(input)
bound3<-data$x[which.max(data$y)]+sd(input)
bound5
bound3
NumA<-sum(input<=peak)
NumB<-sum(input>=peak)
Estimatepeak<-seq(bound5,bound3,length.out=1000)[which.min(err)]
if(NumA>NumB){
  print(paste("Esitmated Peak=",round(peak,3),"please fill up N=",abs(NumA-NumB),"data point in the rigth side",sep=" "))
}else{
  print(paste("Esitmated Peak=",round(peak,3),"please fill up N=",abs(NumA-NumB),"data point in the left side",sep=" "))
}

inputsim<-c(input[input>=Estimatepeak],unlist(lapply(input[input>Estimatepeak],function(x) 2*Estimatepeak-x)))
sdsim<-sd(inputsim)
Max<-max(max(input)-peak,peak-min(input))
hist(input,breaks=100,xlim=c(peak-Max,peak+Max),prob=TRUE)
abline(h=0,lty=2,col="black",lwd=2)
abline(v=peak,lty=2,col="red",lwd=2)
lines(density(input, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
lines(density(rnorm(length(input),peak,sdsim), adjust=2), lty="dotted", col="darkblue", lwd=2) 
legend("topright",legend = c("Estimated peak","Raw fitted density","Fill up density"),lty="dotted",lwd=2,col=c("red","darkgreen","darkblue"),bty = "n")
}

# example: 
# Simulated input data
# set.seed(146)
# input<-rnorm(10000,0,1)
# hist(input,breaks=100,xlim=c(peak-Max,peak+Max),prob=TRUE)
# input<-input[input>-1.5]
# GuassianFillup(input)
