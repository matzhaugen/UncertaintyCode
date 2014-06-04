#dir = "/Users/matzhaugen/GoogleDrive/Research/DSwainCollaboration/SharedBAMS/UncertaintyCode"
#setwd(dir)

dir = "/Users/matzhaugen/GoogleDrive/Research/DSwainCollaboration/SharedBAMS/UncertaintyCodeDaniel"
setwd(dir)

library(ncdf4)
library(nortest)
library(np)
library(fields)
library(forecast)
library(exactRankTests)
library(ggplot2)
library(plyr)
library(perm)
library(DAAG)
library(gregmisc)
library(evd)
library(vcd)
library(fitdistrplus)
library(VGAM)
library(MASS)

#### MATZ's code ###############

# Import data #
# 3 vectors are:
# 1. observational: total_baseline
# 2. preindustrial: get_pi_pr_mod
# 3. historical: get_hist_pr_mod
#preindustrial = get_pi_pr_mod
#historical = get_hist_pr_mod
#Z.observed = total_baseline
#Select a model
#model = 1
#Z.historical = c(get_hist_pr_mod[,,model])
#Z.preindustrial = c(get_pi_pr_mod[,,model])
#Z.historical = Z.historical[!is.na(Z.historical)]
#Z.preindustrial = Z.preindustrial[!is.na(Z.preindustrial)]

load("z.giss.hist.new")
load("z.hadgem.hist.new")
load("z.noresm.hist.new")
load("z.giss.pi.new")
load("z.hadgem.pi.new")
load("z.noresm.pi.new")
load("z.obs.new")

Z.historical = c(hgt.giss.hist,hgt.hadgem.hist,hgt.noresm.hist)
Z.preindustrial = c(hgt.giss.pi,hgt.hadgem.pi,hgt.noresm.pi)
Z.historical = Z.historical[!is.na(Z.historical)]
Z.preindustrial = Z.preindustrial[!is.na(Z.preindustrial)]
Z.observed = hgt.annual.full[1:34]

total_2013 = max(hgt.annual.full)



# Run Bootstrap #
set.seed(3)
B = 10
nmodels = 1
ratios = matrix(0, B, nmodels)
for (i in 1:nmodels) {
	model = i
	#Z.historical = c(get_hist_pr_mod[,,model])
	#Z.preindustrial = c(get_pi_pr_mod[,,model])
	Z.historical = Z.historical[!is.na(Z.historical)]
	Z.preindustrial = Z.preindustrial[!is.na(Z.preindustrial)]
	
	ratios[,i] = output(Z.observed, Z.historical, Z.preindustrial, total_2013, B, plotit=T)
}
#################


#Plotting - ecdf
plot(ecdf(ratios[,1]), main="")
for (i in 2:nmodels) {
	lines(ecdf(ratios[,i]),col=i)
}


hist(ratioPreOverObs[ratioPreOverObs<3], main="")
abline(v=ratioPreOverObs[1])
abline(v=median(ratioPreOverObs),col=2)
legend("topright", c("Point Estimate", "Median"),col=c(1,2), lty=1)
confidenceIntervalRatio = c(quantile(ratioPreOverObs,0.025), quantile(ratioPreOverObs,0.975))#95th confidence interval
title(paste("conf int:",signif(confidenceIntervalRatio[1],3),signif(confidenceIntervalRatio[2],3)))
plot(ecdf(ratioPreOverObs),xlim=c(0,2))
##############
par(mfrow=c(3,3)); for (i in 1:9) {hist(paramsStar[,i])}



#### Functions ####
loglikParetoiii = function(params, obs) {
	-sum(log(dparetoIII(obs,location=params[1],scale=params[2],inequality=params[3])))
}

plotZs = function(obs, params, ...) {
	hist(obs,breaks=10,col="lightgrey",prob=TRUE, ...)
	xmin=range(obs)[1];	xmax=range(obs)[2]
	u = seq(xmin,xmax,length=100)
	lines(u,dparetoIII(u, location =params[1], scale = 	params[2],inequality=params[3]),col="green",lwd=2)
}
fitPareto = function(z) {
	n=10	
	p2Window = 0.5
	Ez = mean(z)
	SDz = sd(z)
	p1 = min(z)-1/2*SDz
	p3 = 0.1
	p2 = Ez/gamma(1-p3)/gamma(1+p3)
	lower = c(p1-SDz,p2*(p2Window), 0.001)
	upper = c(min(z)-SDz/1000, p2*(2-p2Window),0.4)	
	
	optPar = optim(c(p1,p2,p3),loglikParetoiii,obs=Z.observed,method="L-BFGS-B",
		lower = lower, upper=upper)
	
	#update parameters
	p1 = optPar$par[1];
	p2 = optPar$par[1];	
	p3 = optPar$par[1];
	lower = c(p1-SDz,p2*(p2Window), 0.001)
	upper = c(min(z)-SDz/1000, p2*(2-p2Window),0.4)	
	
	locv = seq(lower[1],upper[1],len=n)
	locu = seq(lower[2],upper[2],len=n)
	mins = matrix(0, n,n)
	for (i in 1:n) {
		for (j in 1:n) {
			mins[i,j] = optim(c(locv[i],locu[j],fit$est[3]),loglikParetoiii,obs=z,method="L-BFGS-B", lower=lower,upper=upper)$value
		}
	}
	optind = which(mins == min(mins), arr.ind=T)
	out = optim(c(locv[optind[1]],locu[optind[2]],fit$est[3]),
		loglikParetoiii,obs=z,method="L-BFGS-B", 		
		lower=lower,upper=upper, hessian=T)
}


########### ########### ########### ########### 
#### get Information matrix and return period with confidence intervals #####
########### ########### ########### ########### ########### 
B = 1000
hessian = array(0,dim = c(B,3,3))
paramStar = matrix(0,B,3)
lower = c(-100,50,0.01)
upper = c(15,200,0.5)
nobs = length(Z.observed)
par1 = runif(B,lower[1],upper[1]);par2 = runif(B,lower[2],upper[2]);par3 = runif(B,lower[3],upper[3])
returnPeriod = rep(0,B)
max.dat.ncep = total_2013
par(ask=F)

for (i in 1:B) {
	obsStar = Z.observed[resample(1:nobs,nobs,replace=T)]
	out = optim(c(par1[i],par2[i],par3[i]),loglikParetoiii, obs=obsStar,
	 method="L-BFGS-B", lower=c(5000,50,0.01),upper=c(5500,250,0.30), hessian=T)
	returnPeriod[i] = signif(1/(1 - pparetoIII(max.dat.ncep,out$par[1],out$par[2],out$par[3])), 3)
	hessian[i,,] = out$hess #Hessian
	paramStar[i,] = signif(out$par,4)
	print(paste(i,"Return Period",returnPeriod[i]," Parameters:", paramStar[i,1], paramStar[i,2], paramStar[i,3], signif(out$val,3)))
	#plotZs(obsStar, paramStar[i,], xlim=range(Z.observed),
	#main=paste("Return Period=",returnPeriod[i]))
}
# max(dat.ncep) is the observed data point from 2013 that we want to compute a return period for.
#Return period = 1/(Tail probability of an event) >
# input whatever loc, scal, ineq params from a fit
EReturnPeriod  = median(returnPeriod) #median return period
sdReturnPeriod = sd(returnPeriod)
Ehessian       = apply(hessian,c(2,3),mean)
Epar           = colMeans(paramStar)
covPar         = cov(paramStar) #We note a large magnitude in the off-diagonal of entry 12
confidenceInterval = c(quantile(returnPeriod,0.975), quantile(returnPeriod,0.025)) #95th confidence interval
########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### 

#This function estimates parameter bounds of the pareto iii distribution. 
#Use squared error loss
paretoFirstMoment = function(inequality, scale, firstMoment) {
	output = (firstMoment - scale*gamma(1-inequality)*gamma(1+inequality))^2
	output
}

paretoFirstSecondMoment = function(u, firstMoment, secondMoment, inequality =0.01) {
	output1 = (firstMoment - u[1] - u[2]*gamma(1-inequality)*gamma(1+inequality))^2
	output2 = (secondMoment - u[1]^2 - u[2]^2*gamma(1-2*inequality)*gamma(1+2*inequality) 
					- 2*u[1]*u[2]*gamma(1-inequality)*gamma(1+inequality))^2
	output1 + output2
}
getBounds = function(z) {
	p2Window = 0.5
	p3Window = 0.5
	p2factor = 5
	p3factor = 5	
	Ez = mean(z)
	SDz = sd(z)
	p1 = min(z)-SDz
	p3 = 0.1
	p2 = mean(z-p1)/gamma(1-p3)/gamma(1+p3)
	#Update the first two parameters
	out1 = optim(c(p1,p2),paretoFirstSecondMoment, firstMoment=Ez, secondMoment=mean(z^2))$par
	p1 = out1[1]; p2 = out1[2]
	lower = c(p1-100*SDz,p2/p2factor, 0.001)
	upper = c(min(z)-SDz/1000, p2*p2factor,min(p3factor*p3,0.49))	
	
	optPar = optim(c(p1,p2,p3),loglikParetoiii,obs=z,method="L-BFGS-B",
		lower = lower, upper=upper)
	
	#update parameters
	p1 = optPar$par[1];
	p2 = optPar$par[2];	
	p3 = optPar$par[3];
	print(paste("Optimal bounds: ", p1,p2,p3))
	lower = c(p1-SDz,p2*(p2Window), max(p3Window*p3,0.01))
	upper = c(min(z)-SDz/1000, p2*(2-p2Window),min((2-p3Window)*p3,0.99))	
	list(lower=lower,upper=upper, optimal=c(p1,p2,p3))
}
########### ########### ########### ########### ########### ########### ########### ########### 
######### Get return ratio between historical and preinducstrial using 
######### the observed as a benchmark 
output = function(Z.observed, Z.historical, Z.preindustrial, maxPoint, B=100, plotit=F) {	
	###### 
	boundsObs = getBounds(Z.observed)
	boundsHist = getBounds(Z.historical)
	boundsPre = getBounds(Z.preindustrial)
	
	print("Bounds for the 3 data sets")
	print(boundsObs)
	print(boundsHist)
	print(boundsPre)
	optPar = optim(c(p1,p2,p3),loglikParetoiii,obs=z,method="L-BFGS-B",
		lower = lower, upper=upper)

	nobs = length(Z.observed)
	nhist = length(Z.historical)
	npre = length(Z.preindustrial)
	parObs = cbind(runif(B,boundsObs$lower[1],boundsObs$upper[1]),
					runif(B,boundsObs$lower[2],boundsObs$upper[2]), 
					runif(B,boundsObs$lower[3],boundsObs$upper[3]))
	 	
	parHist = cbind(runif(B,boundsHist$lower[1],boundsHist$upper[1]),
					runif(B,boundsHist$lower[2],boundsHist$upper[2]),
					runif(B,boundsHist$lower[3],boundsHist$upper[3]))
	
	parPre = cbind(runif(B,boundsPre$lower[1],boundsPre$upper[1]),
				runif(B,boundsPre$lower[2],boundsPre$upper[2]),
		   		runif(B,boundsPre$lower[3],boundsPre$upper[3]))

	paramsStar 		= matrix(0,B,9)	
	returnPeriodObs = rep(0,B)
	histAtPobs 		= rep(0,B)
	returnPeriodPre = rep(0,B)
	ratioPreOverObs = rep(0,B)
	PpreAtHistAtPobs= rep(0,B)
	max.dat.ncep = maxPoint
	for (i in 1:B) {
		print(i)
		if (i==1) { #Point estimate
			obsStar = Z.observed
			histStar = Z.historical
			preStar = Z.preindustrial
		} else { # Bootstrap
			obsStar = Z.observed[resample(1:nobs,nobs,replace=T)]
			histStar = Z.historical[resample(1:nhist,nhist,replace=T)]
			preStar = Z.preindustrial[resample(1:npre,npre,replace=T)]
		}
		print(parObs[i,])
		print(parHist[i,])

		out = optim(parObs[i,],loglikParetoiii, obs=obsStar,
		 method="L-BFGS-B", lower=boundsObs$lower,upper=boundsObs$upper, hessian=F)
		 
		returnPeriodObs[i] = signif(1/(1 - pparetoIII(max.dat.ncep,out$par[1],out$par[2],out$par[3])), 3)
		pObs = pparetoIII(max.dat.ncep,out$par[1],out$par[2],out$par[3])
		
		outHist = optim(parHist[i,],loglikParetoiii, obs=histStar,
		 method="L-BFGS-B", lower=boundsHist$lower,upper=boundsHist$upper, hessian=F)
		# print(outHist$par)
		histAtPobs[i] = qparetoIII(pObs, outHist$par[1],outHist$par[2],outHist$par[3])
		# fit preindustrial
		outPre = optim(parPre[i,],loglikParetoiii, obs=preStar,
		 method="L-BFGS-B", lower=boundsPre$lower,upper=boundsPre$upper, hessian=F)
		 #print(outPre$par)
		#	plotZs(preStar, outPre$par, xlim=range(Z.preindustrial))
		# Get return period for preindustrial
		PpreAtHistAtPobs[i] = pparetoIII(histAtPobs[i],outPre$par[1],outPre$par[2],outPre$par[3])
		returnPeriodPre[i] = signif(1/(1 - PpreAtHistAtPobs[i]), 3)
		ratioPreOverObs[i] = returnPeriodPre[i]/returnPeriodObs[i]
		paramsStar[i,] = c(out$par, outHist$par, outPre$par)
		par(mfrow=c(1,3))
		if (plotit) {
			par(ask=T)
			plotZs(obsStar, out$par, xlim=range(Z.observed))
			plotZs(histStar, outHist$par, xlim=range(Z.historical))
			plotZs(preStar, outPre$par, xlim=range(Z.preindustrial))
			print(outPre$par)
			print(outPre$value)
		} else {
			par(ask=F)
		}
	}
	######
	ratioPreOverObs
}

#### END MATZ's code ###########
