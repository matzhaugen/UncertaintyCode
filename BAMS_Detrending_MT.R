library(ncdf4)
# Create Z.observed
# NCEP Reanalysis Data
filename1 = "hgt500_npac_annual_1948.2014.nc"
f1 			   = nc_open(filename1)
hgt_annual		= as.numeric(ncvar_get(f1, "timeseries_annual"))
Z.observed 		= hgt_annual[32:65] # formerly called hgt.annual in code

# Get 2013 observation
dat.ncep <- read.table("hgt.ncep")
dat.ncep <- as.numeric(as.matrix(dat.ncep))
max.dat.ncep = max(dat.ncep)

load("z.giss.hist.new") # New as of June 4, 2014
load("z.hadgem.hist.new") # New as of June 4, 2014
load("z.noresm.hist.new") # New as of June 4, 2014

hgt.giss <- na.omit(hgt.giss.hist)
hgt.noresm <- na.omit(hgt.noresm.hist)
# NA in 27th and final place of hadgem data.
hgt.hadgem <- c(na.omit(hgt.hadgem.hist),NA) # put in final NA
hgt.hadgem <- c(hgt.hadgem[1:26],NA,hgt.hadgem[27:length(hgt.hadgem)]) # put in NA at 27

detrend <- function(vec,n=27){
	N <- length(vec)/n
	vec.list <- list()
	for(i in 1:N){
		#print((i-1)*n + 1:n) # sanity check
		vec.list[[i]] <- vec[(i-1)*n + 1:n]
	}
	fit.list <- lapply(vec.list,function(x){lm(x ~ I(1:n))})
	resids.mat <- sapply(fit.list,function(x){x$resid})
	list("vec"=vec.list,"fit"=fit.list,"resid"=resids.mat)
}

dt.giss <- detrend(vec=hgt.giss)
dt.hadgem <- detrend(vec=hgt.hadgem)
dt.noresm <- detrend(vec=hgt.noresm)

plot.ts(c(dt.giss$resid,unlist(dt.hadgem$resid),dt.noresm$resid))

# Slope coefficients for each realization
sapply(dt.giss$fit,function(x){x$coef[2]})
sapply(dt.hadgem$fit,function(x){x$coef[2]})
sapply(dt.noresm$fit,function(x){x$coef[2]})

# Leave out 2013 observation from fitting the trend to the observed series.
# But: We have to remove the predicted 2013 value using the fitted trend in order to
# compare the 2013 observation to the detrended series.

dt.observed <- detrend(c(Z.observed),n=length(c(Z.observed)))

# Predicted 2013 value using the observed pre-2013 linear trend.
pred.2013 <- summary(dt.observed$fit[[1]])$coef[1,1] + summary(dt.observed$fit[[1]])$coef[2,1]*35
pred.2013
plot.ts(Z.observed)
plot.ts(c(Z.observed,max.dat.ncep))
abline(dt.observed$fit[[1]])
points(35,pred.2013)

c(dt.observed$resid,max.dat.ncep - pred.2013)

Z.giss.dt <- c(dt.giss$resid)
Z.hadgem.dt <- c(dt.hadgem$resid)
Z.noresm.dt <- c(dt.noresm$resid)
Z.observed.dt <- c(dt.observed$resid)
obs.2013.dt <- max.dat.ncep - pred.2013

#save(Z.giss.dt,Z.hadgem.dt,Z.noresm.dt,Z.observed.dt,obs.2013.dt,file="Detrended.RData")