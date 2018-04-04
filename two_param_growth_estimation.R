## load state data from fiadb containing for each county, 
## age and volume for all species groups
## estimate von Bertalanffy growth parameters

rm(list=ls())

library(matrixcalc)
library(ggplot2)
library(plyr)
library(nlWaldTest)


# function to calculate sum of squared residuals (objective function in optimization)

SSQ <- function(beta, x) {
  a <- beta[1]
  b <- beta[2]
  vhat <- rep(0, length(age))
  epsilon <- rep(0, length(age))
  for (i in 1:length(age)) {
    vhat[i] <- a*(1-exp(-b*age[i]))^3
    epsilon[i] <- (vol[i]-vhat[i])^2
  }
  ssq <- sum(epsilon)
  return(ssq)
}

#load data and id tables

dat <- read.csv("Z:/forest_project/bertalanffy/growth_data.csv")

#create uniue id

dat <- subset(dat, VOLCFNET>=0) #four observations dropped
dat$id <- paste(paste(paste(formatC(dat$STATECD,format='d',width=2,flag='0'),
                      formatC(dat$COUNTYCD,format='d',width=3,flag='0'),sep=''),
                      formatC(dat$SPGRPCD,format='d',width=2,flag='0'),sep=''),
                      formatC(dat$FORTYPCD,format='d',width=3,flag='0'),sep='')

#initialize object to store solutions

vbsol <- data.frame(matrix(0,length(unique(dat$id)),6))
vbsol[,1] <-  unique(dat$id)

# loop for each unique id (county,species) to estimate parameters

ptm <- proc.time()

for (i in 1:length(vbsol$X1)) { 
  
  u <- vbsol[i,1]
  fips <- substring(u,1,5)
  species <- substring(u,6,10)
  
  age <- dat$STDAGE[dat$id==u]
  vol <- dat$VOLCFNET[dat$id==u]
  
  # plot(age,vol)
  
  #initial conditions
  
  beta <- c(0,0)
  beta[1] <- 0.5*max(vol)
  beta[2] <- .1*runif(1)
  
  fit <- optim(beta, fn = SSQ, method = "BFGS",hessian = TRUE)
  
  if (is.positive.definite(fit$hessian) & class(try(solve(fit$hessian), silent = T))!='try-error') {
    fit$V <- solve(fit$hessian)  #solve for the inverse hessian
    if (diag(fit$V)[1] >= .00001) {
      vbsol[i,2] <- fit$par[1]
      vbsol[i,3] <- fit$par[2]
      vbsol[i,4] <- fit$V[1,1] #variance of a par
      vbsol[i,5] <- fit$V[2,2] #variance of b par
      vbsol[i,6] <- fit$V[1,2] #covariance of a,b
      
    }
  }
}

  
  #graph results
  
#   age_seq <- (seq(0,max(vol),1))
#   vol_pred <- fit$par[1]*(1-exp(-fit$par[2]*age_seq))^3
#   plot(age,vol)
#   title(unique(paste(georef$county_name[georef$county_fips==fips],
#                      georef$state_abbr[georef$county_fips==fips],
#                      spcodes$common_name[spcodes$spgrpcd==spgrpcd])))
#   lines(age_seq,vol_pred,lty=1,col='red',lwd=3)
  

proc.time() - ptm

saveRDS(vbsol, "Z:/forest_project/bertalanffy/vbsol_raw.rds")

########################################################################################
## Clean Output
########################################################################################

# wald test


confint <- nlConfint(coeff=fit$par, Vcov=fit$V, texts = c('b[1]','b[2]'))


remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.01,.99), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

output <- vbsol
output <- subset(output, X6 > 1.645 & X7 > 1.645)
output <- rename(output, c('X1'='id', 'X2'='a_par', 'X3'='b_par',
                           'X4'='a_se', 'X5'='b_se','X6'='a_t-stat','X7'='b_t-stat',
                           'X8'='ab_corr'))
output$a_par <- remove_outliers(output$a_par)
output <- subset(output, !is.na(a_par))
output$b_par <- remove_outliers(output$b_par)
output <- subset(output, !is.na(b_par))

ggplot(output)+
  geom_histogram(aes(x=a_par),fill = "red", alpha = 0.2)
  
ggplot(output)+
  geom_histogram(aes(x=b_par),fill = "blue", alpha = 0.2)

a <- 79.42315
b <- .02576
x <- seq(0,100,1)
y <- a*(1-exp(-b*x))^3

plot(x,y,type='l')

write.csv(output, 'growth_parameters.csv')





# ################################################################################### Benton & Deschutes OR
# #################################################################################
library(ggplot2)
library(RSQLite)
library(reshape2)
library(directlabels)
setwd('Z:/')
conn <- dbConnect(SQLite(), 'database/main.db')
dbListTables(conn)


d <- read.csv('Z:/growth_data.csv')
g <- dbReadTable(conn, "growth_parameters_old")
sub <- subset(g, id==4100310 | id==4101710)

x=seq(0,75,1)
y=seq(0,150,1)
Benton=sub$a_par[sub$id==4100310]*(1-exp(-sub$b_par[sub$id==4100310]*x))^3
Benton2=sub$a_par[sub$id==4100310]*(1-exp(-sub$b_par[sub$id==4100310]*y))^3
Deschutes=sub$a_par[sub$id==4101710]*(1-exp(-sub$b_par[sub$id==4101710]*x))^3
plot <- data.frame(x,Benton,Deschutes)
plot2 <- data.frame(x,Benton)
plot3 <- data.frame(y,Benton2)
plot <- melt(plot, id.vars = 'x')

ggplot(data=plot2,aes(x=x,y=Benton))+
	geom_line()+
	xlab('Stand Age')+
	ylab('Growing Stock Volume')+
	ggtitle('Douglas Fir Growth Function')+
  geom_vline(aes(xintercept = 75),color='green', linetype="dashed", size=1.2)+
	theme_minimal()

ggplot()+
  geom_line(data=plot3,aes(x=y,y=Benton2))+
  xlab('Stand Age')+
  ylab('Growing Stock Volume')+
  ggtitle('Douglas Fir Growth Function')+
  geom_vline(aes(xintercept = 75),color='green', linetype="dashed", size=1.2)+
  theme_minimal()

direct.label(ggplot(data=plot,aes(x=x,y=value,group=variable,color=variable))+
	geom_line(color='black')+
	xlab('Stand Age')+
	ylab('Growing Stock Volume')+
	ggtitle('Douglas Fir Growth Function')+
	theme_minimal(), list(last.points, hjust = 0.8,
												vjust = 2.8))

bdf_obs <- subset(d, COUNTYCD==3 & STATECD==41 & SPGRPCD==10)

ggplot(data=bdf_obs, aes(x=STDAGE,y=VOLCFNET))+
	geom_point()+
	xlab('Stand Age')+
	ylab('Growing Stock Volume')+
	ggtitle('Douglas Fir Tree Observations')+
	annotate('text', x=100, y=1700, label='n = 1,426', size=5)+
	theme_minimal()




