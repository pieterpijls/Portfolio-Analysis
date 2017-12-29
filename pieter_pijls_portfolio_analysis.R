## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.pos = 'H')

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----  include=FALSE-----------------------------------------------------
rm(list=ls())
set.seed(0387948)

## ----  include=FALSE-----------------------------------------------------
#Load packages
packages <- c("copula","mgcv", "fCopulae", "Ecdat", "fGarch", "MASS", "tidyverse", "ggthemes","plotly","ggplot2","xtable","knitr","statmod", "actuar", "fitdistrplus", "PerformanceAnalytics", "gridExtra", "grid", "scales", "quadprog","quantmod","plyr","reshape","kableExtra","xts","float","lubridate","PortfolioAnalytics")

# This function will load package or if not already installed install it
suppressMessages(packages <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}))

## ----  include=FALSE-----------------------------------------------------
# Load data
load("~/Documents/KU LEUVEN/MFAE/STATISTICS FOR FINANCE&INSURANCE/Assignment 2 Portfolio Analysis/data_returns_387948.RData")
colnames(data)<-c("HON","SUN","STW","DRH","AVI","CBL","MER","GPS","MAS","BLC")
data_backup <- data

## ----fig1, echo=FALSE, fig.cap="\\label{fig:fig1}QQ-plots Daily Stock Returns"----
#QQplots
par(mfrow=c(2,5))
for (i in 1:ncol(data))
{
  qqnorm(data[,i],main=colnames(data[,i]))
  qqline(data[,i])
}

## ----  include=FALSE-----------------------------------------------------
#data for plot
gdata <- as.data.frame(data)
gdata <- cbind(row.names(gdata),gdata)
gdata  <- melt(gdata,id=c("row.names(gdata)"))
gdata[,1]<-as.Date(gdata[,1])
colnames(gdata) <- c("Date","Equity","Return")

## ----fig2, echo=FALSE, fig.cap="\\label{fig:fig2}Daily Stock Returns from 1999-2016"----
#daily returns stocks
p<-ggplot(gdata,aes(x=Date,y=Return)) + geom_point(aes(color=Equity),size=0.5) + scale_x_date(expand = c(0, 0),name="Year",labels = date_format("%Y"), breaks=date_breaks("2 years"),limits=c(as.Date(c("1999-01-04","2016-12-30")))) 
p
#ggplotly(p)

## ---- echo=FALSE---------------------------------------------------------
kable(round(cor(data),2), format = "latex", caption = "Correlation Matrix", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----  include=FALSE-----------------------------------------------------
#take first 126 observations
completedata<-data
data<-data[1:126,]
Returns.xts <- xts(data*100,order.by = index(data)) # multiply by 100 !!

## ----  include=FALSE-----------------------------------------------------
# By the numbers
Returnsmean <- apply(Returns.xts,2,mean)
Returnssd <- apply(Returns.xts,2,sd)
Returnsmedian <- apply(Returns.xts,2,mean)
Returnsqn <- apply(Returns.xts,2,qn) #trade off efficiency and robustness 

## ----  include=FALSE-----------------------------------------------------
#summary stats
Returns <- as.data.frame(data)
Returns$Year <- year(index(data))
Returns$Month <- month(index(data))
Returns <- cbind(row.names(as.data.frame(data)),Returns)
colnames(Returns)<-c("Date",colnames(data),"Year","Month")

Monthlysd <- ddply(Returns,.(Year,Month),colwise(sd))
Monthlysd$Date <- strptime(sprintf("%s-%02d-01",Monthlysd$Year,Monthlysd$Month),"%Y-%m-%d")
msdmelt <- melt(Monthlysd[,-(1:2)],id="Date")
colnames(msdmelt)[2:3] <- c("Equity","std_Return")
ggplot(msdmelt) + geom_line(aes_string(x="Date",y="std_Return",color="Equity"),size=1.5) + ggtitle("Monthly Rolling Standard Deviations of Daily Returns of 10 shares") + theme(panel.background = element_rect(fill="#CCCCEE"), legend.key = element_rect(fill="#CCCCEE"))

## ----  include=FALSE-----------------------------------------------------
# Estimate covariance matrix
Returnscov <- cov(Returns.xts)

## ---- echo=FALSE---------------------------------------------------------
kable(round(Returnscov,2), format = "latex", caption = "Coviarance Matrix", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----error=TRUE, include=FALSE-------------------------------------------
# No short selling - we add an extra constraint e.g. mutual funds 
Amat <- cbind(rep(1,5), Returnsmean[1:10], diag(1,10),diag(-1,10)) # we reconstruct A matrix

N <- length(Returnsmean) 
Amat <- cbind(rep(1,N), Returnsmean, diag(1,nrow=N), diag(-1,nrow=N)) 

# closure constraint
# diag(1,6) 
nrep <- 1000 # we make 1000 replications
muP <- seq(0.001, 0.06, length=nrep) # nrep targets is a sequence from 0,1% to 6% in 1000 steps
sdP <- rep(0, nrep) # 1000 empty to save standard deviation
weights <- matrix(0, nrow=nrep, ncol=10) # 1000 empty matrices were we will save the weights. 

for (i in 1:nrep)
{
  bvec <- c(1, muP[i], rep(0,N), rep(-1,N)) # constraint vector
  result <- solve.QP(Dmat=2*Returnscov, dvec=rep(0,10), Amat=Amat, bvec=bvec,
                     meq=2)  
  # meq=2: all constraints so far are equality constraints
  sdP[i] <- sqrt(result$value)
  weights[i,] <- result$solution
}
# Crashes, why? 
# Because no short selling constraint limits feasible return within range of indivdual equity returns 
# The return cannot be achieved 

## ----  include=FALSE-----------------------------------------------------
# This will work: because we adjust the returns 
muP <- seq(min(Returnsmean)+0.00001, max(Returnsmean)-0.00001, length=nrep)
sdP <- rep(0, nrep) # 1000 empty to save standard deviation
weights <- matrix(0, nrow=nrep, ncol=10) # 1000 empty matrices were we will save the weights. 

for (i in 1:nrep)
{
  bvec <- c(1, muP[i],rep(0,10), rep(-1,N))  # constraint vector
  result <- solve.QP(Dmat=2*Returnscov[1:10,1:10], dvec=rep(0,10), Amat=Amat, bvec=bvec,
                     meq=2)  
  # meq=2: all constraints so far are equality constraints
  sdP[i] <- sqrt(result$value)
  weights[i,] <- result$solution
}

## ----  include=FALSE-----------------------------------------------------
## Efificient frontier and Minimum Variance Portfolio

plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
# Minumum Variance portfolio
MV <- (sdP == min(sdP)) 
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
# Effective portfolios
EFF <- (muP > muP[MV])
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
weights[MV,]

## ----  include=FALSE-----------------------------------------------------
## Risk free portfolio
mufree <- 0.03 # input value of risk-free interest rate

plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
points(0, mufree, cex=2, pch=15) 

## ----  include=FALSE-----------------------------------------------------
## Tangency Portfolio
sharpe <- (muP-mufree)/sdP 
TAN <- (sharpe == max(sharpe)) 
weights[TAN,] 
sdpT <- sdP[TAN]
mupT <- muP[TAN]

## ----  include=FALSE-----------------------------------------------------
# Plot 
plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
points(0, mufree, cex=2, pch=15) 
points(sdP[TAN], muP[TAN], cex=2, pch=16,col='blue') 
lines(c(0,5), mufree+c(0,5)*(muP[TAN]-mufree)/sdP[TAN], lwd=4, lty=1,
      col="blue") 
for (i in 1:ncol(data))
{
  text(Returnssd[i], Returnsmean[i], colnames(data[,i]), cex=1)
}

## ----  include=FALSE-----------------------------------------------------
## Minimum risk contribution portfolio
library(DEoptim)

obj <- function(w) {
      
      if (sum(w)==0) {w <- w + 1e-2}
      
      w <- w/sum(w) #normalizing the weight vector
      
      out <- max(w*(Returnscov%*%w))
      return(out)
}

## ---- results="hide", include=FALSE--------------------------------------
# Adjust itermax if it takes too long to compute
itermax = 100
out <- DEoptim(fn=obj, lower=rep(0,10),upper=rep(1,10),DEoptim.control(itermax=itermax,trace=F))
wMinContr <- out$optim$bestmem/sum(out$optim$bestmem)

## ----  include=FALSE-----------------------------------------------------
muPMRC <- wMinContr%*%Returnsmean
sdPMRC <- sqrt(wMinContr%*%Returnscov%*%wMinContr)

## ----fig3, echo=FALSE, fig.cap="\\label{fig:fig3} Effcient frontier (red), line (blue) connecting the risk-free asset and tangency portfolio (blue dot), the minimum variance portfolio (red triangle), the minimum risk contribution portfolio (green symbol) and 10 individual stocks"----
# FINAL plot with all portfolios stocks and efficient frontier etc
plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
points(0, mufree, cex=2, pch=15) 
points(sdP[TAN], muP[TAN], cex=2, pch=16,col='blue') 
points(sdPMRC, muPMRC, cex=3, pch=18,col='green') 
lines(c(0,5), mufree+c(0,5)*(muP[TAN]-mufree)/sdP[TAN], lwd=4, lty=1,
      col="blue") 
for (i in 1:ncol(data))
{
  text(Returnssd[i], Returnsmean[i], colnames(data[,i]), cex=1)
}

## ----  include=FALSE-----------------------------------------------------
# summary stats portfolio
dt <- data.frame(matrix(ncol = 3, nrow = 3))
colnames(dt) <- c("Minimum Variance","Tangency","Minimum Contribution")
rownames(dt) <- c("Expected Return","Standard Deviation","Sharpe Ratio")
MVsharpe <- (muP[MV]-0.03)/sdP[MV]
TANsharpe <- (muP[TAN]-0.03)/sdP[TAN]
MINCsharpe <- (muPMRC-0.03)/sdPMRC
dt[1,] <- c(muP[MV],muP[TAN],muPMRC)
dt[2,] <- c(sdP[MV],sdP[TAN],sdPMRC)
dt[3,] <- c(MVsharpe,TANsharpe,MINCsharpe)

## ---- echo=FALSE---------------------------------------------------------
kable(dt, format = "latex", caption = "Summary Statistics Portfolios", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----  include=FALSE-----------------------------------------------------
# create table for report
dt <- data.frame(matrix(ncol = 10, nrow = 3))
rownames(dt) <- c("Minimum Variance","Tangency","Minimum Contribution")
colnames(dt) <- c("HON","SUN","STW","DRH","AVI","CBL","MER","GPS","MAS","BLC")
dt[1,] <- round(weights[MV,],2)
dt[2,] <- round(weights[TAN,],2)
dt[3,] <- round(wMinContr,2)

## ---- echo=FALSE---------------------------------------------------------
kable(dt, format = "latex", caption = "Portfolios Weights", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----  include=FALSE-----------------------------------------------------
Returns.xts <- xts(completedata*100,order.by = index(completedata)) # multiply by 100 !!

## ----  include=FALSE-----------------------------------------------------
# Rebalancing with 126 obs
Returns.xts <- cbind(.indexwday(Returns.xts),Returns.xts)

for (i in 1:(nrow(Returns.xts)-1)) # loop over daily returns
  { 
      if (Returns.xts[i,1]==4 && Returns.xts[i+1,1]==1) # if after thursday monday thursday ad end of the week
      { 
        Returns.xts[i,1]<-5
      } else {}
     
  }

Returns.xts[,1][!Returns.xts[,1]==5] <- FALSE # all weekdays which are not friday are zero
Returns.xts[,1][Returns.xts[,1]==5] <- TRUE # all weekdays which are friday or thursday last day week are one

# Create new column Friday which is equal to one if friday or thursday if friday does not exist
colnames(Returns.xts) <- c("Friday",colnames(Returns.xts[,-1]))
Returns <- as.data.frame(Returns.xts)

## ----  include=FALSE-----------------------------------------------------
#Create empty dataframes for returns and weights of different time periodes. 
MVweights   <-  matrix(0,nrow = nrow(Returns.xts)-126,ncol=ncol(Returns)-1)
TANweights  <-  matrix(0,nrow = nrow(Returns.xts)-126,ncol=ncol(Returns)-1)
MINCweights <-  matrix(0,nrow = nrow(Returns.xts)-126,ncol=ncol(Returns)-1)
MVweights[1,]   <- weights[MV,]
TANweights[1,]  <- weights[TAN,]
MINCweights[1,] <- wMinContr

## ----  include=FALSE-----------------------------------------------------
#### MV and TAN portfolio
for (a in 2:(nrow(Returns.xts)-127))
{
  
  if (Returns[a+125,]$Friday==TRUE) # rebalance portfolio and calculate return
  {
    Returnscov <- cov(Returns.xts[a:(a+126),-1])
    Returnsmean <- apply(Returns.xts[a:(a+126),-1],2,mean)
    muP <- seq(min(Returnsmean)+.00001, max(Returnsmean)-0.00001, length=nrep)
    
    Amat <- cbind(rep(1,5), Returnsmean[1:10], diag(1,10)) # we reconstruct A matrix

    for (i in 1:nrep) # compute efficient frontier
            {
              bvec <- c(1, muP[i],rep(0,10)) # constraint vector
              result <- solve.QP(Dmat=2*Returnscov[1:10,1:10], dvec=rep(0,10), Amat=Amat, bvec=bvec,
                                 meq=2)  
              # meq=2: all constraints so far are equality constraints
              sdP[i] <- sqrt(result$value)
              weights[i,] <- result$solution
            }
    
    MV <- (sdP == min(sdP)) 
    MVweights[a,] <- t(weights[MV,])
    
    sharpe <- (muP-mufree)/sdP 
    TAN <- (sharpe == max(sharpe))
    TANweights[a,] <- t(weights[TAN,])
  } else 
  {
    MVweights[a,] <- MVweights[a-1,]
    TANweights[a,] <- TANweights[a-1,]
  }
}

## ----  include=FALSE-----------------------------------------------------
#Next, we create a loop for the minimum variance portfolio. 
for (a in 2:(nrow(Returns.xts)-127))
{
  
   if (Returns[a+125,]$Friday==TRUE) # rebalance portfolio and calculate return
  {
     
    Returnscov <- cov(Returns.xts[a:(a+126),-1]) # remove end of the week column
    
        obj <- function(w) {
          
          if (sum(w)==0) {w <- w + 1e-2}
          
          w <- w/sum(w) #normalizing the weight vector
          
          out <- max(w*(Returnscov%*%w))
          return(out)
        }

    out <- DEoptim(fn=obj, lower=rep(0,10),upper=rep(1,10),DEoptim.control(itermax=itermax,trace=F))
    wMinContr <- out$optim$bestmem/sum(out$optim$bestmem)
     
    MINCweights[a,] <- wMinContr

   } else 
  {
    MINCweights[a,] <- MINCweights[a-1,]
  }
}

## ----  include=FALSE-----------------------------------------------------
#Compare the three portfolios by reporting the following statistics of their weekly ex-post returns:
Returnstest <- as.matrix(as.data.frame(t(Returns[127:nrow(Returns),-1])))
colnames(Returnstest) <- seq(1,ncol(Returnstest),by=1)

MVweights <- as.matrix(as.data.frame(MVweights))
rownames(MVweights) <- seq(1,ncol(Returnstest),by=1)

TANweights <- as.matrix(as.data.frame(TANweights))
rownames(TANweights) <- seq(1,ncol(Returnstest),by=1)

MINCweights <- as.matrix(as.data.frame(MINCweights))
rownames(MINCweights) <- seq(1,ncol(Returnstest),by=1)

MVreturns   <- as.data.frame((MVweights%*%Returnstest)[1,])/100
TANreturns  <- as.data.frame((TANweights%*%Returnstest)[1,])/100
MINCreturns <- as.data.frame((MINCweights%*%Returnstest)[1,])/100

## ----  include=FALSE-----------------------------------------------------
### the sample mean of ex-post returns
MVmean<- apply(MVreturns,2,mean)
TANmean<- apply(TANreturns,2,mean)
MINCmean<- apply(MINCreturns,2,mean)
MVmean
TANmean
MINCmean
#Returnsqn <- apply(Returns.xts,2,qn)

## ----  include=FALSE-----------------------------------------------------
### the sample standard deviation of ex-post returns
MVsd<- apply(MVreturns,2,sd)
TANsd<- apply(TANreturns,2,sd)
MINCsd<- apply(MINCreturns,2,sd)
MVsd
TANsd
MINCsd

## ----  include=FALSE-----------------------------------------------------
###  the information ratio IR = mean/std
MVir <- MVmean/MVsd
TANir <- TANmean/TANsd
MINCir <- MINCmean/MINCsd
MVir
TANir
MINCir

## ----  include=FALSE-----------------------------------------------------
###  the 5% and 1% Value-at-Risk: the loss at the (empirical) 5% and 1% quantile of the ex-post returns
MVvar5 <- quantile(MVreturns[,1],0.05)
TANvar5<- quantile(TANreturns[,1],0.05)
MINCvar5<- quantile(MINCreturns[,1],0.05)
MVvar5
TANvar5
MINCvar5

MVvar1<- quantile(MVreturns[,1],0.01)
TANvar1<- quantile(TANreturns[,1],0.01)
MINCvar1<- quantile(MINCreturns[,1],0.01)
MVvar1
TANvar1
MINCvar1

## ----  include=FALSE-----------------------------------------------------
#make xts object of returns and plus one
MVreturns   <- xts((1+(MVreturns)),order.by = index(completedata[127:nrow(completedata)]))
TANreturns  <- xts((1+(TANreturns)),order.by = index(completedata[127:nrow(completedata)]))
MINCreturns <- xts((1+(MINCreturns)),order.by = index(completedata[127:nrow(completedata)]))

## ----  include=FALSE-----------------------------------------------------
# Cumulative returns
MVcumreturns   <- MVreturns
TANcumreturns  <- TANreturns
MINCcumreturns <- MINCreturns

for (i in 2:nrow(MVreturns))
{
  MVcumreturns[i,]    <- MVcumreturns[i-1,]%*%MVreturns[i,]
  TANcumreturns[i,]   <- TANcumreturns[i-1,]%*%TANreturns[i,]
  MINCcumreturns[i,]  <- MINCcumreturns[i-1,]%*%MINCreturns[i,]
}

## ----fig5, echo=FALSE, fig.cap="\\label{fig:fig5}Portfolio wealth evolution starting with an inital wealth of 1000$ with a rolling window with six months"----
# Plot ewealth evolution
PfReturns <- cbind(MVcumreturns,TANcumreturns,MINCcumreturns)*1000
PfReturns <- data.frame(date=index(PfReturns),coredata(PfReturns))
colnames(PfReturns) <- c("Date","MV","TAN","MIN CONTR")
PfReturns <- melt(PfReturns,id="Date")
colnames(msdmelt) <- c("MV","TAN","MIN CONTR")

p <- ggplot(PfReturns, aes(x = Date,y=value)) + geom_line(aes(color=variable)) 
p

## ----  include=FALSE-----------------------------------------------------
####  the asset turnover
MVat   <- MVreturns
TANat  <- TANreturns
MINCat <- MINCreturns

for (i in 2:nrow(MVweights))
{
  MVat[i,] <- sum(abs(MVweights[i,]-MVweights[i-1,]))
  TANat[i,] <- sum(abs(TANweights[i,]-TANweights[i-1,]))
  MINCat[i,] <- sum(abs(MINCweights[i,]-MINCweights[i-1,]))
}

MVat <- subset(MVat, MVat[,1] > 0)
TANat <- subset(TANat, TANat[,1] > 0)
MINCat <- subset(MINCat, MINCat[,1] > 0)

MVat <- mean(MVat[,1])
TANat <- mean(TANat[,1])
MINCat <- mean(MINCat[,1])

## ----  include=FALSE-----------------------------------------------------
# Create table of statistics
dt <- data.frame(matrix(ncol = 3, nrow = 7))
colnames(dt) <- c("Minimum Variance","Tangency","Minimum Contribution")
rownames(dt) <- c("Cumulative Return","Expected Return","Standard Deviation","Information Ratio","Value-at-Risk 5%","Value-at-Risk 1%","Asset Turnover")
dt[1,] <- c(MVcumreturns[nrow(MVcumreturns),],TANcumreturns[nrow(MVcumreturns),],MINCcumreturns[nrow(MVcumreturns),])
dt[2,] <- c(MVmean,TANmean,MINCmean)
dt[3,] <- c(MVsd,TANsd,MINCsd)
dt[4,] <- c(MVir,TANir,MINCir)
dt[5,] <- c(MVvar5,TANvar5,MINCvar5)
dt[6,] <- c(MVvar1,TANvar1,MINCvar1)
dt[7,]<- c(MVat,TANat,MINCat)

## ---- echo=FALSE---------------------------------------------------------
kable(dt, format = "latex", caption = "Summary Statistics Portfolios (rebalancing with six months)", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----  include=FALSE-----------------------------------------------------
# Rolling window with three months
# We do exactly the same expect 63 observations are taken into account 
# Reload data
completedata <- data_backup
data <- data[1:63,]
Returns.xts <- xts(data*100,order.by = index(data)) # multiply by 100 !!

## ----  include=FALSE-----------------------------------------------------
# By the numbers
Returnsmean <- apply(Returns.xts,2,mean)
Returnssd <- apply(Returns.xts,2,sd)
Returnsmedian <- apply(Returns.xts,2,mean)
Returnsqn <- apply(Returns.xts,2,qn) #trade off efficiency and robustness 

## ----  include=FALSE-----------------------------------------------------
Returns <- as.data.frame(data)
Returns$Year <- year(index(data))
Returns$Month <- month(index(data))
Returns <- cbind(row.names(as.data.frame(data)),Returns)
colnames(Returns)<-c("Date",colnames(data),"Year","Month")

Monthlysd <- ddply(Returns,.(Year,Month),colwise(sd))
Monthlysd$Date <- strptime(sprintf("%s-%02d-01",Monthlysd$Year,Monthlysd$Month),"%Y-%m-%d")
msdmelt <- melt(Monthlysd[,-(1:2)],id="Date")
colnames(msdmelt)[2:3] <- c("Equity","std_Return")
ggplot(msdmelt) + geom_line(aes_string(x="Date",y="std_Return",color="Equity"),size=1.5) + ggtitle("Monthly Rolling Standard Deviations of Daily Returns of 10 shares") + theme(panel.background = element_rect(fill="#CCCCEE"), legend.key = element_rect(fill="#CCCCEE"))

## ----  include=FALSE-----------------------------------------------------
# Estimate covariance matrix
Returnscov <- cov(Returns.xts)

## ---- include=FALSE------------------------------------------------------
kable(round(Returnscov,2), format = "latex", caption = "Coviarance Matrix", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----error=TRUE, include=FALSE-------------------------------------------
# No short selling - we add an extra constraint e.g. mutual funds 
Amat <- cbind(rep(1,5), Returnsmean[1:10], diag(1,10),diag(-1,10)) # we reconstruct A matrix

N <- length(Returnsmean) 
Amat <- cbind(rep(1,N), Returnsmean, diag(1,nrow=N), diag(-1,nrow=N)) 

# closure constraint
# diag(1,6) 
nrep <- 1000 # we make 1000 replications
muP <- seq(0.001, 0.06, length=nrep) # nrep targets is a sequence from 0,1% to 6% in 1000 steps
sdP <- rep(0, nrep) # 1000 empty to save standard deviation
weights <- matrix(0, nrow=nrep, ncol=10) # 1000 empty matrices were we will save the weights. 

for (i in 1:nrep)
{
  bvec <- c(1, muP[i], rep(0,N), rep(-1,N)) # constraint vector
  result <- solve.QP(Dmat=2*Returnscov, dvec=rep(0,10), Amat=Amat, bvec=bvec,
                     meq=2)  
  # meq=2: all constraints so far are equality constraints
  sdP[i] <- sqrt(result$value)
  weights[i,] <- result$solution
}
# Crashes, why? 
# Because no short selling constraint limits feasible return within range of indivdual equity returns 
# The return cannot be achieved 

## ----  include=FALSE-----------------------------------------------------
# This will work: because we adjust the returns 
muP <- seq(min(Returnsmean)+0.00001, max(Returnsmean)-0.00001, length=nrep)
sdP <- rep(0, nrep) # 1000 empty to save standard deviation
weights <- matrix(0, nrow=nrep, ncol=10) # 1000 empty matrices were we will save the weights. 

for (i in 1:nrep)
{
  bvec <- c(1, muP[i],rep(0,10), rep(-1,N))  # constraint vector
  result <- solve.QP(Dmat=2*Returnscov[1:10,1:10], dvec=rep(0,10), Amat=Amat, bvec=bvec,
                     meq=2)  
  # meq=2: all constraints so far are equality constraints
  sdP[i] <- sqrt(result$value)
  weights[i,] <- result$solution
}

## ----  include=FALSE-----------------------------------------------------
## Efificient frontier and Minimum Variance Portfolio

plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
# Minumum Variance portfolio
MV <- (sdP == min(sdP)) 
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
# Effective portfolios
EFF <- (muP > muP[MV])
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
weights[MV,]

## ----  include=FALSE-----------------------------------------------------
## Risk free portfolio

mufree <- 0.03 # input value of risk-free interest rate

plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
points(0, mufree, cex=2, pch=15) 

## ----  include=FALSE-----------------------------------------------------
## Tangency Portfolio
sharpe <- (muP-mufree)/sdP 
TAN <- (sharpe == max(sharpe)) 
weights[TAN,] 
sdpT <- sdP[TAN]
mupT <- muP[TAN]

## ----  include=FALSE-----------------------------------------------------
plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
points(0, mufree, cex=2, pch=15) 
points(sdP[TAN], muP[TAN], cex=2, pch=16,col='blue') 
lines(c(0,5), mufree+c(0,5)*(muP[TAN]-mufree)/sdP[TAN], lwd=4, lty=1,
      col="blue") 
for (i in 1:ncol(data))
{
  text(Returnssd[i], Returnsmean[i], colnames(data[,i]), cex=1)
}

## ----  include=FALSE-----------------------------------------------------
## Minimum risk contribution portfolio
library(DEoptim)

obj <- function(w) {
      
      if (sum(w)==0) {w <- w + 1e-2}
      
      w <- w/sum(w) #normalizing the weight vector
      
      out <- max(w*(Returnscov%*%w))
      return(out)
}

## ---- results="hide", include=FALSE--------------------------------------
itermax = 100
out <- DEoptim(fn=obj, lower=rep(0,10),upper=rep(1,10),DEoptim.control(itermax=itermax,trace=F))
wMinContr <- out$optim$bestmem/sum(out$optim$bestmem)

## ----  include=FALSE-----------------------------------------------------
muPMRC <- wMinContr%*%Returnsmean
sdPMRC <- sqrt(wMinContr%*%Returnscov%*%wMinContr)

## ----fig, include=FALSE, fig.cap="\\label{fig:fig} Effcient frontier (red), line (blue) connecting the risk-free asset and tangency portfolio (blue dot), the minimum variance portfolio (red triangle), the minimum risk contribution portfolio (green symbol) and 10 individual stocks"----
plot(sdP, muP, type="l", xlim=c(0,6), ylim=c(-0.4,0.4), lty=3)
points(sdP[MV], muP[MV], cex=2, pch=17, col='red') 
lines(sdP[EFF], muP[EFF], type="l", xlim=c(0,5), ylim=c(0,0.3),
      lwd=3, col="red") 
points(0, mufree, cex=2, pch=15) 
points(sdP[TAN], muP[TAN], cex=2, pch=16,col='blue') 
points(sdPMRC, muPMRC, cex=3, pch=18,col='green') 
lines(c(0,5), mufree+c(0,5)*(muP[TAN]-mufree)/sdP[TAN], lwd=4, lty=1,
      col="blue") 
for (i in 1:ncol(data))
{
  text(Returnssd[i], Returnsmean[i], colnames(data[,i]), cex=1)
}

## ----  include=FALSE-----------------------------------------------------
dt <- data.frame(matrix(ncol = 3, nrow = 3))
colnames(dt) <- c("Minimum Variance","Tangency","Minimum Contribution")
rownames(dt) <- c("Expected Return","Standard Deviation","Sharpe Ratio")
MVsharpe <- (muP[MV]-0.03)/sdP[MV]
TANsharpe <- (muP[TAN]-0.03)/sdP[TAN]
MINCsharpe <- (muPMRC-0.03)/sdPMRC
dt[1,] <- c(muP[MV],muP[TAN],muPMRC)
dt[2,] <- c(sdP[MV],sdP[TAN],sdPMRC)
dt[3,] <- c(MVsharpe,TANsharpe,MINCsharpe)

## ---- include=FALSE------------------------------------------------------
kable(dt, format = "latex", caption = "Summary Statistics Portfolios", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----  include=FALSE-----------------------------------------------------
dt <- data.frame(matrix(ncol = 10, nrow = 3))
rownames(dt) <- c("Minimum Variance","Tangency","Minimum Contribution")
colnames(dt) <- c("HON","SUN","STW","DRH","AVI","CBL","MER","GPS","MAS","BLC")
dt[1,] <- round(weights[MV,],2)
dt[2,] <- round(weights[TAN,],2)
dt[3,] <- round(wMinContr,2)

## ---- include=FALSE------------------------------------------------------
kable(dt, format = "latex", caption = "Portfolios Weights", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

## ----  include=FALSE-----------------------------------------------------
Returns.xts <- xts(completedata*100,order.by = index(completedata)) # multiply by 100 !!

## ----  include=FALSE-----------------------------------------------------
Returns.xts <- cbind(.indexwday(Returns.xts),Returns.xts)

for (i in 1:(nrow(Returns.xts)-1)) # loop over daily returns
  { 
      if (Returns.xts[i,1]==4 && Returns.xts[i+1,1]==1) # if after thursday monday thursday ad end of the week
      { 
        Returns.xts[i,1]<-5
      } else {}
     
  }

Returns.xts[,1][!Returns.xts[,1]==5] <- FALSE # all weekdays which are not friday are zero
Returns.xts[,1][Returns.xts[,1]==5] <- TRUE # all weekdays which are friday or thursday last day week are one

colnames(Returns.xts) <- c("Friday",colnames(Returns.xts[,-1]))
Returns <- as.data.frame(Returns.xts)

## ----  include=FALSE-----------------------------------------------------
#Create empty dataframes for returns and weights of different time periodes. 
MVweights   <-  matrix(0,nrow = nrow(Returns)-63,ncol=ncol(Returns)-1)
TANweights  <-  matrix(0,nrow = nrow(Returns)-63,ncol=ncol(Returns)-1)
MINCweights <-  matrix(0,nrow = nrow(Returns)-63,ncol=ncol(Returns)-1)
MVweights[1,]   <- weights[MV,]
TANweights[1,]  <- weights[TAN,]
MINCweights[1,] <- wMinContr

## ----  include=FALSE-----------------------------------------------------
#### MV and TAN portfolio
for (a in 2:(nrow(Returns.xts)-65))
{
  
  if (Returns[a+63,]$Friday==TRUE) # rebalance portfolio and calculate return
  {
    Returnscov <- cov(Returns.xts[a:(a+63),-1])
    Returnsmean <- apply(Returns.xts[a:(a+63),-1],2,mean)
    muP <- seq(min(Returnsmean)+.00001, max(Returnsmean)-0.00001, length=nrep)
    
    Amat <- cbind(rep(1,5), Returnsmean[1:10], diag(1,10)) # we reconstruct A matrix

    for (i in 1:nrep) # compute efficient frontier
            {
              bvec <- c(1, muP[i],rep(0,10)) # constraint vector
              result <- solve.QP(Dmat=2*Returnscov[1:10,1:10], dvec=rep(0,10), Amat=Amat, bvec=bvec,
                                 meq=2)  
              # meq=2: all constraints so far are equality constraints
              sdP[i] <- sqrt(result$value)
              weights[i,] <- result$solution
            }
    
    MV <- (sdP == min(sdP)) 
    MVweights[a,] <- t(weights[MV,])
    
    sharpe <- (muP-mufree)/sdP 
    TAN <- (sharpe == max(sharpe))
    TANweights[a,] <- t(weights[TAN,])
  } else 
  {
    MVweights[a,] <- MVweights[a-1,]
    TANweights[a,] <- TANweights[a-1,]
  }
}

## ----  include=FALSE-----------------------------------------------------
#Next, we create a loop for the minimum variance portfolio. 
for (a in 2:(nrow(Returns.xts)-65))
{
  
   if (Returns[a+63,]$Friday==TRUE) # rebalance portfolio and calculate return
  {
     
    Returnscov <- cov(Returns.xts[a:(a+63),-1]) # remove end of the week column
    
        obj <- function(w) {
          
          if (sum(w)==0) {w <- w + 1e-2}
          
          w <- w/sum(w) #normalizing the weight vector
          
          out <- max(w*(Returnscov%*%w))
          return(out)
        }

    out <- DEoptim(fn=obj, lower=rep(0,10),upper=rep(1,10),DEoptim.control(itermax=itermax,trace=F))
    wMinContr <- out$optim$bestmem/sum(out$optim$bestmem)
     
    MINCweights[a,] <- wMinContr

   } else 
  {
    MINCweights[a,] <- MINCweights[a-1,]
  }
}

## ----  include=FALSE-----------------------------------------------------
#Compare the three portfolios by reporting the following statistics of their weekly ex-post returns:
Returnstest <- as.matrix(as.data.frame(t(Returns[64:nrow(Returns),-1])))
colnames(Returnstest) <- seq(1,ncol(Returnstest),by=1)

MVweights <- as.matrix(as.data.frame(MVweights))
rownames(MVweights) <- seq(1,ncol(Returnstest),by=1)

TANweights <- as.matrix(as.data.frame(TANweights))
rownames(TANweights) <- seq(1,ncol(Returnstest),by=1)

MINCweights <- as.matrix(as.data.frame(MINCweights))
rownames(MINCweights) <- seq(1,ncol(Returnstest),by=1)

MVreturns   <- as.data.frame((MVweights%*%Returnstest)[1,])/100
TANreturns  <- as.data.frame((TANweights%*%Returnstest)[1,])/100
MINCreturns <- as.data.frame((MINCweights%*%Returnstest)[1,])/100

## ----  include=FALSE-----------------------------------------------------
### the sample mean of ex-post returns
MVmean<- apply(MVreturns,2,mean)
TANmean<- apply(TANreturns,2,mean)
MINCmean<- apply(MINCreturns,2,mean)
MVmean
TANmean
MINCmean
#Returnsqn <- apply(Returns.xts,2,qn)

## ----  include=FALSE-----------------------------------------------------
### the sample standard deviation of ex-post returns
MVsd<- apply(MVreturns,2,sd)
TANsd<- apply(TANreturns,2,sd)
MINCsd<- apply(MINCreturns,2,sd)
MVsd
TANsd
MINCsd

## ----  include=FALSE-----------------------------------------------------
###  the information ratio IR = mean/std

MVir <- MVmean/MVsd
TANir <- TANmean/TANsd
MINCir <- MINCmean/MINCsd
MVir
TANir
MINCir

## ----  include=FALSE-----------------------------------------------------
###  the 5% and 1% Value-at-Risk: the loss at the (empirical) 5% and 1% quantile of the ex-post returns
MVvar5 <- quantile(MVreturns[,1],0.05)
TANvar5<- quantile(TANreturns[,1],0.05)
MINCvar5<- quantile(MINCreturns[,1],0.05)
MVvar5
TANvar5
MINCvar5

MVvar1<- quantile(MVreturns[,1],0.01)
TANvar1<- quantile(TANreturns[,1],0.01)
MINCvar1<- quantile(MINCreturns[,1],0.01)
MVvar1
TANvar1
MINCvar1

## ----  include=FALSE-----------------------------------------------------
#make xts object of returns and plus one
MVreturns   <- xts((1+(MVreturns)),order.by = index(completedata[64:nrow(completedata)]))
TANreturns  <- xts((1+(TANreturns)),order.by = index(completedata[64:nrow(completedata)]))
MINCreturns <- xts((1+(MINCreturns)),order.by = index(completedata[64:nrow(completedata)]))

## ----  include=FALSE-----------------------------------------------------
MVcumreturns   <- MVreturns
TANcumreturns  <- TANreturns
MINCcumreturns <- MINCreturns

for (i in 2:nrow(MVreturns))
{
  MVcumreturns[i,]    <- MVcumreturns[i-1,]%*%MVreturns[i,]
  TANcumreturns[i,]   <- TANcumreturns[i-1,]%*%TANreturns[i,]
  MINCcumreturns[i,]  <- MINCcumreturns[i-1,]%*%MINCreturns[i,]
}

## ----fig10, echo=FALSE, fig.cap="\\label{fig:fig10}Portfolio wealth evolution starting with an inital wealth of 1000$ with a rolling window with three months"----
PfReturns <- cbind(MVcumreturns,TANcumreturns,MINCcumreturns)*1000
PfReturns <- data.frame(date=index(PfReturns),coredata(PfReturns))
colnames(PfReturns) <- c("Date","MV","TAN","MIN CONTR")
PfReturns <- melt(PfReturns,id="Date")
colnames(msdmelt) <- c("MV","TAN","MIN CONTR")

p <- ggplot(PfReturns, aes(x = Date,y=value)) + geom_line(aes(color=variable)) 
p

## ----  include=FALSE-----------------------------------------------------
####  the asset turnover
MVat   <- MVreturns
TANat  <- TANreturns
MINCat <- MINCreturns

for (i in 2:nrow(MVweights))
{
  MVat[i,] <- sum(abs(MVweights[i,]-MVweights[i-1,]))
  TANat[i,] <- sum(abs(TANweights[i,]-TANweights[i-1,]))
  MINCat[i,] <- sum(abs(MINCweights[i,]-MINCweights[i-1,]))
}

MVat <- subset(MVat, MVat[,1] > 0)
TANat <- subset(TANat, TANat[,1] > 0)
MINCat <- subset(MINCat, MINCat[,1] > 0)

MVat <- mean(MVat[,1])
TANat <- mean(TANat[,1])
MINCat <- mean(MINCat[,1])

## ----  include=FALSE-----------------------------------------------------
dt <- data.frame(matrix(ncol = 3, nrow = 7))
colnames(dt) <- c("Minimum Variance","Tangency","Minimum Contribution")
rownames(dt) <- c("Cumulative Return","Expected Return","Standard Deviation","Information Ratio","Value-at-Risk 5%","Value-at-Risk 1%","Asset Turnover")
dt[1,] <- c(MVcumreturns[nrow(MVcumreturns),],TANcumreturns[nrow(MVcumreturns),],MINCcumreturns[nrow(MVcumreturns),])
dt[2,] <- c(MVmean,TANmean,MINCmean)
dt[3,] <- c(MVsd,TANsd,MINCsd)
dt[4,] <- c(MVir,TANir,MINCir)
dt[5,] <- c(MVvar5,TANvar5,MINCvar5)
dt[6,] <- c(MVvar1,TANvar1,MINCvar1)
dt[7,]<- c(MVat,TANat,MINCat)

## ---- echo=FALSE---------------------------------------------------------
kable(dt, format = "latex", caption = "Summary Statistics Portfolios (rebalancing with three months)", booktabs = T) %>% kable_styling(latex_options = c("striped", "hold_position"))

