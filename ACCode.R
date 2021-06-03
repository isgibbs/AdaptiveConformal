### This file contains two functions for running adaptive conformal inference in order to reproduce Figures 1, 2, 4, 5, 6, and 7 in https://arxiv.org/abs/2106.00170. 


library(quantreg)
library(rugarch)

### Main method for forming election night predictions of county vote totals as in Figure 2
runElectionNightPred <- function(Y,X,alpha,gamma,tinit = 500,splitSize = 0.75,updateMethod="Simple",momentumBW=0.95){
  T <- length(Y)
  ## Initialize data storage variables
  alphaTrajectory <- rep(alpha,T-tinit)
  adaptErrSeq <-  rep(0,T-tinit)
  noAdaptErrorSeq <-  rep(0,T-tinit)
  alphat <- alpha
  for(t in tinit:T){
    ### Split data into training and calibration set
    trainPoints <- sample(1:(t-1),round(splitSize*(t-1)))
    calpoints <- (1:(t-1))[-trainPoints]
    Xtrain <- X[trainPoints,]
    Ytrain <- Y[trainPoints]
    XCal <- X[calpoints,]
    YCal <- Y[calpoints]
    
    ### Fit quantile regression on training setting
    lqrfitUpper <- rq(Ytrain ~ Xtrain, tau=1-alpha/2)
    lqrfitLower <- rq(Ytrain ~ Xtrain, tau=alpha/2)
    
    ### Compute conformity score on calibration set and on new data example
    predLowForCal <- cbind(rep(1,nrow(XCal)),XCal)%*%coef(lqrfitLower)
    predUpForCal <- cbind(rep(1,nrow(XCal)),XCal)%*%coef(lqrfitUpper)
    scores <- sapply(1:length(YCal),function(x){max(YCal[x]-predUpForCal[x],predLowForCal[x]-YCal[x])})
    qUp <- sum(coef(lqrfitUpper)*c(1,X[t,]))
    qLow <- sum(coef(lqrfitLower)*c(1,X[t,]))
    newScore <- max(Y[t]-qUp,qLow-Y[t])
    
    ## Compute errt for both methods
    confQuantNaive <- quantile(scores,1-alpha)
    noAdaptErrorSeq[t-tinit+1] <- as.numeric(confQuantNaive < newScore)
    
    if(alphat >=1){
      adaptErrSeq[t-tinit+1] <- 1
    }else if (alphat <=0){
      adaptErrSeq[t-tinit+1] <- 0
    }else{
      confQuantAdapt <- quantile(scores,1-alphat)
      adaptErrSeq[t-tinit+1] <- as.numeric(confQuantAdapt < newScore)
    }
    
    ## update alphat
    alphaTrajectory[t-tinit+1] <- alphat
    if(updateMethod=="Simple"){
      alphat <- alphat + gamma*(alpha-adaptErrSeq[t-tinit+1])
    }else if(updateMethod=="Momentum"){
      w <- rev(momentumBW^(1:(t-tinit+1)))
      w <- w/sum(w)
      alphat <- alphat + gamma*(alpha - sum(adaptErrSeq[1:(t-tinit+1)]*w))
    }

    if(t %% 100 == 0){
      print(sprintf("Done %i time steps",t))
    }
  }
  return(list(alphaTrajectory,adaptErrSeq,noAdaptErrorSeq))
}

### Main method for forming volatility predictions as in Figure 1
garchConformalForcasting <- function(returns,alpha,gamma,lookback=1250,garchP=1,garchQ=1,startUp = 100,verbose=FALSE,updateMethod="Simple",momentumBW = 0.95){
  T <- length(returns)
  startUp <- max(startUp,lookback)
  garchSpec <- ugarchspec(mean.model=list(armaOrder = c(0, 0),include.mean=FALSE),variance.model=list(model="sGARCH",garchOrder=c(1,1)),distribution.model="norm")
  alphat <- alpha
  ### Initialize data storage variables
  errSeqOC <- rep(0,T-startUp+1)
  errSeqNC <- rep(0,T-startUp+1)
  alphaSequence <- rep(alpha,T-startUp+1)
  scores <- rep(0,T-startUp+1)
  
  for(t in startUp:T){
    if(verbose){
      print(t)
    }
    ### Fit garch model and compute new conformity score
    garchFit <- ugarchfit(garchSpec, returns[(t-lookback+1):(t-1) ],solver="hybrid")
    sigmaNext <- sigma(ugarchforecast(garchFit,n.ahead=1))
    scores[t-startUp + 1] <- abs(returns[t]^2- sigmaNext^2)/sigmaNext^2
    
    recentScores <- scores[max(t-startUp+1 - lookback + 1,1):(t-startUp)]
    
    ### compute errt for both methods
    errSeqOC[t-startUp+1] <- as.numeric(scores[t-startUp + 1] > quantile(recentScores,1-alphat))
    errSeqNC[t-startUp+1] <- as.numeric(scores[t-startUp + 1] > quantile(recentScores,1-alpha))
    
    ### update alphat
    alphaSequence[t-startUp+1] <- alphat
    if(updateMethod=="Simple"){
      alphat <- alphat + gamma*(alpha - errSeqOC[t-startUp+1])
    }else if(updateMethod=="Momentum"){
      w <- rev(momentumBW^(1:(t-startUp+1)))
      w <- w/sum(w)
      alphat <- alphat + gamma*(alpha - sum(errSeqOC[1:(t-startUp+1)]*w))
    }
    if(t %% 100 == 0){
      print(sprintf("Done %g steps",t))
    }
  }
  
  return(list(alphaSequence,errSeqOC,errSeqNC))
}
