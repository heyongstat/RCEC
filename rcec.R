#----------------------------------------------------------------------------------------
#  Robust Covariance Estimation for Compositional Data
#  Input:
#           x ------ n x p composition data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding 
#        corr ------ correlation estimation based on adaptive thresholding
#----------------------------------------------------------------------------------------
rcec <- function(x, nFoler = 5, soft = 1){
  p <- ncol(x)
  clrX <- log(x) - rowSums(log(x)) %*%matrix(1,1,p) / p
  rcecPred <- RobustAdaptThresoldCov(clrX, soft = soft)
  sigma <- rcecPred$sigma
  corr <- rcecPred$corr
  return(list(sigma = sigma, corr = corr))
}
#------------------------------------------------------------------------------------------------
#  Robust Covariance Estimation for Compositional Data with the Latent Varibale Observable
#  Input:
#           x ------ n x p composition data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding 
#        corr ------ correlation estimation based on adaptive thresholding
#------------------------------------------------------------------------------------------------
oracle <- function(x, nFoler = 5, soft = 1){
  p <- ncol(x)
  rcecPred <- RobustAdaptThresoldCov(x, soft = soft)
  sigma <- rcecPred$sigma
  corr <- rcecPred$corr
  return(list(sigma = sigma, corr = corr))
}
#----------------------------------------------------------------------------------------
#
# Median of Means (MOM) estimate of the covriance matrix
#
#
#----------------------------------------------------------------------------------------
MOM<- function(x){
  n <- nrow(x)
  p <- ncol(x)
  M <- ceiling((2+L)*log(p)) 
  d <- floor(n/M)
  L <- 1
  groupsize <- rep(d,M)+c(rep(1,n-M*d),rep(0,(d+1)*M-n))
  grouplabel <- rep(1:M,groupsize)
  index <- sample(1:n)
  datagroup <- list()
  W <- list()
  mumat <- matrix(0,M,p)
  for(k in 1:M)
  {
    datagroup[[k]] <-(x[index[grouplabel==k],,drop=F])
    W[[k]] <- matrix(0,p,p)
    mumat[k,]=apply(datagroup[[k]],2,mean)
    for (i in 1:p){
      for (j in i:p){
        W[[k]][i,j]<-sum(datagroup[[k]][,i]*datagroup[[k]][,j])/groupsize[k]
        W[[k]][j,i]<-W[[k]][i,j] 
      }
    }  
  }
  mu<-apply(mumat,2,median)
  Gamma<-matrix(apply(sapply(W,c),1,median),p,p)-mu%*%t(mu)
  return(Gamma)
}


#----------------------------------------------------------------------------------------
#  Robust Adaptive thresholding estimation of cov(x)
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding 
#        corr ------ correlation estimation based on adaptive thresholding 
#----------------------------------------------------------------------------------------


RobustAdaptThresoldCov <- function(x, nFolder = 5, soft = 1){
  n <- nrow(x)
  p <- ncol(x)
  # Set the grid for the choice of tuning parameter
  nGrid <- 100
  cov <- MOM(x)
  grid <- 10*rep(1:nGrid)/nGrid
  # Multi-folder cross validation
  part <- 1 + sample(c(1:n))%%nFolder
  error <- matrix(0, nFolder, nGrid)
  for (i in 1:nFolder){
    xTest <- x[which(part == i),]
    xTrain <- x[which(part != i),]
    covTest <- MOM(xTest)
    covTrain <- MOM(xTrain)
    for (j in 1:nGrid){
      sigmaTrain <- RobustAdaptThreshold(covTrain,grid[j],soft,xTrain)
      error[i,j] <- (norm(sigmaTrain-covTest, "F"))
    }
  }
  errorSum <- colSums(error)
  lambda <- grid[which(errorSum == min(errorSum))][1]
  sigma <- RobustAdaptThreshold(cov,lambda,soft,x)
  corr <- diag(diag(sigma)^(-0.5))%*%sigma%*%diag(diag(sigma)^(-0.5))
  return(list(sigma = sigma, corr = corr))
}



#----------------------------------------------------------------------------------------
#  Apply adaptive thresholding to the robust covariance estimate
#  Input:
#           cov ------ p x p covariance matrix
#        lambda ------ tuning parameter
#          soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#         sigma ------ p x p matrix, robust adaptive thresholding result
#----------------------------------------------------------------------------------------

RobustAdaptThreshold <- function(cov,lambda,soft,x){
  n <- nrow(x)
  p <- ncol(x)
  covOffDiag <- cov - diag(diag(cov))
  covDiaghalf <- (diag(cov))^0.5
  sigmaTmp <- abs(covOffDiag) - lambda*((log(p)/n)^0.5)*(covDiaghalf%*%t(covDiaghalf)-diag(diag(cov)))
  sigmaTmp[which(sigmaTmp < 0)] <- 0
  if (soft == 1){
    sigma <- diag(diag(cov)) + sigmaTmp*sign(covOffDiag)
  }else{
    sigma <- cov
    sigma[which(sigmaTmp < 1e-10)] <- 0
    sigma <- sigma + diag(diag(cov))
  }
  return(sigma)
}




