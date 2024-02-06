pcaCV <- function(y, x, as, Bs, nfold, gp)
{
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  folds <- cvFolds(n, K = nfold, R = 1, type = "random")
  CVmat <- as.matrix(expand.grid(as, Bs))
  CVmat <- cbind(CVmat, NA)
  colnames(CVmat) <- c("as", "Bs", "PE")
  CVmat <- CVmat[CVmat[,2] >= CVmat[,1],]
  
  data <- cbind(y, x)
  
  for(i in 1:dim(CVmat)[1]){
    a <- CVmat[,1][i]
    B <- CVmat[,2][i]
    prediction_mat <- numeric()
    
    for(f in 1:nfold){
      try({
        dtrain <- data[folds$which!=f,]
        dtest <- data[folds$which==f,]
        
        yf <- dtrain[,1]
        xf <- dtrain[,-1]
        yft <- dtest[,1]
        xft <- dtest[,-1]
        
        trainmod <- fpcr(y=yf, x=xf, nbf = B, ncomp = a, gp = gp, CV = FALSE)
        predictions <- predict_fpcr(object = trainmod, xnew = xft)
        prediction_mat[f] <- mean((yft - predictions)^2)
      }, silent = TRUE)
    }
    CVmat[,3][i] <- mean(sort(prediction_mat)[1:(length(prediction_mat)*0.85)])
  }
  
  optimaB <- CVmat[which.min(CVmat[,3]),][1:2]
  optima <- optimaB[1]
  optimB <- optimaB[2]
  
  return(list(a = optima, B = optimB))
}
