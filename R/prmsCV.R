prmsCV <- function(y, x, as, Bs, nfold, probp1, hampelp2, hampelp3, numit,
                   prec, gp, type = c("classical", "robust"))
{
  
  type <- match.arg(type)
  
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
      
      BS.sol <- getAmat(data = xf, nbf = B, gp = gp)
      evalbase <- BS.sol$evalbase
      sinp_mat <- BS.sol$sinp_mat
      xfd <- BS.sol$Amat
      trainmod <- fpls_sub(y = yf, x = xfd, probp1 = probp1, hampelp2 = hampelp2,
                           hampelp3=hampelp3, numit=numit, prec=prec, a = a, type = type)
      coef <- evalbase %*% t(solve(sinp_mat)) %*% trainmod$bhat
      predictions <- xft %*% coef * (gp[2] - gp[1]) + trainmod$b0hat
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