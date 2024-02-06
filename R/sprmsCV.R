sprmsCV <- function(y, x, as, etas, Bs, nfold, probp1, hampelp2, hampelp3,
                    numit, prec, gp, type = c("classical", "robust"))
{
  
  upper.trim.mean <- function(x,trim) {
    x <- sort(x) 
    mean(x[1:floor(length(x)*(1-trim))])
  }
  
  
  type <- match.arg(type)
  
  n <- dim(x)[1]
  p <- dim(x)[2]

  folds <- cvFolds(n, K = nfold, R = 1, type = "random")
  CVmat <- as.matrix(expand.grid(as, etas, Bs))
  CVmat <- cbind(CVmat, NA)
  colnames(CVmat) <- c("as", "etas", "Bs", "PE")
  CVmat <- CVmat[CVmat[,3] >= CVmat[,1],]
  
  data <- cbind(y, x)
  
  for(i in 1:dim(CVmat)[1]){
    a <- CVmat[,1][i]
    eta <- CVmat[,2][i]
    B <- CVmat[,3][i]
    prediction_mat <- numeric()
    
    for(f in 1:nfold){
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
      trainmod <- sfpls_sub(y = yf, x = xfd, a = a, eta = eta, probp1 = probp1, hampelp2 = hampelp2,
                            hampelp3=hampelp3, numit=numit, prec=prec, type = type)
      coef <- evalbase %*% t(solve(sinp_mat)) %*% trainmod$bhat
      predictions <- xft %*% coef * (gp[2] - gp[1]) + trainmod$b0hat
      if(type == "classical"){
        prediction_mat[f] <- mean((yft - predictions)^2)
      }else if(type == "robust"){
        prediction_mat[f] <- upper.trim.mean((yft - predictions)^2, trim = 0.2)
      }
    }
    if(type == "classical"){
      CVmat[,4][i] <- mean(sort(prediction_mat)[1:(length(prediction_mat)*0.85)])
    }else if(type == "robust"){
      CVmat[,4][i] <- mean(prediction_mat)
    }
  }
  
  optimaB <- CVmat[which.min(CVmat[,4]),][1:3]
  optima <- optimaB[1]
  optime <- optimaB[2]
  optimB <- optimaB[3]
  
  return(list(a = optima, eta = optime, B = optimB))
}