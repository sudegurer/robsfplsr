fpls_sub <- function(y, x, a, probp1, hampelp2, hampelp3, numit, prec, type)
{
  
  data <- datam <- cbind(y, x)
  data <- as.matrix(data)
  n <- nrow(data)
  q <- ncol(data)
  rnames <- rownames(data)
  rownames(data) <- 1:n
  
  datamc <- daprpr(datam, type = type)
  datac <- attr(datamc,"Center")
  datas <- attr(datamc,"Scale")
  y0 <- datam[,1]
  ys <- datamc[,1]
  ns <-nrow(datamc)
  qs <- ncol(datamc)
  ps <- qs - 1
  
  if(type == "classical"){
    res.nipls <- nipls(data=datamc, a=a)
    b <- coef(res.nipls)
  }else if(type == "robust"){
    zerows <- vector(length=0)
    wx <- sqrt(apply(datamc[,2:qs]^2, 1, sum))
    wx <- wx/median(wx)
    wy <- abs(datamc[,1])
    
    if (length(wy)/2>sum(wy==0)){ # not too many zeros
      wy <- wy/(1.4826*median(wy))
    } else{
      wy <- wy/(1.4826*median(wy[wy!=0]))
    }
    
    probct <- qnorm(probp1)
    hampelb <- qnorm(hampelp2)
    hampelr <- qnorm(hampelp3)
    wx[which(wx <= probct)] <- 1 
    wx[which(wx > probct & wx <= hampelb)] <- probct/abs(wx[which(wx > probct & wx <= hampelb)])
    wx[which(wx > hampelb & wx <= hampelr)] <- probct*(hampelr-abs(wx[which(wx > hampelb & wx <= hampelr)]))/
      (hampelr -hampelb)*1/abs(wx[which(wx > hampelb & wx <= hampelr)])
    wx[which(wx > hampelr)] <- 0
    wy[which(wy <= probct)] <- 1 
    wy[which(wy > probct & wy <= hampelb)] <- probct/abs(wy[which(wy > probct & wy <= hampelb)])
    wy[which(wy > hampelb & wy <= hampelr)] <- probct*(hampelr-abs(wy[which(wy > hampelb & wy <= hampelr)]))/
      (hampelr -hampelb)*1/abs(wy[which(wy > hampelb & wy <= hampelr)])
    wy[which(wy > hampelr)] <- 0 
    
    
    w <- wx * wy
    if(any(w<1e-6)){
      w0 <- which(w<1e-6)
      w <- replace(w,list=w0,values=1e-6)
      we <- w
    } else {
      wxe <- wx
      wye <- wy
      we <- w
    }
    dataw <- as.data.frame(datamc * sqrt(we))
    loops <- 1
    rold <- 10^-5
    difference <- 1
    while ((difference > prec) && loops < numit) {    
      res.nipls <- nipls(data=dataw,a=a)
      yp <- fitted(res.nipls)
      r <- datamc[,1] - yp
      b <- coef(res.nipls)
      Tpls <- res.nipls$scores/sqrt(we)
      if (length(r)/2>sum(r==0)){ 
        r <- abs(r)/(1.4826*median(abs(r))) 
      } else{
        r <- abs(r)/(1.4826*median(abs(r[r!=0])))
      }
      
      dt <- daprpr(Tpls, type = type)
      wtn <- sqrt(apply(dt^2, 1, sum))
      wtn <- wtn/median(wtn)
      
      probct <- qnorm(probp1)
      hampelb <- qnorm(hampelp2)
      hampelr <- qnorm(hampelp3)
      wye <- r
      wye[which(r <= probct)] <- 1
      wye[which(r > probct & r <= hampelb)] <- probct/abs(r[which(r > probct & r <= hampelb)])
      wye[which(r > hampelb & r <= hampelr)] <- probct*(hampelr-abs(r[which(r > hampelb & r <= hampelr)]))/
        (hampelr -hampelb)*1/abs(r[which(r > hampelb & r <= hampelr)])
      wye[which(r > hampelr)] <- 0
      wye <- as.numeric(wye)
      
      probct <- qchisq(probp1,a)
      hampelb <- qchisq(hampelp2, a)
      hampelr <- qchisq(hampelp3, a)
      wte <- wtn
      wte[which(wtn <= probct)] <- 1 
      wte[which(wtn > probct & wtn <= hampelb)] <- probct/abs(wtn[which(wtn > probct & wtn <= hampelb)])
      wte[which(wtn > hampelb & wtn <= hampelr)] <- probct*(hampelr-abs(wtn[which(wtn > hampelb & wtn <= hampelr)]))/
        (hampelr -hampelb)*1/abs(wtn[which(wtn > hampelb & wtn <= hampelr)])
      wte[which(wtn > hampelr)] <- 0
      
      difference <- abs(sum(b^2) - rold)/rold
      rold <- sum(b^2)
      we <- wye * wte
      if(any(we<1e-6)){
        w0 <- which(we<1e-6)
        we <- replace(we,list=w0,values=1e-6)
        zerows <- unique(c(zerows,as.numeric(names(w0))))
      }
      
      if(length(zerows)>=(n/2)){
        break
      }
      dataw <- as.data.frame(datamc * sqrt(we))
      loops <- loops + 1
    }
    if (difference > prec){
      warning(paste("Method did not converge. The scaled difference between norms of the coefficient vectors is ", 
                    round(difference, digits=4)))
    }
    
    w <- we
    w[zerows] <- 0 
    wt <- wte
    wt[zerows] <- 0
    wy <- wye
    wy[zerows] <- 0 
  }
  
  P <- res.nipls$loadings
  W <- res.nipls$W
  R <- res.nipls$R
  Tpls <- scale(datam[,2:qs],center=datac[2:qs],scale=datas[2:qs]) %*% R 

  X0 <- data[,2:q]
  Xs <- daprpr(X0, type = type)
  Xss <- attr(Xs, "Scale")
  coef <- datas[1]/Xss*b
  if(type == "classical"){
    intercept <- mean(data[,1] - data[,2:q]%*%coef)
  }else{
    intercept <- median(data[,1] - data[,2:q]%*%coef)
  }
  
  return(list(bhat = coef, b0hat = intercept, b = b, P = P, W = W, R = R, T = Tpls))
  
}