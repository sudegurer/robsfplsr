fplsR <- function(y, x, gp = NULL, a = NULL, B = NULL, probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
                  numit = 100, prec = 0.01, type = c("classical", "robust"), nfold=10, CV = TRUE,
                  as = 1:5, Bs = c(4, 5, 8, 10))
{

  if(!is.matrix(x))
    stop("Error!! x must be a matrix!")
  if(!is.null(gp)){
  if(length(gp) != dim(x)[2])
    stop("Error!! The number of columns of x must be equal to the length of grid points!")
  }
  if(!type %in% c("classical", "robust"))
    stop("Error!! type must be one of the followings: classical or robus !")

  type <- match.arg(type)
  if(is.null(gp))
    gp <- seq(0, 1, length.out = dim(x)[2])

  if(CV == TRUE){
    optimod <- prmsCV(y=y, x=x, as = as, Bs = Bs, nfold = nfold, probp1 = probp1, hampelp2 = hampelp2,
                      hampelp3=hampelp3, numit=numit, prec=prec, gp=gp, type = type)
    a <- optimod$a
    B <- optimod$B
  }

  BS.sol <- getAmat(data = x, nbf = B, gp = gp)
  xfd <- BS.sol$Amat
  evalbase <- BS.sol$evalbase
  sinp_mat <- BS.sol$sinp_mat

  fitmodel <- fpls_sub(y=y, x=xfd, a=a, probp1=probp1, hampelp2=hampelp2,
                       hampelp3=hampelp3, numit=numit, prec=prec, type=type)


  coef <- evalbase %*% t(solve(sinp_mat)) %*% fitmodel$bhat
  fits <- x %*% coef * (gp[2]-gp[1]) + fitmodel$b0hat
  residuals <- y - fits

  return(list(fitted.values = fits, coef = coef, intercept = fitmodel$b0hat,
              residuals = residuals, gp = gp))
}
