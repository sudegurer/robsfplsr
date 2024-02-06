dgp <- function(n, nknots, norder, domain = c(0, 1), snr, simind, out.p = 0)
{

  if(out.p < 0 | out.p > 1)
    stop("Error!! Outlier percentage must be between 0 and 1!")

  inner.prod <- function(f,basis,j){

    rng <- getbasisrange(basis)
    knots <- c(rng[1],basis$params,rng[2])
    nbasis <- basis$nbasis
    norder <- basis$nbasis - length(knots) + 2

    a <- rng[1]
    if(j-norder > 0) a = knots[j-norder+1]

    b <- rng[2]
    if (j <= nbasis-norder) b = knots[j+1]

    bfun <- function(t){
      mat <- eval.basis(t,basis)
      z <- t(mat[,j])
    }

    y <- integrate(function(t) {f(t)*bfun(t)},a,b)
    y$value
  }

  beta_fun <- function(t, ii){
    if(ii == 1){bf = 3*t+exp(t^2)*cos(3*pi*t)+1}
    else if(ii == 2){bf = sin(4*pi*t)*exp(-10*t^2)}
    else{print("model does not exit")}
    return(bf)
  }

  knots <- seq(domain[1], domain[2], length.out = nknots)
  nbasis <- nknots + norder - 2
  basis <- create.bspline.basis(knots, nbasis, norder)

  cMat1 <- matrix(rnorm(n*nbasis),n,nbasis)
  xfd <- fd(coef=t(cMat1),basisobj=basis)
  tobs    = seq(domain[1], domain[2], length.out = 501)
  x <- t(eval.fd(tobs, xfd))

  G1 <- matrix(0,nbasis,1)
  beta.func <- function(t){beta_fun(t,ii=simind)}
  for(j in 1:nbasis) G1[j] = inner.prod(beta.func,basis,j)
  y0 <- cMat1 %*% G1

  eps0 <- sd(y0)

  y <- y0 + rnorm(n, mean = 0, sd = eps0/sqrt(snr))
  beta0 = apply(as.matrix(tobs),1,beta_fun,ii=simind)

  indx <- 1:n
  if(out.p > 0){
    cMat1.out <- matrix(rchisq(n*nbasis, 2),n,nbasis)
    xfd.out <- fd(coef=t(cMat1.out),basisobj=basis)
    x.out <- t(eval.fd(tobs, xfd.out))

    y.out <- y0 + rnorm(n, mean = 0.75, sd = eps0/sqrt(snr))

    out.indx <- sample(indx, round(n*out.p), replace = F)
    x[out.indx,] <- x.out[out.indx,]
    y[out.indx] <- y.out[out.indx]
  }
  return(list(X = x, Y = y, b = beta0))
}
