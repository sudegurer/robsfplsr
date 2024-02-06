fpcr <- function(y, x, nbf = NULL, ncomp = NULL, gp = NULL, nfold=10, CV = TRUE,
                 as = 1:5, Bs = c(4, 5, 8, 10, 15, 20))
{

  if(!is.matrix(x))
    stop("Error!! x must be a matrix!")
  if(!is.null(nbf)){
    if(!nbf >= 4)
      stop("Error!! The number of basis functions (nbasis) must be greater than three to apply cubic B-spline basis expansion!")
  }
  if(!is.null(ncomp)){
    if(!ncomp > 0)
      stop("Error!! The number of principal components (ncomp) must be greater than or equal to one!")
  }
  if(!is.null(gp)){
    if(length(gp) != dim(x)[2])
      stop("Error!! The number of columns of x must be equal to the length of grid points!")
  }

  if(is.null(gp))
    gp <- seq(0, 1, length.out = dim(x)[2])

  if(CV == TRUE){
    optimod <- pcaCV(y = y, x = x, as = as, Bs = Bs, nfold = nfold, gp = gp)
    nbf <- optimod$B
    ncomp <- optimod$a
  }

  pca.results <- getPCA(data = x, nbasis = nbf, ncomp = ncomp, gp = gp)
  x.scores <- pca.results$PCAscore
  model <- lm(y~x.scores)

  coeffs <- as.matrix(model$coefficients)
  fitteds <- cbind(1, x.scores) %*% coeffs
  resids <- y - fitteds

  evalbase <- pca.results$evalbase
  PCAcoef <- pca.results$PCAcoef

  f.coeff = evalbase %*% (PCAcoef$coefs %*% as.matrix(coeffs[-1]))

  model.details <- list()
  model.details$gp <- gp
  model.details$PCAcoef <- pca.results$PCAcoef
  model.details$evalbase <- pca.results$evalbase
  model.details$bs_basis <- pca.results$bs_basis



  return(list(y = y, x = x, f.coeff = f.coeff, fitted.values = fitteds, residuals = resids,
              coeffs = coeffs, model.details = model.details))
}

