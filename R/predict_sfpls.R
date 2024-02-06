predict_sfpls <- function(object, xnew)
{

  if(!is.matrix(xnew))
    stop("Error!! xnew must be a list!")

  coef <- object$coef
  b0 <- object$intercept
  gp <- object$gp
  predictions <- xnew %*% coef * (gp[2] - gp[1]) + b0

  return(predictions)
}
