\name{predict_fpcr}
\alias{predict_fpcr}
\title{Prediction for a scalar-on-function linear regression model based on functional principal component analysis}
\description{This function is used to make prediction for a new set of functional predictor based upon a fitted scalar-on-function linear regression model in the output of \code{\link{fpcr}}.}
\usage{
predict_fpcr(object, xnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{fpcr}}.}
  \item{xnew}{A matrix consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{fpcr}}.}
}
\value{A vector of predicted values of the scalar response variable for the given new functional predictor \code{xnew}.}
\author{Sude Gurer, Han Lin Shang, Abhijit Mandal, and Ufuk Beyaztas}
\examples{
nknots <- 50
norder <- 4
snr <- 5
n <- 100
ntest = 5000
domain <- c(0,1)
simind <- 2

data <- dgp(n=n, nknots=nknots, norder=norder, domain = c(0, 1),
snr=snr, simind=simind, out.p = 0.1)
y  = data$Y
x = data$X

data.test <- dgp(n=ntest, nknots=nknots, norder=norder,
domain = c(0, 1), snr=snr, simind=simind, out.p = 0)
ytest <- data.test$Y
xtest <- data.test$X

model.fpcr <- fpcr(y, x)
predictions.fpcr <- predict_fpcr(object = model.fpcr, xnew = xtest)
}
