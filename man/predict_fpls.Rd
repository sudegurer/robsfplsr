\name{predict_fpls}
\alias{predict_fpls}
\title{Prediction for a scalar-on-function linear regression model based on both classical and robust functional partial least squares regression}
\description{This function is used to make prediction for a new functional predictor based upon a fitted scalar-on-function linear regression model in the output of \code{\link{fplsR}}.
}
\usage{
predict_fpls(object, xnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{fplsR}}.}
  \item{xnew}{A matrix consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{fplsR}}.}
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

model.fplsR <- fplsR(y, x, type = "classical")
model.RfplsR <- fplsR(y, x, type = "robust")

predictions.fplsR <- predict_fpls(object = model.fplsR, xnew = xtest)
predictions.RfplsR <- predict_fpls(object = model.RfplsR, xnew = xtest)
}
