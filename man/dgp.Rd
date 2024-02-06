\name{dgp}
\alias{dgp}
\title{Generate functional data for the scalar-on-function regression model (sparse or non-sparse)}
\description{
This function is used to simulate data for the (sparse and non-sparse) scalar-on-function regression model \deqn{
Y = \beta_0 + \int X(t) \beta(t) dt + \epsilon,} where \eqn{Y} denotes the scalar response, \eqn{X(t)} denotes the functional predictor, \eqn{\beta(t)} denotes the regression coefficient function, and \eqn{\epsilon} is the error process.
}
\usage{
dgp(n, nknots, norder, domain = c(0, 1), snr, simind, out.p = 0)
}
\arguments{
  \item{n}{An integer, specifying the number of observations for each variable to be generated.}
  \item{nknots}{An integer, denoting the number of cubic B-spline basis functions defined on the domain.}
  \item{norder}{An integer, denoting the degree of B-spline basis.}
  \item{domain}{A vector with two elements, denoting the starting and end points of the fine grid. Default is [0,1].}
  \item{snr}{An integer, denoting the signal-to-ratio. Default value is 5.}
  \item{simind}{An integer. If 1, then the non-sparse functional dataset is generated. If 2, then, sparse functional dataset is generated.}
  \item{out.p}{An integer between 0 and 1, denoting the outlier percentage in the generated data.}
}
\details{
In the data generation process, first, the  functional predictor is simulated at 501 equally spaced points within the interval [0,1] based on the following process: \deqn{
X_i(t) = \sum a_{i,j} B_j(t),} where \eqn{B_j(t)} represents cubic B-spline basis functions defined on 50 equally spaced knots over [0, 1], and the coefficients \eqn{a_{ij}} are sampled from a standard normal distribution. If simind = 1, then, the regression coefficient function is generated as follows: \deqn{\beta(t) = 3t e^{t^2} \cos(3 \pi t) + 1}. If simind = 2, then, the regression coefficient function is generated as follows: \deqn{
\beta(t) = \sin(4 \pi t) e^{-10 t^2}} with a zero sub-region occurring after the discrete point 368 within the 501 equally spaced points in the interval [0, 1].
}
\value{
\item{X}{A matricex containing the observations of simulated functional predictor variable.}
\item{Y}{A vector containing the observations of simulated scalar response variable.}
\item{b}{A vector containing the generated regression coefficient function.}
}
\author{
Sude Gurer, Han Lin Shang, Abhijit Mandal, and Ufuk Beyaztas
}
\examples{
library(fda.usc)

nknots <- 50
norder <- 4
snr <- 5
n <- 100
domain <- c(0,1)
simind <- 2

data <- dgp(n=n, nknots=nknots, norder=norder, domain = c(0, 1), snr=snr, simind=simind, out.p = 0)
y  = data$Y
x = data$X
b <- data$b

fx <- fdata(x, argvals = seq(0, 1, length.out = 501))
plot(fx, lty = 1, ylab = "", xlab = "Grid point",
     main = expression(X(t)), mgp = c(2, 0.5, 0), ylim = range(fx))
fb <- fdata(b, argvals = seq(0, 1, length.out = 501))
plot(fb, lty = 1, ylab = "", xlab = "Grid point",
     main = expression(beta(t)), mgp = c(2, 0.5, 0), ylim = range(fb))
}
