\name{sfplsR}
\alias{sfplsR}
\title{Sparse functional partial least squares regression}
\description{This function is used to perform both classical and robust scalar-on-function regression model \deqn{
Y = \beta_0 + \int X(t) \beta(t) dt + \epsilon,} based on sparse functional partial least squares decomposition of the functional predictor.}
\usage{
sfplsR(y, x, gp = NULL, a = NULL, B = NULL, probp1 = 0.95, hampelp2 = 0.975,
hampelp3 = 0.999, numit = 100, prec = 0.01, type = c("classical", "robust"),
nfold=10, CV = TRUE, as = 1:5, etas = c(0, 0.1, 0.3, 0.5, 0.7, 0.9),
Bs = c(4, 5, 8, 10))
}
\arguments{
  \item{y}{A vector containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{x}{A matrix containing the observations of functional predictor \eqn{X(t)}.}
  \item{gp}{A vector containing the grid points of the functional predictor.}
  \item{a}{An integer denoting the number of functional partial least squares components to be computed for the functional predictor.}
  \item{B}{An integer denoting the number of B-spline basis expansion functions. Default is NULL.}
  \item{probp1}{A numeric value used to determine the first outlier cutoff point for the weights.}
  \item{hampelp2}{A numeric value used to determine the first outlier cutoff point for the weights.}
  \item{hampelp3}{A numeric value used to determine the third outlier cutoff point for the weights.}
  \item{numit}{n integer value defining the maximum iteration used to achieve convergence.}
  \item{prec}{A numeric value used for the precision of the coefficient estimate.}
  \item{type}{Method type used to estimate the scalar-on-function linear regression model. Possibilities are "classical" and "robust".}
  \item{nfold}{An integer denoting the number of folds used in the k-fold cross validation. Default value is 10.}
  \item{CV}{Logical. If TRUE, then nfold cross-validation is used to find optimum values of \code{a}, and \code{B}. If FALSE, then the specified \code{a} and \code{B} values are used in the model.}
  \item{as}{A vector containing the candidate elements for the \code{a}.}
  \item{etas}{A vector containing the candidate elements for the sparsity parameter.}
  \item{Bs}{A vector containing the candidate elements for the \code{B}.}
}
\details{If \code{type = "classical"}, then, the sparse NIPALS algorithm of Lee et. al. (2011) is used to obtain functional partial least squares regression components.

If \code{type = "robust"}, then, the robust sparse partial least squares regression algorithm Hoffmann et al. (2015) is used to obtain functional partial least squares regression components.
}
\value{
\item{fitted.values}{A vector containing the residuals.}
\item{coef}{A vector containing estimated regression coefficient function by the functional partial least squares regression.}
\item{intercept}{intercept}{A numeric value containing the estimated intercept.}
\item{residuals}{A vector containing the residuals.}
\item{gp}{A vector containing the grid points.}
}
\references{
D. Lee and W. Lee and Y. Lee and Y.Pawitan (2011), "Sparse partial least-squares regression and its applications to high-throughput data analysis", \emph{Chemometrics and Intelligent Laboratory Systems} \bold{109}(1), 1–8.

I. Hoffmann and S.Serneels and P. Filzmoser and C. Croux (2015), "Sparse partial robust M regression", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{149}, 50–59.
}
\author{
Sude Gurer, Han Lin Shang, Abhijit Mandal, and Ufuk Beyaztas
}
\examples{
nknots <- 50
norder <- 4
snr <- 5
n <- 100
domain <- c(0,1)
simind <- 2

data <- dgp(n=n, nknots=nknots, norder=norder, domain = c(0, 1),
snr=snr, simind=simind, out.p = 0.1)
y  = data$Y
x = data$X
b <- data$b

model.sfplsR <- sfplsR(y, x, type = "classical")
model.RsfplsR <- sfplsR(y, x, type = "robust")
}
