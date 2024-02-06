\name{fpcr}
\alias{fpcr}
\title{Functional principal component regression}
\description{This function is used to perform classical scalar-on-function regression model \deqn{
Y = \beta_0 + \int X(t) \beta(t) dt + \epsilon,} based on functional principal component decomposition of the functional predictor.}
\usage{
fpcr(y, x, nbf = NULL, ncomp = NULL, gp = NULL, nfold=10, CV = TRUE,
as = 1:5, Bs = c(4, 5, 8, 10, 15, 20))
}
\arguments{
  \item{y}{A vector containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{x}{A matrix containing the observations of functional predictor \eqn{X(t)}.}
  \item{nbf}{An integer denoting the number of B-spline basis expansion functions. Default is NULL.}
  \item{ncomp}{An integer denoting the number of functional principal components to be computed for the functional predictor.}
  \item{gp}{A vector containing the grid points of the functional predictor.}
  \item{nfold}{An integer denoting the number of folds used in the k-fold cross validation. Default value is 10.}
  \item{CV}{Logical. If TRUE, then, nfold cross-validation is used to find optimum values of \code{nbf}, and \code{ncomp}. If FALSE, then the specified \code{nbf} and \code{ncomp} values are used in the model.}
  \item{as}{A vector containing the candidate elements for the \code{ncomp}.}
   \item{Bs}{A vector containing the candidate elements for the \code{nbf}.}
}
\details{When performing a scalar-on-function regression model based on the functional principal component analysis, first, the functional predictor \eqn{X(t)} is decomposed by the functional principal component analysis method: \deqn{X(t) = \bar{X}(t) + \sum_{k=1}^{K} a_k \phi_k(t),} where \eqn{\bar{X}(t)} is the mean function, \eqn{\phi_k(t)} is the weight function, and \eqn{a_k = \int (X(k) - \bar{X}(k)) \phi_k(t)} is the principal component score for the functional predictor.

If \code{CV = TRUE}, then, a two-dimensional grid search (\code{as} and \code{Bs} vectors) with 10-fold cross-validation is used to determine optimum values of the number of principal components (\code{ncomp}) and the number of B-spline basis expansion functions (\code{nbf}).
}
\value{
\item{y}{Scalar response variable.}
\item{x}{Functional predictor variable.}
\item{f.coeff}{A vector containing estimated regression coefficient function by the functional principal component regression.}
\item{fitted.values}{A vector containing the residuals.}
\item{coeffs}{A vector containing the estimated model parameters in the finite-dimensional space.}
\item{model.details}{A list object containing model details, such as grid points and functional principal components.}
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

model.fpcr <- fpcr(y, x)
}

