\name{compute_alphas}
\alias{compute_alphas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of alphas in the asymptotic covariance matrix of the deflation-based FastICA estimator
}
\description{
Using the estimates of the independent components, the function computes for a given set of nonlinearities the quantities (alphas), which determine the choices of the nonlinearities and in which order they are used in the adaptive deflation-based FastICA method.
}
\usage{
compute_alphas(Z, gs=gf, dgs=dgf, name=gnames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{a numeric matrix of the estimated independent components, which should be standardized so that the mean is zero and the covariance matrix is the identity matrix.}
  \item{gs}{a vector of functions containing the nonlinearities.}
  \item{dgs}{a vector of functions containing the first derivatives of the nonlinearities.}
  \item{name}{a vector of strings containing the names of the nonlinearities.}
}
\details{See the references.}

\value{
A matrix where the ith row gives the estimates of alphas for the ith nonlinearity and the jth column corresponds to the jth component of \code{Z}.}

\references{
\cite{Hyvarinen, A. and Oja, E. (1997), A fast fixed-point algorithm for independent component analysis, \emph{Neural Computation}, vol. 9, 1483--1492.}

\cite{Nordhausen, K., Ilmonen, P., Mandal, A., Oja, H. and Ollila, E. (2011), Deflation-based FastICA reloaded,
in \emph{Proc. "19th European Signal Processing Conference 2011 (EUSIPCO 2011)",} Barcelona, 1854--1858.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014), Deflation-based FastICA with adaptive choices of nonlinearities, \emph{IEEE Transactions on Signal Processing}, 62(21), 5716--5724.}
}
\author{Jari Miettinen}

\seealso{\link{fICA}, \link{nonlinearities}, \link[JADE]{FOBI}, \link[JADE]{k_JADE}}

\examples{
A <- matrix(rnorm(9),3,3)
s1 <- rt(1000,6)
s2 <- rexp(1000,1)
s3 <- runif(1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

Sest <- fICA(X,method="def")$S

compute_alphas(Sest, gs=gf[1:3], dgs=dgf[1:3], name=gnames[1:3])
}

\keyword{multivariate}
