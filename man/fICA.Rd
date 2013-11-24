\name{fICA}
\alias{fICA}
\alias{fICA.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Symmetric and Deflation-based FastICA Methods for Independent Component Analysis
}
\description{
The symmetric and deflation-based FastICA methods for the independent component problem. The function estimates the unmixing matrix using a single contrast function G.
}
\usage{
fICA(X, g=NULL, dg=NULL, init=NULL, method="sym", eps=1e-06, maxiter=100)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix. Missing values are not allowed.}
  \item{g}{the first derivative of the contrast function G. Also called the nonlinearity.}
  \item{dg}{the second derivative of the contrast function G.}
  \item{init}{a numeric matrix for the initial value of the algorithm. If \code{NULL}, the identity matrix is used.}
  \item{method}{symmetric (\code{"sym"}) or deflation-based (\code{"def"})}.
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
}

\details{If deflation-based fastICA is to be used we recommend the use of \link{adapt_fICA} either with fixed contrast function or adaptively selected contrast functions
since \link{adapt_fICA} is affine equivariant and extracts the sources in an optimal order.

}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{g}{nonlinearity used.}
  \item{method}{symmetric or deflation-based.}
  \item{S}{estimated source components standardized to have mean 0 and unit variances.}
}
\references{
\cite{Hyvarinen, A. and Oja, E. (1997), A fast fixed-point algorithm for independent component analysis, \emph{Neural Computation}, vol. 9, 1483--1492.}
}

\author{Jari Miettinen}

\seealso{\link{adapt_fICA}}



\examples{
# creating some toy data
A<- matrix(rnorm(9),3,3)
s1 <- rt(1000,6)
s2 <- rexp(1000,1)
s3 <- runif(1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

res1<-fICA(X,method="def")
coef(res1)
plot(res1) 
require(JADE)
MD(coef(res1),A)

}

\keyword{multivariate}
