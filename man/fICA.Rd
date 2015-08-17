\name{fICA}
\alias{fICA}
\alias{fICA.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Symmetric and Deflation-based FastICA Methods for Independent Component Analysis
}
\description{
The symmetric and deflation-based FastICA methods for the independent component problem. The function estimates the unmixing matrix using a single nonlinearity g.
}
\usage{
fICA(X, g="tanh", dg=NULL, G=NULL, init=NULL, method="sym2", inR=TRUE, 
     eps=1e-06, maxiter=1000)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix. Missing values are not allowed.}
  \item{g}{the nonlinearity, tanh by default, see details.}
  \item{dg}{the first derivative of the nonlinearity, see details.}
  \item{G}{the integral function of the nonlinearity, see details.}
  \item{init}{a numeric matrix for the initial value of the algorithm}
  \item{method}{squared symmetric (\code{"sym2"}), symmetric (\code{"sym"}) or deflation-based (\code{"def"})}.
  \item{inR}{a logical which indicates whether R or C is used for computations, see details.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
}

\details{
The deflation-based FastICA estimate depends on the initial value of the unmixing matrix. This is due to the fact that the algorithm does not always find the global maximum of the objective function, but it may stop to local maxima as well. Hence, the deflation-based FastICA estimate is affine equivariant only if the initial value is affine equivariant. Therefore, we recommend the use of \link{reloaded_fICA} or \link{adapt_fICA}.

The standard nonlinearities can be chosen by \code{g="pow3"}, \code{g="tanh"} or \code{g="gaus"}. These estimates can be computed either in R or in C. In order to use some other nonlinearity, one has to give the nonlinearity and its derivative as functions (and the integral function when squared symmetric FastICA is computed), and the computations have to be performed in R.
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

\cite{Miettinen, J., Nordhausen, K., Oja, H., Taskinen, S. and Virta, J. (2015), The squared symmetric FastICA estimator, submitted. }
}

\author{Jari Miettinen}

\seealso{\link{adapt_fICA}, \link{reloaded_fICA}, \link{nonlinearities}}


\examples{
# creating some toy data
A<- matrix(rnorm(9),3,3)
s1 <- rt(1000,6)
s2 <- rexp(1000,1)
s3 <- runif(1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

# tanh is the default nonlinearity
res1<-fICA(X,method="sym")
coef(res1)
plot(res1) 
require(JADE)
MD(coef(res1),A)

# nonlinearity pow3 is chosen as follows
res2<-fICA(X,g="pow3",method="sym")
coef(res2)
require(JADE)
MD(coef(res2),A)

# nonlinearity from gf is chosen as follows
res3<-fICA(X,g=gf[[6]],dg=dgf[[6]],method="sym")
coef(res3)
require(JADE)
MD(coef(res3),A)

}

\keyword{multivariate}
