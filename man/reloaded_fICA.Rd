\name{reloaded_fICA}
\alias{reloaded_fICA}
\alias{reloaded_fICA.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reloaded Deflation-based FastICA Method for Independent Component Analysis
}
\description{
The reloaded and the affine equivariant deflation-based FastICA method for the independent component problem. The function estimates the rows of the unmixing matrix one by one, either in asymptotically optimal order, or so that the objective function is globally maximized at each stage.
}
\usage{
reloaded_fICA(X, g="tanh", dg=NULL, G=NULL, kj=0, method="alpha",
             inR=TRUE, eps=1e-06, maxiter=100)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix. Missing values are not allowed.}
  \item{g}{the nonlinearity, tanh by default, see details.}
  \item{dg}{the first derivative of the nonlinearity, see details.}
  \item{G}{the integral of the nonlinearity, see details.}
  \item{kj}{defines the initial estimate of the unmixing matrix, see details.}
  \item{method}{"alpha" or "G", see details.}
  \item{inR}{a logical which indicates whether R or C is used for computations, see details.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
}
\details{ The algorithm first finds initial estimates of the sources. The method to find the estimates is decided by the choice of the argument \code{kj}. If the value of \code{kj} is an integer between 1 and number of the sources, then the method is \code{kj}-JADE, otherwise it is FOBI. For the meaning of the value \code{kj} used as \code{kj}-JADE, see the help for \code{\link[JADE]{k_JADE}}.

The function has the argument \code{method}, which determines the extraction order of the components. If \code{method="alpha"}, the components are found in asymptotically optimal order, that is, the reloaded deflation-based FastICA estimate is computed. If \code{method="G"}, the objective function is globally maximized at each stage. 

The standard nonlinearities can be chosen by \code{g="pow3"}, \code{g="tanh"} or \code{g="gaus"}. These estimates can be computed either in R or in C. In order to use some other nonlinearity, one has to give the nonlinearity and its derivative (and integral if \code{method="G"})  as functions, and the computations have to be performed in R.
}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{g}{nonlinearity used.}
  \item{alphas}{the statistics for the choice of the nonlinearities.}
  \item{init_est}{method that was used for the initial estimate (FOBI or k-JADE).}
  \item{S}{estimated source components standardized to have mean 0 and unit variances.}
}
\references{
\cite{Hyvarinen, A. and Oja, E. (1997), A fast fixed-point algorithm for independent component analysis, \emph{Neural Computation}, vol. 9, 1483--1492.}

\cite{Nordhausen, K., Ilmonen, P., Mandal, A., Oja, H. and Ollila, E. (2011), Deflation-based FastICA reloaded,
in \emph{Proc. "19th European Signal Processing Conference 2011 (EUSIPCO 2011)",} Barcelona, 1854--1858.}
}
\author{Jari Miettinen}

\seealso{\link{fICA}, \link{nonlinearities}, \link[JADE]{FOBI}, \link[JADE]{k_JADE}}

\examples{
# creating some toy data
A<- matrix(rnorm(9),3,3)
s1 <- rt(1000,6)
s2 <- rexp(1000,1)
s3 <- runif(1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

# tanh is the default nonlinearity
res1<-reloaded_fICA(X)
coef(res1)
plot(res1) 
require(JADE)
MD(coef(res1),A)

# nonlinearity pow3 and method "G" is chosen as follows
res2<-reloaded_fICA(X,g="pow3",method="G")
coef(res2)
require(JADE)
MD(coef(res2),A)

}

\keyword{multivariate}
