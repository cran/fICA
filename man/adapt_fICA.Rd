\name{adapt_fICA}
\alias{adapt_fICA}
\alias{adapt_fICA.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptive Deflation-based FastICA Method for Independent Component Analysis
}
\description{
The adaptive deflation-based FastICA method for the independent component problem. The function estimates the unmixing matrix by finding the best nonlinearity, for each component separately, from a set of nonlinearities.
}
\usage{
adapt_fICA(X, gs=gf, dgs=dgf, name=gnames, kj=0, inR=TRUE, 
           eps=1e-06, maxiter=100)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix. Missing values are not allowed.}
  \item{gs}{a list of functions containing the nonlinearities.}
  \item{dgs}{a list of functions containing the first derivatives of the nonlinearities.}
  \item{name}{a list of strings containing the names of the nonlinearities.}
  \item{kj}{defines the initial estimate of the unmixing matrix, see details.}
  \item{inR}{a logical which indicates whether R or C is used for computations. If FALSE, the default set of nonlinearities \link{gf} is used.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
}
\details{ The algorithm first finds initial estimates of the sources. The method to find the estimates is decided by the choice of the argument \code{kj}. If the value of \code{kj} is an integer between 1 and number of the sources, then the method is \code{kj}-JADE, otherwise it is FOBI. 
For the meaning of the value \code{kj} used as \code{kj}-JADE, see the help for \code{\link[JADE]{k_JADE}}.

}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{gs}{nonlinearities that were available.}
  \item{used_gs}{nonlinearities, in order of appearance, that were used. The last row of the unmixing matrix follows directly from the other rows, and hence no nonlinearity is connected to it.}
  \item{alphas}{the statistics for the choice of the nonlinearities.}
  \item{init_est}{method that was used for the initial estimate (FOBI or k-JADE).}
  \item{S}{estimated source components standardized to have mean 0 and unit variances.}
}
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

res1<-adapt_fICA(X, inR=FALSE)
res1
coef(res1)
plot(res1) 

require(JADE)
MD(coef(res1),A)

# changing the set of candidate nonlinearities

g <- function(x){x^2}
dg <- function(x){2*x}
gf_new <- c(gf[-c(5,8,10)],g)
dgf_new <- c(dgf[-c(5,8,10)],g)
gnames_new <- c(gnames[-c(5,8,10)],"skew")

res2<-adapt_fICA(X, gs=gf_new, dgs=dgf_new, name=gnames_new)
res2
MD(coef(res2),A)

# reloaded FastICA using tanh

res3<-adapt_fICA(X, gs=gf[2], dgs=dgf[2], name=gnames[2])
res3
MD(coef(res3),A)
}

\keyword{multivariate}
