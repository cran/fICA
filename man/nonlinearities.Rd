\name{nonlinearities}
\alias{nonlinearities}
\alias{gf}
\alias{dgf}
\alias{gnames}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set of Nonlinearities for Adaptive Deflation-based FastICA Method
}
\description{
The default set of nonlinearities with their first derivatives and names used in \code{\link{adapt_fICA}}.
}
\usage{
gf

dgf

gnames
}
%- maybe also 'usage' for other objects documented here.
\details{ 
 The set of nonlinearities includes both well-known functions (\emph{pow3}, \emph{tanh} and \emph{gaus}) and the ones suggested in Miettinen et al. (2014).
 
  
 The object \code{gf} contains the nonlinearities which are:
\tabular{lll}{ 
gf[[1]] \tab pow3  \tab \eqn{x^3} \cr
gf[[2]] \tab tanh  \tab \eqn{tanh(x)} \cr
gf[[3]] \tab gaus  \tab \eqn{exp(-(x)^2/2)} \cr
gf[[4]] \tab lt0.6 \tab \eqn{(x+0.6)_-^2} \cr
gf[[5]] \tab rt0.6 \tab \eqn{(x-0.6)_+^2} \cr
gf[[6]] \tab bt    \tab \eqn{(x)_+^2-(x)_-^2} \cr
gf[[7]] \tab bt0.2 \tab \eqn{(x-0.2)_+^2-(x+0.2)_-^2} \cr
gf[[8]] \tab bt0.4 \tab \eqn{(x-0.4)_+^2-(x+0.4)_-^2} \cr
gf[[9]] \tab bt0.6 \tab \eqn{(x-0.6)_+^2-(x+0.6)_-^2} \cr
gf[[10]] \tab bt0.8 \tab \eqn{(x-0.8)_+^2-(x+0.8)_-^2} \cr
gf[[11]] \tab bt1.0 \tab \eqn{(x-1.0)_+^2-(x+1.0)_-^2} \cr
gf[[12]] \tab bt1.2 \tab \eqn{(x-1.2)_+^2-(x+1.2)_-^2} \cr
gf[[13]] \tab bt1.4 \tab \eqn{(x-1.4)_+^2-(x+1.4)_-^2} \cr
gf[[14]] \tab bt1.6 \tab \eqn{(x-1.6)_+^2-(x+1.6)_-^2} 
}
The object \code{dgf} the corresponding first derivatives in the same order and gnames the names of them.

For skew sources \emph{lt0.6}  and \emph{rt0.6} combined are more efficient than the commonly used \emph{skew}. The rest of the functions are useful for example for sources with multimodal density functions.  

The user can easily submit a own set or modify the set suggested here. See the example below and the examples in \code{\link{adapt_fICA}}.
}
\references{
\cite{Hyvarinen, A. and Oja, E. (1997), A fast fixed-point algorithm for independent component analysis, \emph{Neural Computation}, vol. 9, 1483--1492.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014), Deflation-based FastICA with adaptive choices of nonlinearities, to appear in \emph{IEEE Transactions on Signal Processing}.}
}
\author{Jari Miettinen}

\seealso{\link{adapt_fICA}}

\examples{
# leaving out functions from the default set and adding a new function  
g <- function(x){x^2}
dg <- function(x){2*x}
gf_new <- c(gf[-c(6,8,10)],g)
dgf_new <- c(dgf[-c(6,8,10)],g)
gnames_new <- c(gnames[-c(6,8,10)],"skew")
}

\keyword{utilities}
