\name{sl-package}
\alias{sl-package}
\docType{package}
\title{
Synthetic likelihood approach to fitting dynamically complex models.
}
\description{ This package contains functions for facilitating inference about ecological dynamic models using the synthetic 
likelihood approach. It also contains example data series, and code implementing model fits. 

Functions \code{\link{nlar}}, \code{\link{sl.acf}} and \code{\link{order.dist}} are designed for fast computation of statistics for use in 
synthetic likelihoods. \code{\link{robust.vcov}} is designed for numerically and statistically robust estimation of the covariance matrix 
for the synthetic likelihood. \code{\link{get.trans}} is for getting transformations of statistics to normality, while 
\code{\link{MVN.check}} is for checking the multivariate normality of statistics approximation. 

Example real analyses are provided in \code{\link{bupalus}} and \code{\link{blowfly}}. {\code{\link{ricker}}}, \code{\link{stork}} 
and \code{\link{logistic}}, provide further examples with simulated data. An example of the construction of a robust version of the 
synthetic likelihood is provided in \code{\link{logistic}}, which also looks at generating from one model and fitting to another. 
\code{\link{stork}} is an example in which demographic stochasticity can cause extinction.

Several models are provided for Nicholson's blowflies, largely at the instigation of referees. \code{\link{blowfly.ll}} implements 
the model from section 1.2.3 of the supplementary material (A better alternative model) of the cited paper. Other models are coded in \code{\link{dsbf.ll}} 
(all stochasticity is demographic --- poor fit) and \code{\link{desbf.ll}} (basically the same as \code{\link{blowfly.ll}}, but with explicit 
demographic stochasticity). To the level of approximation that we are really working here, the model in  \code{\link{blowfly.ll}} is probably
the most reasonable, in my view.

For a complete list of functions type \code{library(help=sl)}.

}
\author{
Simon Wood <simon.wood@r-project.org>

}
\references{
Wood, S.N. (2010) Nicholson's Blowflies Revisited Again: Statistical inference on the edge of chaos.
}
\keyword{ package }
\keyword{models} 



\examples{
## see examples for ricker, logistic and blowfly
}
