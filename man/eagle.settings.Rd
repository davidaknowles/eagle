\name{eagle.settings}
\alias{eagle.settings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default list of settings for EAGLE
}
\description{
Creates a default list of settings for EAGLE
}
\usage{
s <- eagle.settings()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{debug}{output debugging information? Also performs additional checks. }
       \item{max.iterations}{maximum iterations of EM to run}
       \item{convergence.tolerance}{consider EM to have converged if the change in the lower bound is less than this}
       \item{random.effect.variance}{initial random effect variance}
       \item{learn.rev}{whether to learn the overdispersion hyperparameters a,b}
       \item{rev.model}{rev.model: one of c("global","regression","local","local.regression"). global: single fixed random effect variance across all exonic SNPs (not recommended). regression: random effect variance is log-linear in the log average read depth (not recommended). local: random effect variance at exonic SNP s is v_s ~ InverseGamma(a,b) [recommended]. local.regression: log(rev)  ~ N(mt+c,v) [gives similar results to "local" but not as well calibrated]}
       \item{rep.global.shape}{the shape parameter (a in the manuscript) in the default "local" overdispersion prior model}
       \item{rep.global.rate}{the rate parameter (b in the manuscript) in the default "local" overdispersion prior model}
       \item{rep.slope}{m in the local.regression overdispersion model log(rev)  ~ N(mt+c,v)}
       \item{rep.intercept}{c in the local.regression overdispersion model log(rev)  ~ N(mt+c,v)}
       \item{rep.rep}{1/v in the local.regression overdispersion model log(rev)  ~ N(mt+c,v)}
       \item{coeff.regulariser}{whether to regularize the regression coefficents with a term -coef.reg |beta|^2}
       \item{return.aux.variables}{whether to return the auxiliary variables g (only used for debugging)}
       \item{storeAllCoeffs}{whether to store coefficients throughout the EM algorithm (for debugging)}
       \item{null.first}{whether to run the null or alternative model first. Only matters if learning the overdispersion hyperparamters, in which case these are learnt on whichever is run first, and fixed for the second. }
       \item{rerunFirst}{whether to rerun the first model after learning hyperparameters (not used)}
       \item{traceEvery}{traceEvery: how often to output convergence info to stdout}
       \item{learnRepRep}{whether to learn v in the local.regression model}
       \item{learnBetas}{whether to learn the regression coefficients}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
David A. Knowles
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
