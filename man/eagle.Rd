\name{eagle}
\alias{eagle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run the full EAGLE procedure. 
}
\description{
The usual usage is to 1) fit the null model across exonic SNPs, while learning the overdispersion hyperparameters, 2) fit the alternative models holding these hyperparameters fixed 3) perform likelihood ratio tests for each locus 
}
\usage{
eagle(alt,n,xFull,xNull,max.its=1000,tol=10.0,debug=F,learn.rev=T,rev=1.0,...
	traceEvery=1,rev.model="global",null.first=T,coeff.reg=0.0,...
	return.aux.variables=F,storeAllCoeffs=F,repRep=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alt}{ list (over exonic SNPs) of alternative read counts}
  \item{n}{ list (over exonics SNPs) of total read counts}
  \item{xFull}{list of design matrices for the alternative hypothesis (e.g. including environment)}
  \item{xNull}{list of design matrices for the null hypothesis} 
       \item{max.its}{maximum iterations of EM to run}
       \item{tol}{consider EM to have converged if the change in the lower bound is less than this}
	    \item{debug}{output debugging information? Also performs additional checks.}
       \item{learn.rev}{whether to learn the overdispersion hyperparameters a,b}
	   \item{rev}{initial random effect variance}
	    \item{traceEvery}{traceEvery: how often to output convergence info to stdout}
       \item{rev.model}{rev.model: one of c("global","regression","local","local.regression"). global: single fixed random effect variance across all exonic SNPs (not recommended). regression: random effect variance is log-linear in the log average read depth (not recommended). local: random effect variance at exonic SNP s is v_s ~ InverseGamma(a,b) [recommended]. local.regression: log(rev)  ~ N(mt+c,v) [gives similar results to "local" but not as well calibrated]}
	     \item{nullFirst}{whether to run the null or alternative model first. Only matters if learning the overdispersion hyperparamters, in which case these are learnt on whichever is run first, and fixed for the second. }
		\item{coeff.reg}{whether to regularize the regression coefficents with a term -coef.reg |beta|^2}
		\item{return.aux.variables}{whether to return the auxiliary variables g (only used for debugging)}
       \item{storeAllCoeffs}{whether to store coefficients throughout the EM algorithm (for debugging)}

       \item{repRep}{1/v in the local.regression overdispersion model log(rev)  ~ N(mt+c,v)}
       \item{rerunFirst}{whether to rerun the first model after learning hyperparameters (not used)}

       \item{learnRepRep}{whether to learn v in the local.regression model}
       \item{learnBetas}{whether to learn the regression coefficients}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The returned list is the same as for \code{\link{eagle.helper}}: 
	\item{p.values}{$p$-values for the likelihood ratio tests. }
	\item{q.values}{Corresponding $q$-values calculated using Benjamini-Hochberg FDR.}
	\item{res.full}{Results for the full (alternative) design matrices. See \code{\link{eagle.vem}} for details.}
	\item{res.null}{Results for the null design matrices.}
	\item{settings}{Settings used, see \code{\link{eagle.settings}}}
	\item{timeFirst}{Computation time for first run.}
	\item{timeSecond}{Computation time for second run.}
	\item{res.firstOld}{If the first run was re-run after fitting hyperparameters, then these are the results of that initial run.}
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

