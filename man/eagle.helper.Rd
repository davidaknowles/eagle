\name{eagle.helper}
\alias{eagle.helper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function to run the EAGLE procedure. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eagle.helper(alt,n,xFull,xNull,s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alt}{ list (over exonic SNPs) of alternative read counts}
  \item{n}{ list (over exonics SNPs) of total read counts}
  \item{xFull}{list of design matrices for the alternative hypothesis (e.g. including environment)}
  \item{xNull}{list of design matrices for the null hypothesis} 
  \item{s}{list of settings. Create using \code{\link{eagle.settings}}() and then customize. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The returned list is the same as for \code{\link{eagle}}: 
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
%%  ~~who you are~~
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

