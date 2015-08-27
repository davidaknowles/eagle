\name{eagle.vem}
\alias{eagle.vem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run the variational EM algorithm of EAGLE. 
}
\description{
Run the variational EM algorithm of EAGLE. The specific design matrices determine whether this is for the null or alternative models. 
}
\usage{
eagle.vem(alt,n,x,settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alt}{ list (over exonic SNPs) of alternative read counts}
  \item{n}{ list (over exonics SNPs) of total read counts}
  \item{x}{list of design matrices for each exonic SNP}
  \item{settings}{list of settings. Create using \code{\link{eagle.settings}}() and then customize. }
}
\details{
}
\value{
  \item{coeffs}{List of regression coefficients at each exonic SNP.}
  \item{mlPerIteration}{(Approximate) marginal likelihood calculated at every EM iteration}
  \item{log.likelihoods}{Per locus marginal likelihoods (used for the likelihood ratio test)}
  \item{rep.shape}{Shape of the variational posterior on the random effect precision at each locus}
  \item{rep.shape}{Rate of the variational posterior on the random effect precision at each locus}
  \item{rep.global.shape}{Shape of the prior on the random effect precisions ($a$ parameter in the manuscript)}
  \item{rep.global.rate}{Rate of the prior on the random effect precisions ($a$ parameter in the manuscript)}
}
\references{
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
