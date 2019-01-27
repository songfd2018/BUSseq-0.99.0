\name{summary.BUSseqfits}
\alias{summary.BUSseqfits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a \code{BUSseqfits} Object
}
\description{
Call the function to summarize a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
\method{summary}{BUSseqfits}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.
}
\item{...}{not used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Summarize the results from the output of \code{BUSseq_MCMC}.
}
\references{
Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effect Correction, Manuscript.
}
\author{
Fangda Song
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(BUSseq)

#"BUSseqfits_example" is a sample BUSseqfits object
class(BUSseqfits_example)
summary(BUSseqfits_example)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
