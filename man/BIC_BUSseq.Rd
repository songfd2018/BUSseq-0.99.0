\name{BIC_BUSseq}
\alias{BIC_BUSseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain BIC from a \code{BUSseqfits} Object
}
\description{
The function gives the Bayesian Informtion Criterion (BIC) value for a \code{BUSseqfits} object.
}
\usage{
BIC_BUSseq(BUSseqfits_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{BIC_val}{The BIC value for a BUSseqfits object.}
}
\references{
Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effects Correction, Manuscript.
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
Example_BIC <- BIC_BUSseq(BUSseqfits_example)

}