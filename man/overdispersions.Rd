\name{overdispersions}
\alias{overdispersions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the Overdispersion Parameters from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated overdispersion parameters of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
overdispersions(BUSseqfits_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{phi.est}{The estimated overdispersion parameters, a B by G matrix, where phi_{bg} is the overdispersion parameter of gene g in batch b.}
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
phi.est <- overdispersions(BUSseqfits_example)

}