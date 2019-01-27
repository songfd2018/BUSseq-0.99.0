\name{location_batch_effects}
\alias{location_batch_effects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the Location Batch Effects from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated location batch effects of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
location_batch_effects(BUSseqfits_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{nu.est}{The estimated location batch effects, a B by G matrix, where [b,g] element is the location batch effect on gene g in the batch b compared with the first batch. Note that the first column is zero as the first batch is taken as the reference batch without batch effects.}
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
nu.est <- location_batch_effects(BUSseqfits_example)

}
