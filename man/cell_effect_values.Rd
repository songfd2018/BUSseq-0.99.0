\name{cell_effect_values}
\alias{cell_effect_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the cell-specific size effects from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated cell-specific size effects of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
cell_effect_values(BUSseqfits_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained by the function \code{BUSseq_MCMC}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{delta.est}{The estimated cell-specific size effects, a list with length equal to the number of batch, whose b-th element is a vector of cell-specific effect in the b-th batch. Note that the first element in each vector is zero as the first cell in each batch is taken as the reference cell.}
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
delta.est <- cell_effect_values(BUSseqfits_example)

}
