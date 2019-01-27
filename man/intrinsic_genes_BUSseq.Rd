\name{intrinsic_genes_BUSseq}
\alias{intrinsic_genes_BUSseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the Intrinsic Gene Indices from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated intrinsic gene indices of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
intrinsic_genes_BUSseq(BUSseqfits_obj, fdr_threshold = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
  \item{fdr_threshold}{The false discovery rate level we want to control. The default is 0.05.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{D.est}{The intrinsic gene indices. The return is a vector.}
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
D.est <- intrinsic_genes_BUSseq(BUSseqfits_example)

}