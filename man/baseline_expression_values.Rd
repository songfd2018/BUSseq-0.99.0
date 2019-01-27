\name{baseline_expression_values}
\alias{baseline_expression_values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the Log-Scale Baseline Expression levels from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated log-scale baseline expression levels from a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
baseline_expression_values(BUSseqfits_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{A BUSseqfits object obtained from the function \code{BUSseq_MCMC}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{alpha.est}{The estimated log-scale baseline expression levels, a G-dimensional vector whose g-th element is the estimated log-scale mean gene expression level of gene g in the first cell type.}
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
alpha.est <- baseline_expression_values(BUSseqfits_example)

}
