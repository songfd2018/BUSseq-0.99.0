\name{celltypes}
\alias{celltypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the Cell-type Indicators from a \code{BUSseqfits} Object
}
\description{
The function gives the cell-type indicators of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
celltypes(BUSseqfits_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object output by the function \code{BUSseq_MCMC}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{w.est}{The estimated cell-type indicators, an R list with length equal to the batch number. The b-th element is the estimated cell-type indicator vector for cells in batch b.}
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
celltypes_est <- celltypes(BUSseqfits_example)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
