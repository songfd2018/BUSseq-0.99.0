\name{celltype_mean_expression}
\alias{celltype_mean_expression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Obtain the Cell-Type-Specific Mean Expression Levels from a \code{BUSseqfits} Object
}
\description{
  The function gives the estimated cell-type-specific mean expression levels of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
  celltype_mean_expression(BUSseqfits_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
%\details{
  %%  ~~ If necessary, more details than the description above ~~
    %}
\value{
  \item{mu.est}{The estimated cell-type-specific mean expression levels, a G by K matrix, whose [g,k] element is the mean expression levels of cell type k on gene g.}
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
  mu.est <- celltype_mean_expression(BUSseqfits_example)
  
}
