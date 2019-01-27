\name{print.CountData}
\alias{print.CountData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Print the Output by the Function \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}
}
\description{
  The function prints a \code{CountData} Object.
}
\usage{
  \method{print}{CountData}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A \code{CountData} object obtained from the function \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
  }
  \item{\dots}{not used.}
}
%\details{
  %%  ~~ If necessary, more details than the description above ~~
    %}
\value{
  Print the results from the output by the function \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
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
  Example_CountData_raw <- raw_read_counts(BUSseqfits_example)
  print(Example_CountData_raw)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
