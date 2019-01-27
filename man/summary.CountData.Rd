\name{summary.CountData}
\alias{summary.CountData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a \code{CountData} Object
}
\description{
Call the function to summarize a \code{CountData} object, the output of \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
}
\usage{
\method{summary}{CountData}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{CountData} object obtained from the function \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
}
\item{...}{not used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Summarize the results from the output of \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
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
summary(Example_CountData_raw)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
