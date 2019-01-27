\name{heatmap_data_BUSseq}
\alias{heatmap_data_BUSseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot the Heatmap of the Log-scale Read Count Data for a \code{CountData} Object
}
\description{
  Plot the heatmap of the log-scale read count data across multiple batches, and then save the resulting images in the user's directory as "png" format. 
}
\usage{
heatmap_data_BUSseq(CountData_obj, gene_set = NULL, project_name="BUSseq_heatmap", 
                          image_dir = NULL, color_key_seq = NULL, 
                          image_width = 1440, image_height = 1080)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CountData_obj}{
    A \code{CountData} object obtained from the function \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
  }
  \item{gene_set}{
    A vector of gene indices indicating the gene set of interest to display in the heatmap. The default is all genes. We also recommend displaying the intrinsic genes only by setting this argument as \code{intrisic_genes_BUSseq(BUSseqfits_obj)}.
  }
  \item{project_name}{
    A string to name the "png" image. By default, the figure is named as "BUSseq_heatmap_log1p_data.png."
  }
  \item{image_dir}{
    A directory to store the gnereated heatmap. The default is to create a folder called "image" in the current directory and save there.
  }
  \item{color_key_seq}{
    A numeric vector indicating the splitting points for binning log-scale read counts into colors. The default is to space the color key points equally between the minimum and maximum of the log-scale read count data.
  }
  \item{image_width}{
    The width in pixels of the graphical device to plot. The default is 1440 px.
  }
  \item{image_height}{
    The height in pixels of the graphical device to plot. The default is 1080 px.
  }
}
\details{
  To cope with the zeros in the count data, we take the transformation log(1+x) on all count data, which corresponds to the R function \code{log1p()} instead of \code{log()}.
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
\dontrun{
  library(BUSseq)
  
  #Disply the read count data for a sample BUSseqfits object
  Example_CountData_raw <- raw_read_counts(BUSseqfits_example)
  class(Example_CountData_raw)
  heatmap_data_BUSseq(Example_CountData_raw)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
