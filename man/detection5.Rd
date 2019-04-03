\name{detection5}
\alias{detection5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detection of start regions
}
\description{

}
\usage{
 detection5(d,wmax,flanks,seuil_net,genome_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
data.frame with the following column names : 5end, str, sum_freq
The 5end column has to be the first column of the data.fame. It contains 5' positions of mapped reads
The str column contains corresponding strand
The sum_freq columns contains the number of mapped reads with the corresponding 5' position and strand
}
  \item{wmax}{
Maximal size of the window to compute local enrichment. Has to be an odd number and is equal to 2*maximal width+1.
}
  \item{flanks}{
Size of the neighboring region to compute local enrichment. A vector containing multiple sizes can also be.
}
  \item{seuil_net}{
Enrichment threshold to identify start regions
}
  \item{genome_size}{
Genome size
}
}
\details{
See APERO's publication for more details
}
\value{
A start region data.frame described with central position, positionnal uncertainty, strand, and number of starting reads
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Simon Leonard
}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
