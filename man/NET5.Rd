\name{NET5}
\alias{NET5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  NET5 function used in the first module of the APERO programm dedicated to the identification of 5’ end RNAs.
}
\description{
Compute the local enrichment value as described in the paper  
}
\usage{
NET5(test,window_size,flanks_size,genome_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test}{
data.frame with the following column names : 5end, str, sum_freq
    
The 5end column has to be the first column of the data.fame. It contains 5' positions of mapped reads

The str column contains corresponding strand

The sum_freq columns contains the number of mapped reads with the corresponding 5' position and strand
}
  \item{window_size}{
Size of the window to compute local enrichment. Has to be an odd number and is equal to 2w+1 corresponding to the local enrichment formula described in the paper. 
}
  \item{flanks_size}{
Size of the neighboring region to compute local enrichment. Corresponds to f in the local enrichment formula described in the paper.
}
  \item{genome_size}{
Genome size
}
}
\details{
See APERO's publication for more details
}
\value{
The test data.frame with an additional column containing the local enrichment value. 
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

