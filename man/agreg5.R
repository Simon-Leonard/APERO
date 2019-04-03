\name{agreg5}
\alias{agreg5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function used to merge overlapping start regions 
}
\description{
Take detected start regions and mapped read position to merge overlapping start regions
}
\usage{
agreg5(tot2,ff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tot2}{
Start regions data.frame with the following column names : 5end (has to be the first column), PU, str and sumfreq.
    
5end : central position of the start region. 
    
PU : positional uncertainty of the region.
    
str : strand of the region
    
sumfreq : number of starting reads in the region
}
  \item{ff}{
Mapped read coordinates data.frame with the following column names : 5end (has to be the first column), str, sum_freq.
    
5end : 5' positions of mapped reads

str : strain

sum_freq : number of mapped reads with the corresponding 5' position and strand
}
}
\details{
See APERO's publication for more details
}
\value{
A data.frame containing the merged start regions with the same columns than the tot2 data.frame
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
