\name{APERO_end_detection}
\alias{APERO_end_detection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Second module of the APERO program which localizes the 3' end of each RNA.
}
\description{
Starting from an identified Transcript Start Site (TSS), the method consists in iteratively extending the identified RNA by locating the 3' end of the longest fragments (if their number is statistically significant).
This end position can then either be the actual 3' end of the transcript, or a mere intermediate point or spontaneous breaking site within the transcript; to distinguish between these two scenarios, the program computes the number of fragments that overlap this position and were not yet counted in the previous iteration.
If this number is significant (as compared to the number of fragments starting at the 5' end, i.e. the expression strength of this sRNA), the sRNA is further extended toward the 3' end of these new fragments and the operation is repeated; otherwise, the iteration stops and this point is defined as the 3' end.
The tendency of the algorithm to increase the transcript size is set by a threshold parameter (Fmin), which an optimal value can be suggested automatically from the analysis of the dataset provided by the user.
}
\usage{
APERO_end_detection(work_dir = getwd(), start_table, mTEX_bam, pTEX_bam = NA, ptt_file = NULL, readthrough_proportion = 0.01, Fmin = NA, thread_number = 1, genome_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{work_dir}{
Path to the bam file(s) used for detection
}
  \item{start_table}{
Output of the first module of APERO i.e. output of the APERO_start_detection function
}
  \item{mTEX_bam}{
Name of the bam file not treated with TEX
}
  \item{pTEX_bam}{
Name of the bam file treated with TEX. Replace by NA if you do not have TEX treated library (default value)
}
  \item{paired_end_data}{
Is your data paired-end ? Set to FALSE if you have single end data (default : TRUE)
}
  \item{ptt_file}{
Name of the annotation table which is then used to classify each TSS with respect to neighboring or overlapping CDS. It has to be already imported in R and respect the ptt format
Here is an example of a required ptt file for RNA annotation:

Escherichia coli str. K-12 substr. MG1655, complete genome. - 1..4641652
    
4140 proteins
    
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
    
190..255	+	21	16127995	thrL	b0001	-	-	thr operon leader peptide
    
337..2799	+	820	16127996	thrA	b0002	-	COG0527E,COG0527E	Bifunctional aspartokinase/homoserine dehydrogenase 1
    
2801..3733	+	310	16127997	thrB	b0003	-	COG0083E,COG0083E	homoserine kinase
    
3734..5020	+	428	16127998	thrC	b0004	-	COG0498E,COG0498E	L-threonine synthase
}
  \item{readthrough_proportion}{
Proportion of reads considered as readthrough
}
  \item{Fmin}{
Threshold parameter controling the tendency to increase the transcript size. If this value is set to NA, a default value is set automatically (equal to the first quartile of the ratio distribution).
}
  \item{thread_number}{
Number of thread to run the program
}
  \item{genome_size}{
Genome size
}
}
\details{
See APERO's publication for more details
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
