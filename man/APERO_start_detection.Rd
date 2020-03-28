\name{APERO_start_detection}
\alias{APERO_start_detection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
First module of the APERO programm which is dedicated to the identification and mapping of transcriptional start regions, i.e., the 5' ends of putative RNAs.
}
\description{
Inspired by previous algorithms specifically dedicated to Transcript Start Site (TSS) detection, start regions are defined as regions exhibiting a local enrichment in read starts compared to neighbor regions. 
Two common difficulties are taken into account: 1- due to wobbling of the RNA Polymerase, start region peaks are not always defined sharply, but can extend over several nucleotides; 2- adjacent peaks distant of several nucleotides are sometimes observed, reflecting alternate 5' starts of the same transcript and preventing the identification of individual peaks. 

The analysis procedure was designed to handle these issues and optimize the precision, by varying the scanning window sizes, and favoring the most well-defined start region. 
The user can play on two major parameters : the maximal accepted width w of a start region peak, and the minimal distance between independent start regions (i.e., the spatial resolution of the detection method, d). 
Two minors paramets can also be adjusted : The minimum enrichment cutoff (we advise to set 1/2d) and the minimum number of reads needed to detect a start region.
}
\usage{
APERO_start_detection(work_dir = getwd(), bam_name, ptt_file = NULL, wmax, min_dist, enrichment, min_read_number, genome_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{work_dir}{
Path to the bam file used for detection
}
  \item{bam_name}{
Name of the bam file used for detection. It has to contain the ".bam" extension. Example :"my_bam_file.bam"
}
  \item{paired_end_data}{
Is your data paired-end ? Set to FALSE if you have single end data (default : TRUE)
}
  \item{ptt_file}{
Name of the annotation table which is then used to classify each start region with respect to neighboring or overlapping CDS. It has to be already imported in R and respect the ptt format
Here is an example of a required ptt file for RNA annotation:

Escherichia coli str. K-12 substr. MG1655, complete genome. - 1..4641652
    
4140 proteins
    
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
    
190..255	+	21	16127995	thrL	b0001	-	-	thr operon leader peptide
    
337..2799	+	820	16127996	thrA	b0002	-	COG0527E,COG0527E	Bifunctional aspartokinase/homoserine dehydrogenase 1
    
2801..3733	+	310	16127997	thrB	b0003	-	COG0083E,COG0083E	homoserine kinase
    
3734..5020	+	428	16127998	thrC	b0004	-	COG0498E,COG0498E	L-threonine synthase
}
  \item{wmax}{
Maximal accepted width w of a start region peak. It has to be higher or equal to 0
}
  \item{min_dist}{
Minimal distance between independent start regions. it has to be higher than 0 
}
  \item{enrichment}{
Minimum enrichment cutoff to detect a start region
}
  \item{min_read_number}{
Minimum number of reads to detect a start region
}
  \item{genome_size}{
Genome size
}
}
\details{
See APERO's publication for more details
}
\value{
A table representing detected start regions with their positions, strands, window widths, number of reads and annotations (if an annotation table is provided)
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
