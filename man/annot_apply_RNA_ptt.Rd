\name{annot_apply_RNA_ptt}
\alias{annot_apply_RNA_ptt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotation function of the APERO package
}
\description{
RNAs or start regions are classified according to their localization regarding to CDS and strand from which they are transcribed.
Orphan sRNAs are transcribed from intergenic regions. 
Transcribed from the same strand of the CDS are: P RNAs/start region corresponding to primary RNAs/start region located in the 250 nucleotides (nt) upstream the CDS without overlapping it; 5’-UTR and 3’-UTR RNAs/start region overlapping start and stop codon respectively. 
Transcribed from the opposite strand of the CDS are:  Ai corresponding to antisense RNAs/start region transcribed inside CDS; Div corresponding to RNAs/start region starting in the 250 nucleotides (nt) upstream the CDS; 5’Ai and 3’Ai RNAs/start region overlapping start and stop codon respectively.
}
\usage{
annot_apply_RNA_ptt(data,ref,genome_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
6 columns data.frame. 1st column = "Position"; 6th column = "Positional.Uncertainty"; One column (2nd, 3rd, 4th or 5th) has to be "str" (for strand)
}
  \item{ref}{
ptt-like annotation data.frame.
Example :

Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
    
190..255	+	21	16127995	thrL	b0001	-	-	thr operon leader peptide
337..2799	+	820	16127996	thrA	b0002	-	COG0527E,COG0527E	Bifunctional aspartokinase/homoserine dehydrogenase 1
2801..3733	+	310	16127997	thrB	b0003	-	COG0083E,COG0083E	homoserine kinase
3734..5020	+	428	16127998	thrC	b0004	-	COG0498E,COG0498E	L-threonine synthase
}
  \item{genome_size}{
Genome size
}
}
\details{
See APERO's publication for more details
}
\value{
Returns an 8 columns table : The input data merged with a "Class" column (short annotation) and a "Comment" column (detailed annotation) 
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
