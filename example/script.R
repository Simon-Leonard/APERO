# Load Packages :
library(Rsamtools)
library(reshape2)
library(snowfall)

library(APERO)

# Load annotation file
ptt=read.table("E:/APERO/NC_016810.ptt",sep="\t",skip=2,header=T,stringsAsFactors = F)

# 5'end detection
res=APERO_start_detection(work_dir = "E:/APERO/", bam_name = "SRR2149882_position_1_to_10kb.bam", ptt_file =ptt
                          , wmax = 10, min_dist = 10, enrichment = 0.1, min_read_number = 0, genome_size = 4878012)


#3'end detection
res2=APERO_end_detection(work_dir = "E:/APERO/", start_table = res, mTEX_bam = "SRR2149882_position_1_to_10kb.bam"
                         , readthrough_proportion = 0.01, Fmin=NA, thread_number = 8, genome_size = 4878012
                         , ptt_file=ptt)
