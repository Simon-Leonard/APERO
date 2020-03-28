APERO_start_detection <-
function(work_dir=getwd(),bam_name, paired_end_data=TRUE, ptt_file=NULL,wmax,min_dist,enrichment,min_read_number,genome_size){
  #type help(APERO_start_detection) to see function help
  options(digits=20)
  setwd(work_dir)
  
  if (wmax<0){return("Please use a positive value for wmax")}
  if(min_dist<1){return("d has to be higher than 1")}
  if(genome_size<0){return("Genome size has to be a positive value")}
  
  
  ### BAM import and processing----
  a=scanBam(bam_name,
             param=ScanBamParam(what=c("pos","flag","mpos","qname","strand","cigar"), 
                                flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                                 hasUnmappedMate=FALSE,
                                )),
             maxMemory=500,)
    pos=a[[1]]$pos
    flag=a[[1]]$flag
    mpos=a[[1]]$mpos
    qname=a[[1]]$qname
    str=a[[1]]$strand
    cigar=a[[1]]$cigar
    bam=data.frame(pos=pos,mpos=mpos,qname=qname,str=str,flag=flag,cigar=cigar)
    rm(a,pos,flag,mpos,qname,str,cigar)
  
if (paired_end_data==FALSE){
  bam$qwidth=apply(bam,1,qw) ### Use read position and cigar to calculate and return read qwidth
  bam$right_end_read=bam$pos+bam$qwidth-1
  bam$'5end'=ifelse(bam$str=="+", bam$pos, bam$right_end_read)
  # bam$'3end'=ifelse(bam$str=="+", bam$right_end_read, bam$pos)
  
  x=bam[,c("5end","str")]
  rm(bam)
  
  x$agreg=paste(x$'5end',x$str,sep=".")
  y=data.frame(table(x$agreg))
  y=data.frame(colsplit(y$Var1,"\\.",c("5end","str")),sum_freq=y$Freq)
  y=y[,c(1,3,2)]
  colnames(y)=c("5end","sum_freq","str")
  
  d=y[order(y$'5end',decreasing = F),]
  rm(y,x)
  
}else{
  FLAG=c(81,83,97,99,145,147,161,163)
  bam=bam[which(bam$flag%in%FLAG),]
  rm(FLAG)
  
  bam$qwidth=apply(bam,1,qw) ### Use read position and cigar to calculate and return read qwidth
  ### Use Apply because it crashes with snowfall
  if (sum(table(bam$qname)==2)!=length(table(bam$qname))){bam=bam[bam$qname %in% (names(table(bam$qname)[table(bam$qname)==2])),]}
  
  
  ap=bam[which(bam$str=="+"),]# flag=97-99
  an=bam[which(bam$str=="-"),] #flag=81-83
  rm(bam)
  
  ap$right_end_read=ap$pos+ap$qwidth-1
  an$right_end_read=an$pos+an$qwidth-1
  
  p=match(ap$qname,an$qname,nomatch=0)
  ap$right_end_mate=an$right_end_read[p] 
  n=match(an$qname,ap$qname,nomatch=0)
  an$right_end_mate=ap$right_end_read[n]
  
  
  r1=rbind(ap,an)
  rm(ap,an,p,n)
  r1p=r1[which(r1$str=="+"),]
  r1n=r1[which(r1$str=="-"),]
  r1p$lg=r1p$right_end_mate-r1p$pos+1
  r1n$lg=r1n$mpos-r1n$right_end_read-1
  colnames(r1p)=c("5end","mate3end","qname","str","flag","cigar","qwidth","3end","mate5end","TLEN")
  colnames(r1n)=c("3end","mate5end","qname","str","flag","cigar","qwidth","5end","mate3end","TLEN")
  
  r1n=r1n[,c(8,9,3,4,5,6,7,1,2,10)]
  r1=rbind(r1p,r1n)
  r1=r1[which(r1$flag<128),] ### We only keep rows with forward reads
  rm(r1p,r1n)
  
  max_insert_size=500
  
  tp=r1[which(r1$str=="+"),]
  tp=tp[which(tp$'mate3end'>=tp$'5end' & tp$'mate5end'-tp$'5end'<max_insert_size),]
  tn=r1[which(r1$str=="-"),]
  tn=tn[which(tn$'5end'>=tn$'mate3end' & tn$'5end'-tn$'mate5end'<max_insert_size),]
  
  t2=rbind(tp,tn)
  rm(tp,tn)
  
  
  x=t2[,c(1,4)]
  x$agreg=paste(x$'5end',x$str,sep=".")
  y=data.frame(table(x$agreg))
  y=data.frame(colsplit(y$Var1,"\\.",c("5end","str")),sum_freq=y$Freq)
  y=y[,c(1,3,2)]
  colnames(y)=c("5end","sum_freq","str")
  
  d=y[order(y$'5end',decreasing = F),]
  rm(y,x)
}
    
  
  ### Start detection ----
  flanks=c(1:min_dist)
  dem=detection5(d,wmax=(2*wmax)+1,flanks=flanks,seuil_net=enrichment,genome_size=genome_size)
  # Function to detect start regions
  # type help(detection5) to see function help
  
  dem=dem[which(dem$freq>min_read_number),]
  
  ### annotation ----
  if (!is.null(ptt_file)){
    dem$C1=NA
    dem$C2=NA
    dem=dem[,c(1,3:6,2)]

    dem=annot_apply_RNA_ptt(dem,ptt_file,genome_size)
    # Function to annotate start regions
    # type help(annot_apply_RNA_ptt) to see function help
    
    dem$C1=NULL
    dem$C2=NULL
  }
  
  ###----
  print("DONE")
  return(dem)
}
