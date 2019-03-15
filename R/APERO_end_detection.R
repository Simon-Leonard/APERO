APERO_end_detection <-
function(work_dir=getwd(),start_table,mTEX_bam,pTEX_bam=NA,ptt_file=NULL,readthrough_proportion=0.01,Fmin=NA,thread_number=1,genome_size){
  options(digits=20)
  setwd(work_dir)
  
  if (thread_number<1){return("The number of thread has to be higher than 1")}
  if(readthrough_proportion<0){return("Readthrough proportion has to be a positive value")}
  if(readthrough_proportion>1){return("Readthrough proportion has to be lower than 1")}
  if(!is.na(Fmin)){if(Fmin<0){return("Fmin has to be a positive value")}}
  if(genome_size<0){return("Genome size has to be a positive value")}
  
  # Functions ----
  couv99_mtex_seul=function(f,pop,frac){
    options(digits=20)
    pu=as.numeric(as.character(f[3]))
    pos=as.numeric(as.character(f[2]))
    lg=as.numeric(as.character(f[8]))
    brin=as.character(f[4])
    if (brin=="+"){
      po=pop[pop$`5end` %in% (floor(pos)-5):(pos+lg-1) & pop$str==brin & pop$mate5end<=(pos+lg-1),]
      f[10]=sum(po$freq)/(max(po$mate5end)-min(po$'5end')+1)
      
      po=pop[pop$`5end` %in% (ceiling(pos)+6):(pos+lg-1) & pop$str==brin & pop$mate5end>(pos+lg-1),]
      
      
      if (sum(po$freq)>0){
        p=rep(po$mate5end,po$freq)
        p=p[order(p,decreasing = T)]
        if(sum(po$freq)==1){
          f[12]=po$mate5end
          f[11]=sum(po$freq)
        }else{
          pp=p[(ceiling(frac*length(p))+1):length(p)]
          f[12]=max(pp)
          f[11]=length(pp)
        }
      }else {
        f[12]=0
        f[11]=0
      }
      
      f[13]=ifelse(nrow(po)>0,sum(po$freq[po$mate5end==f[12]]),0)
      
    }else {
      po=pop[pop$`5end` %in% (ceiling(pos)+5):(pos-lg+1) & pop$str==brin & pop$mate5end>=(pos-lg+1),]
      f[10]=sum(po$freq)/(max(po$'5end')-min(po$mate5end)+1)
      
      po=pop[pop$`5end` %in% (floor(pos)-6):(pos-lg+1) & pop$str==brin & pop$mate5end<(pos-lg+1),]
      
      
      if (sum(po$freq)>0){
        p=rep(po$mate5end,po$freq)
        p=p[order(p,decreasing = F)]
        if(sum(po$freq)==1){
          f[12]=po$mate5end
          f[11]=sum(po$freq)
        }else{
          pp=p[(ceiling(frac*length(p))+1):length(p)]
          f[12]=min(pp)
          f[11]=length(pp)
        }
      }else {
        f[12]=0
        f[11]=0
      }
      
      f[13]=ifelse(nrow(po)>0,sum(po$freq[po$mate5end==f[12]]),0)
    }
    return(f)
  }

  R2_99pour=function(f,pop,tex=c("oui","non"),frac){
    pos=as.numeric(as.character(f[2]))
    pu=as.numeric(as.character(f[3]))
    brin=as.character(f[4])
    if (tex=="oui"){po=pop[pop$`5end` %in% (pos-pu):(pos+pu) & pop$str==brin,]}
    if (tex=="non"){po=pop[pop$`5end` %in% (floor(pos)-5):(ceiling(pos)+5) & pop$str==brin,]}
    
    
    if (sum(po$freq)>0){
      
      p=rep(po$mate5end,po$freq)
      
      if (brin=="+"){
        p=p[order(p,decreasing = T)]
        if(sum(po$freq)==1){f[6]=po$mate5end
        }else{
          pp=p[(ceiling(frac*length(p))+1):length(p)]
          f[6]=max(pp)
        }
        
        
      }else {
        p=p[order(p,decreasing = F)]
        if(sum(po$freq)==1){f[6]=po$mate5end
        }else{
          pp=p[(ceiling(frac*length(p))+1):length(p)]
          f[6]=min(pp)
        }
        
      }
      
    }
    
    
    return(f)
  }
  
  recursive3=function(tab,h,pop,frac,nb_coeur,seuil_ratio,dem_TEX=c("OUI","NON")){
    T1<-Sys.time()
    h1=tab[tab$rap>seuil_ratio,] # 3'end extension
    if(nrow(h1)==0){print("Nothing else to extend");return(list(NULL,tab))}
    h2=tab[tab$rap<=seuil_ratio,] # No 3'end extension
    
    if(dem_TEX=="OUI"){
      max_freq_dem_int=apply(h1,1,max_dem_int,dem=h)
      h2b=h1[max_freq_dem_int>h1$sumfreq,] # No 3'end extension due to the presence of a bigger start between the 5' and the 3'end
      h1=h1[max_freq_dem_int<=h1$sumfreq,] # 3'end extension
      if(nrow(h1)==0){print("Nothing to extend in the input");return(list(NULL,tab))}
    }else{h2b=NULL}
    
    
    h1$lg_agreg=ifelse(h1$str=="+",h1$fin_couple_chev-h1$Position+1,h1$Position-h1$fin_couple_chev+1)
    
    h1$nb_int=0
    h1$nb_chev=0
    h1$fin_couple_chev=0
    h1$nb_fin_chev=0
    
    if(nrow(h1)>20){
      library(snowfall)
      cl=makeCluster(nb_coeur, type = "SOCK")
      vglobal=list(pop)
      clusterExport(cl, "vglobal", envir = environment())
      da=h1
      res<-parApply(cl, MARGIN = 1,X = da, FUN = couv99_mtex_seul,pop=pop, frac)
      stopCluster(cl)
      
      
      res2=data.frame(t(res))
    }else{res2=data.frame(t(apply(X = h1,MARGIN = 1,FUN = couv99_mtex_seul,pop=pop,frac)))}
    
    
    
    
    res2$Position=as.numeric(as.character(res2$Position))
    res2$fin_couple_chev=as.numeric(as.character(res2$fin_couple_chev))
    res2$Position=as.numeric(as.character(res2$Position))
    res2$nb_int=as.numeric(as.character(res2$nb_int))
    res2$nb_chev=as.numeric(as.character(res2$nb_chev))
    res2$rap=as.numeric(as.character(res2$rap))
    res2$sumfreq=as.numeric(as.character(res2$sumfreq))
    
    res2$rap=res2$nb_chev/res2$nb_int
    T2<-Sys.time() 
    print(difftime(T2,T1))
    return(list(res2,rbind(h2,h2b)))
  }

  
  # -----
  nb_coeur=round(thread_number) ### Number of threads
  frac=readthrough_proportion
  
  d = start_table
  MTEX=mTEX_bam
  PTEX=pTEX_bam
  
  
  d$ID_Transcrit=1:nrow(d)
  d=data.frame(d$ID_Transcrit,d$Position,d$Positional.Uncertainty,d$str,d$freq)
  colnames(d)=gsub("d.","",colnames(d))
  colnames(d)[3]="PU"
  d$R2_99pour=0
  
  
  
  
  if (!is.na(PTEX)){
    {a=scanBam(MTEX,
               param=ScanBamParam(what=c("pos","flag","mpos","qname","strand","cigar"), 
                                  flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                                   hasUnmappedMate=FALSE,
                                  )),
               maxMemory=500)
    pos=a[[1]]$pos
    flag=a[[1]]$flag
    mpos=a[[1]]$mpos
    qname=a[[1]]$qname
    str=a[[1]]$strand
    cigar=a[[1]]$cigar
    bam=data.frame(pos=pos,mpos=mpos,qname=qname,str=str,flag=flag,cigar=cigar)
    rm(a,pos,flag,mpos,qname,str,cigar)
    
    FLAG=c(81,83,97,99,145,147,161,163)
    bam=bam[which(bam$flag%in%FLAG),]
    rm(FLAG)
    
    bam$qwidth=apply(bam,1,qw) ### Apply with snowfall crash
    if (sum(table(bam$qname)==2)!=length(table(bam$qname))){bam=bam[bam$qname %in% (names(table(bam$qname)[table(bam$qname)==2])),]}
    
    
    ap=bam[which(bam$str=="+"),]# flag=97-99
    an=bam[which(bam$str=="-"),] # flag=81-83
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
    r1=r1[which(r1$flag<128),] ### only keep forward reads
    rm(r1p,r1n)
    
    max_insert_size=500
    
    tp=r1[which(r1$str=="+"),]
    tp=tp[which(tp$'mate3end'>=tp$'5end' & tp$'mate5end'-tp$'5end'<max_insert_size),]
    tn=r1[which(r1$str=="-"),]
    tn=tn[which(tn$'5end'>=tn$'mate3end' & tn$'5end'-tn$'mate5end'<max_insert_size),]
    
    t2=rbind(tp,tn)
    rm(tp,tn)
    
    # table faster than aggregate
    x=t2[,c(1,4,9)]
    x$agreg=paste(x$'5end',x$mate5end,x$str,sep=".")
    y=data.frame(table(x$agreg))
    y=data.frame(colsplit(y$Var1,"\\.",c("5end","mate5end","str")),freq=y$Freq)
    colnames(y)[1]="5end"
    y=y[,c(1,2,4,3)]
    
    
    y=y[order(y$"5end",decreasing = F),]
    
    MTEX=y
    
    rm(y,x)
    
    }
    {a=scanBam(PTEX,
               param=ScanBamParam(what=c("pos","flag","mpos","qname","strand","cigar"), 
                                  flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                                   hasUnmappedMate=FALSE,
                                  )),
               maxMemory=500)
      pos=a[[1]]$pos
      flag=a[[1]]$flag
      mpos=a[[1]]$mpos
      qname=a[[1]]$qname
      str=a[[1]]$strand
      cigar=a[[1]]$cigar
      bam=data.frame(pos=pos,mpos=mpos,qname=qname,str=str,flag=flag,cigar=cigar)
      rm(a,pos,flag,mpos,qname,str,cigar)
      
      FLAG=c(81,83,97,99,145,147,161,163)
      bam=bam[which(bam$flag%in%FLAG),]
      rm(FLAG)
      
      bam$qwidth=apply(bam,1,qw) ### Apply with snowfall crash
      if (sum(table(bam$qname)==2)!=length(table(bam$qname))){bam=bam[bam$qname %in% (names(table(bam$qname)[table(bam$qname)==2])),]}
      
      
      ap=bam[which(bam$str=="+"),]#flag=97-99
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
      r1=r1[which(r1$flag<128),] ### only keep forward reads
      rm(r1p,r1n)
      
      max_insert_size=500
      
      tp=r1[which(r1$str=="+"),]
      tp=tp[which(tp$'mate3end'>=tp$'5end' & tp$'mate5end'-tp$'5end'<max_insert_size),]
      tn=r1[which(r1$str=="-"),]
      tn=tn[which(tn$'5end'>=tn$'mate3end' & tn$'5end'-tn$'mate5end'<max_insert_size),]
      
      t2=rbind(tp,tn)
      rm(tp,tn)
      
      # table faster than aggregate
      x=t2[,c(1,4,9)]
      x$agreg=paste(x$'5end',x$mate5end,x$str,sep=".")
      y=data.frame(table(x$agreg))
      y=data.frame(colsplit(y$Var1,"\\.",c("5end","mate5end","str")),freq=y$Freq)
      colnames(y)[1]="5end"
      y=y[,c(1,2,4,3)]
      
      
      y=y[order(y$"5end",decreasing = F),]
      
      PTEX=y
      
      rm(y,x)
      
    }
    
    
    
    dm=d
    dp=d
    
    ###### MINUS TEX
    dm$R2_99pour=0
    T1<-Sys.time()
    cl=makeCluster(nb_coeur, type = "SOCK")
    vglobal<-list(MTEX)
    clusterExport(cl, "vglobal", envir = environment())
    res<-parApply(cl, MARGIN = 1,X = dm, FUN = R2_99pour,pop=MTEX,tex="non",frac)
    stopCluster(cl)
    T2<-Sys.time() 
    difftime(T2,T1)
    res2=data.frame(t(res))
    res2$ID_Transcrit=(as.character(res2$ID_Transcrit))
    res2$Position=as.numeric(as.character(res2$Position))
    res2$PU=as.numeric(as.character(res2$PU))
    res2$freq=as.numeric(as.character(res2$freq))
    res2$R2_99pour=as.numeric(as.character(res2$R2_99pour))
    res2$str=(as.character(res2$str))
    dm=res2
    dm$R2_99pour[dm$R2_99pour==0]=dm$Position[dm$R2_99pour==0]
    
    ###### PLUS TEX
    dp$R2_99pour=0
    T1<-Sys.time()
    cl=makeCluster(nb_coeur, type = "SOCK")
    vglobal<-list(MTEX)
    clusterExport(cl, "vglobal", envir = environment())
    res<-parApply(cl, MARGIN = 1,X = dp, FUN = R2_99pour,pop=PTEX,tex="oui",frac)
    stopCluster(cl)
    T2<-Sys.time() 
    difftime(T2,T1)
    res2=data.frame(t(res))
    res2$ID_Transcrit=(as.character(res2$ID_Transcrit))
    res2$Position=as.numeric(as.character(res2$Position))
    res2$PU=as.numeric(as.character(res2$PU))
    res2$freq=as.numeric(as.character(res2$freq))
    res2$R2_99pour=as.numeric(as.character(res2$R2_99pour))
    res2$str=(as.character(res2$str))
    dp=res2
    
    res2=d
    
    rp=res2[res2$str=="+",]
    rn=res2[res2$str=="-",]
    dpp=dp[dp$str=="+",]
    dpn=dp[dp$str=="-",]
    dmp=dm[dm$str=="+",]
    dmn=dm[dm$str=="-",]
    
    rp$R2_99pour=pmax(dpp$R2_99pour,dmp$R2_99pour)
    rn$R2_99pour=pmin(dpn$R2_99pour,dmn$R2_99pour)
    
    res2=rbind(rp,rn)
    res2=res2[order(res2$Position,decreasing = F),]
    
    rm(d,dp,PTEX,res,dm,dpp,dpn,dmp,dmn,rp,rn)
    
    m=res2
    rm(res2)
    
    m$lg_R2_99=ifelse(m$str=="+",m$R2_99pour-m$Position+1,m$Position-m$R2_99pour+1)
    m$lg_agreg=m$lg_R2_99
    m$sumfreq=apply(m,1,sumfreq)
    
    f=m
    rm(m)
    
    f$nb_int=0
    f$nb_chev=0
    f$fin_couple_chev=0
    f$nb_fin_chev=0
    
    T1<-Sys.time()
    cl=makeCluster(nb_coeur, type = "SOCK")
    vglobal<-list(MTEX)
    clusterExport(cl, "vglobal", envir = environment())
    res<-parApply(cl, MARGIN = 1,X = f, FUN = couv99_mtex_seul,pop=MTEX, frac)
    stopCluster(cl)
    T2<-Sys.time() 
    difftime(T2,T1)
    
    h=data.frame(t(res))
    rm(f,res)
    h$ID_Transcrit=(as.character(h$ID_Transcrit))
    h$Position=as.numeric(as.character(h$Position))
    h$PU=as.numeric(as.character(h$PU))
    h$str=(as.character(h$str))
    h$freq=(as.character(h$freq))
    h$R2_99pour=as.numeric(as.character(h$R2_99pour))
    h$lg_R2_99=as.numeric(as.character(h$lg_R2_99))
    h$lg_agreg=as.numeric(as.character(h$lg_agreg))
    h$sumfreq=as.numeric(as.character(h$sumfreq))
    h$nb_int=as.numeric(as.character(h$nb_int))
    h$nb_chev=as.numeric(as.character(h$nb_chev))
    h$fin_couple_chev=as.numeric(as.character(h$fin_couple_chev))
    h$nb_fin_chev=as.numeric(as.character(h$nb_fin_chev))
    
    h$rap=h$nb_chev/h$nb_int
    seuil_ratio_chev_int=ifelse(is.na(Fmin),summary(h$rap)[2],Fmin)
    
    print(paste("Fmin =",seuil_ratio_chev_int))
    
    vglobal=list(MTEX)
    list1=recursive3(h,h,pop=MTEX,frac,nb_coeur,seuil_ratio=seuil_ratio_chev_int,dem_TEX="OUI")
    list1[[2]]$nb_allong=0
    
    fin=list1[[2]]
    compt=1
    while(is.null(list1[[1]])==FALSE){
      list1=recursive3(list1[[1]],h,MTEX,frac,nb_coeur,seuil_ratio=seuil_ratio_chev_int,dem_TEX="OUI")
      
      if(nrow(list1[[2]])>0){
        list1[[2]]$nb_allong=compt
        fin=rbind(fin,list1[[2]])
      }
      compt=compt+1 
    }
    
    
    fi=fin[,-(6:7)]
    fi$lg=as.numeric(as.character(fi$lg_agreg))
    fi$Position=as.numeric(as.character(fi$Position))
    fi$PU=as.numeric(as.character(fi$PU))
    fi$ID_Transcrit=as.character(fi$ID_Transcrit)
    fi$freq=(as.character(fi$freq))
    
    fi=fi[order(fi$Position,decreasing = F),]
    
  }else{
    rm(PTEX)
    
    {a=scanBam(MTEX,
               param=ScanBamParam(what=c("pos","flag","mpos","qname","strand","cigar"), 
                                  flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                                   hasUnmappedMate=FALSE,
                                  )),
               maxMemory=500)
      pos=a[[1]]$pos
      flag=a[[1]]$flag
      mpos=a[[1]]$mpos
      qname=a[[1]]$qname
      str=a[[1]]$strand
      cigar=a[[1]]$cigar
      bam=data.frame(pos=pos,mpos=mpos,qname=qname,str=str,flag=flag,cigar=cigar)
      rm(a,pos,flag,mpos,qname,str,cigar)
      
      FLAG=c(81,83,97,99,145,147,161,163)
      bam=bam[which(bam$flag%in%FLAG),]
      rm(FLAG)
      
      bam$qwidth=apply(bam,1,qw) ### Apply with snowfall crash
      if (sum(table(bam$qname)==2)!=length(table(bam$qname))){bam=bam[bam$qname %in% (names(table(bam$qname)[table(bam$qname)==2])),]}
      
      
      ap=bam[which(bam$str=="+"),]# flag=97-99
      an=bam[which(bam$str=="-"),] # flag=81-83
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
      r1=r1[which(r1$flag<128),] ### only keep forward reads
      rm(r1p,r1n)
      
      max_insert_size=500
      
      tp=r1[which(r1$str=="+"),]
      tp=tp[which(tp$'mate3end'>=tp$'5end' & tp$'mate5end'-tp$'5end'<max_insert_size),]
      tn=r1[which(r1$str=="-"),]
      tn=tn[which(tn$'5end'>=tn$'mate3end' & tn$'5end'-tn$'mate5end'<max_insert_size),]
      
      t2=rbind(tp,tn)
      rm(tp,tn)
      
      # table faster than aggregate
      x=t2[,c(1,4,9)]
      x$agreg=paste(x$'5end',x$mate5end,x$str,sep=".")
      y=data.frame(table(x$agreg))
      y=data.frame(colsplit(y$Var1,"\\.",c("5end","mate5end","str")),freq=y$Freq)
      colnames(y)[1]="5end"
      y=y[,c(1,2,4,3)]
      
      
      y=y[order(y$"5end",decreasing = F),]
      
      MTEX=y
      
      rm(y,x)
      
    }
    
    
    d$R2_99pour=0
    ###### MINUS TEX
    T1<-Sys.time()
    cl=makeCluster(nb_coeur, type = "SOCK")
    vglobal<-list(MTEX)
    clusterExport(cl, "vglobal", envir = environment())
    res<-parApply(cl, MARGIN = 1,X = d, FUN = R2_99pour,pop=MTEX,tex="non",frac)
    stopCluster(cl)
    T2<-Sys.time() 
    difftime(T2,T1)
    
    res2=data.frame(t(res))
    rm(d,res)
    
    res2$ID_Transcrit=(as.character(res2$ID_Transcrit))
    res2$Position=as.numeric(as.character(res2$Position))
    res2$PU=as.numeric(as.character(res2$PU))
    res2$freq=as.numeric(as.character(res2$freq))
    res2$R2_99pour=as.numeric(as.character(res2$R2_99pour))
    res2$str=(as.character(res2$str))
    
    
    m=res2
    rm(res2)
    
    m$lg_R2_99=ifelse(m$str=="+",m$R2_99pour-m$Position+1,m$Position-m$R2_99pour+1)
    m$lg_agreg=m$lg_R2_99
    m$sumfreq=apply(m,1,sumfreq)
    

    f=m
    rm(m)
    f$nb_int=0
    f$nb_chev=0
    f$fin_couple_chev=0
    f$nb_fin_chev=0
    
    T1<-Sys.time()
    cl=makeCluster(nb_coeur, type = "SOCK")
    vglobal<-list(MTEX)
    clusterExport(cl, "vglobal", envir = environment())
    res<-parApply(cl, MARGIN = 1,X = f, FUN = couv99_mtex_seul,pop=MTEX, frac)
    stopCluster(cl)
    T2<-Sys.time() 
    difftime(T2,T1)
    
    h=data.frame(t(res))
    rm(f,res)
    h$ID_Transcrit=(as.character(h$ID_Transcrit))
    h$Position=as.numeric(as.character(h$Position))
    h$PU=as.numeric(as.character(h$PU))
    h$str=(as.character(h$str))
    h$freq=(as.character(h$freq))
    h$R2_99pour=as.numeric(as.character(h$R2_99pour))
    h$lg_R2_99=as.numeric(as.character(h$lg_R2_99))
    h$lg_agreg=as.numeric(as.character(h$lg_agreg))
    h$sumfreq=as.numeric(as.character(h$sumfreq))
    h$nb_int=as.numeric(as.character(h$nb_int))
    h$nb_chev=as.numeric(as.character(h$nb_chev))
    h$fin_couple_chev=as.numeric(as.character(h$fin_couple_chev))
    h$nb_fin_chev=as.numeric(as.character(h$nb_fin_chev))
    
    
    h$rap=h$nb_chev/h$nb_int
    seuil_ratio_chev_int=ifelse(is.na(Fmin),summary(h$rap)[2],Fmin)
    
    print(paste("Fmin =",seuil_ratio_chev_int))
    
    vglobal=list(MTEX)
    list1=recursive3(h,h,pop=MTEX,frac,nb_coeur,seuil_ratio=seuil_ratio_chev_int,dem_TEX="NON")
    list1[[2]]$nb_allong=0
    
    fin=list1[[2]]
    compt=1
    while(is.null(list1[[1]])==FALSE){
      list1=recursive3(list1[[1]],h,MTEX,frac,nb_coeur,seuil_ratio=seuil_ratio_chev_int,dem_TEX="NON")
      
      if(nrow(list1[[2]])>0){
        list1[[2]]$nb_allong=compt
        fin=rbind(fin,list1[[2]])
      }
      compt=compt+1 
    }
    
    
    fi=fin[,-(6:7)]
    fi$lg=as.numeric(as.character(fi$lg_agreg))
    fi$Position=as.numeric(as.character(fi$Position))
    fi$PU=as.numeric(as.character(fi$PU))
    fi$ID_Transcrit=as.character(fi$ID_Transcrit)
    fi$freq=(as.character(fi$freq))
    
    fi=fi[order(fi$Position,decreasing = F),]
  }
  
    #Annotation ----
   if (!is.null(ptt_file)) {
    f = fi[, c(2, 1, 4, 5, 14, 3)]
    colnames(f)[6] = "Positional.Uncertainty"
    
    save_position=f$Position
    save_pu=f$Positional.Uncertainty
    id_order=f$ID_Transcrit
    
    f$Position = ifelse(f$str == "+", ((f$Position - f$Positional.Uncertainty) + 
                                         (f$Position + f$lg - 1))/2, ((f$Position + f$Positional.Uncertainty) + 
                                                                        (f$Position - f$lg + 1))/2)
    f$Positional.Uncertainty = ifelse(f$str == "+", ((f$Position + 
                                                        f$lg - 1) - (f$Position - f$Positional.Uncertainty))/2, 
                                      ((f$Position - f$lg + 1) - (f$Position + f$Positional.Uncertainty))/2)
    fii = annot_apply_RNA_ptt(f, ptt, genome_size)


    fin = fi[, c(1:3, 14, 4, 5, 13, 12)]
    colnames(fin)[3] = "Positional.Uncertainty"
    colnames(fin)[7] = "iteration_nb"
    colnames(fin)[8] = "last_Ftsse"
    
    fin$Class=fii$Class[match(fin$ID_Transcrit,fii$ID_Transcrit)]
    fin$Comment=fii$Comment[match(fin$ID_Transcrit,fii$ID_Transcrit)] 
  }
  else {
    fin = fi[, c(1:3, 14, 4, 5, 13, 12)]
    colnames(fin)[3] = "Positional.Uncertainty"
    colnames(fin)[7] = "iteration_nb"
    colnames(fin)[8] = "last_Ftsse"
  }
  #----
  print("DONE")
  return(fin)
}
