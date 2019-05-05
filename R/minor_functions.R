###
### This file contains minor functions used is the main functions
###

qw=function(b){ 
    # Use read position and cigar to calculate and return read qwidth
    let=as.numeric(gregexpr("[A-Z]",as.character(b[6]))[[1]])
    let=c(0,let,nchar(as.character(as.character(b[6]))))
    le=as.numeric(gregexpr("[M|D|N|X|=]",as.character(b[6]))[[1]])
    sum=0
    for(i in 1:length(le)){
      sum=sum+as.numeric(substring(as.character(b[6]),let[max(which(let<le[i]))]+1,le[i]-1))
    }
    return(sum)
  }


  cov5<-function(data,genome_size,str=c("+","-")) { 
    # Make a coverage matrix of the detected 5'ends
    m=rep(0,genome_size)
    if(sum(data$str==str)>0){
      data=data[data$str==str,]
      for (i in 1:nrow(data)) {
        for (j in (data[i,1]-data[i,4]):(data[i,1]+data[i,4])){
          m[j]=1/((2*data[i,2])+1)
        }
      }
    }
    
    return(m)
  }


  sumfreq=function(vec){
    # vec[5] contains numbers separated by "_"
    # sumfreq(vec) return the sum of this numbers
    if(grepl("_",vec[5])){
      a=as.numeric(gregexpr("_",vec[5])[[1]])
      a=c(0,a,nchar(vec[5])+1)
      b=c()
      for(i in 1:(length(a)-1)){
        b=c(b,substring(vec[5],a[i]+1,a[i+1]-1))
      }
      return(sum(as.numeric(b)))
    }else{return(as.numeric(vec[5]))}
  }


  max_dem_int=function(f,dem){ 
    # Find the biggest start between the actual start et the actual end
    options(digits=20)
    pu=as.numeric(as.character(f[3]))
    pos=as.numeric(as.character(f[2]))
    lg=as.numeric(as.character(f[14]))
    brin=as.character(f[4])
    # fin=as.numeric(as.character(f[18]))
    fin=ifelse(brin=="+",pos+lg-1,pos-lg+1)
    if (brin=="+"){
      de=dem[dem$Position > ceiling(pos) & dem$Position < fin & dem$str==brin & dem$Position!=pos,]
      if (nrow(de)>0){
        return(max(as.numeric(as.character(de$sumfreq))))
      }else{return(0)}
    }else {
      de=dem[dem$Position < floor(pos) & dem$Position > fin & dem$str==brin & dem$Position!=pos,]
      if (nrow(de)>0){
        return(max(as.numeric(as.character(de$sumfreq))))
      }else{return(0)}
    } 
  }


R2_99pour=function(f,pop,tex=c("yes","no"),frac){ 
  # Take a RNA start and return the longuest end
  # frac=X means "do not take the X% of longuest reads to prevent readthrough
  pos=as.numeric(as.character(f[2]))
  pu=as.numeric(as.character(f[3]))
  brin=as.character(f[4])
  if (tex=="yes"){po=pop[pop$`5end` %in% (pos-pu):(pos+pu) & pop$str==brin,]}
  if (tex=="no"){po=pop[pop$`5end` %in% (floor(pos)-5):(ceiling(pos)+5) & pop$str==brin,]}
  
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


couv99_mtex_seul=function(f,pop,frac){
  # Take one RNA defined as a start position and an end position and calculate :
  #   - Number of reads inside the RNA
  #   - Number of reads overlapping the RNA end
  #   - End position if we consider reads overlapping the RNA end
  #   - Number of reads ending à this end position
  
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


recursive3=function(tab,h,pop,frac,nb_coeur,seuil_ratio,dem_TEX=c("yes","no")){
   # Take a dataframe containing RNA start and end positions
   # Return RNA coordinates with extended positions (under condition)
  
  T1<-Sys.time()
  h1=tab[tab$rap>seuil_ratio,] # 3'end extension
  if(nrow(h1)==0){print("Nothing else to extend");return(list(NULL,tab))}
  h2=tab[tab$rap<=seuil_ratio,] # No 3'end extension
  
  if(dem_TEX=="yes"){
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
