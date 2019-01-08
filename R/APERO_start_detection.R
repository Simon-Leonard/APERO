APERO_start_detection <-
function(work_dir=getwd(),bam_name,ptt_file=NULL,wmax,min_dist,enrichment,min_read_number,genome_size){
  options(digits=20)
  setwd(work_dir)
  
  if (wmax<0){return("Please use a positive value for wmax")}
  if(min_dist<1){return("d has to be higher than 1")}
  if(genome_size<0){return("Genome size has to be a positive value")}

  
  ### Functions ----
  library("Rsamtools")
  library(reshape2)
  qw=function(b){
    let=as.numeric(gregexpr("[A-Z]",as.character(b[6]))[[1]])
    let=c(0,let,nchar(as.character(as.character(b[6]))))
    le=as.numeric(gregexpr("[M|D|N|X|=]",as.character(b[6]))[[1]])
    sum=0
    for(i in 1:length(le)){sum=sum+as.numeric(substring(as.character(b[6]),let[max(which(let<le[i]))]+1,le[i]-1))}
    return(sum)
  }
  NET5=function(test,window_size,flanks_size,genome_size){# <---flank---[window]---flank--->
    
    colnames(test)[1]="5end"
    test$`5end`=as.numeric(as.character(test$`5end`))
    test$`sum_freq`=as.numeric(as.character(test$`sum_freq`))
    dp=test[which(test$str=="+"),]
    dn=test[which(test$str=="-"),]
    
    M=data.frame(P=rep(0,genome_size),N=rep(0,genome_size))
    M$P[dp$'5end']=dp$sum_freq
    M$N[dn$'5end']=dn$sum_freq
    
    if (window_size>1){
      www=(window_size-1)/2
      a=c()
      b=c()
      mp=which(M$P==0)
      mn=which(M$N==0)
      for (x in 1:www){
        spp=mp[which((mp+x)%in%dp$'5end'[which(dp$sum_freq!=0)])]
        spn=mp[which((mp-x)%in%dp$'5end'[which(dp$sum_freq!=0)])]
        if (sum(spn%in%spp)>=1){
          sp=c(spp,spn[-which(spn%in%spp)])
        }else{sp=c(spp,spn)}
        if (sum(sp%in%a)>=1){
          sp=sp[-which(sp%in%a)]
        }
        a=c(a,sp)
        if (sum(a%in%dp$'5end')>=1){
          if (sum(a%in%dp$'5end')==length(a)) {
            a=c()
          }else{ a=a[-which(a%in%dp$'5end')]}
        }
        
        snp=mn[which((mn+x)%in%dn$'5end'[which(dn$sum_freq!=0)])]
        snn=mn[which((mn-x)%in%dn$'5end'[which(dn$sum_freq!=0)])]
        if (sum(snn%in%snp)>=1){
          sn=c(snp,snn[-which(snn%in%snp)])
        }else{sn=c(snp,snn)}
        if (sum(sn%in%b)>=1){
          sn=sn[-which(sn%in%b)]
        }
        b=c(b,sn)
        if (sum(b%in%dn$'5end')>=1){
          if (sum(b%in%dn$'5end')==length(b)){
            b=c()
          }else{b=b[-which(b%in%dn$'5end')]}
        }
      }
      
      extension=0
      if(length(a)>0){
        mp=data.frame('5end'=a,sum_freq=0,str="+")
        extension=1
      }else{mp=NULL}
      if(length(b)>0){
        mn=data.frame('5end'=b,sum_freq=0,str="-")
        extension=1
      }else{mn=NULL}
      
      if(extension==1){
        m=rbind(mn,mp)
        c=ncol(test)-ncol(m)
        m=data.frame(m,matrix(NA,nrow = nrow(m),ncol=c))
        colnames(m)=colnames(test)
        test=rbind(test,m)
        test=test[order(test$"5end",decreasing=FALSE),]
      }  
    }
    
    test$ws=(window_size-1)/2
    test$fs=flanks_size
    i=which(colnames(test)=="ws")
    j=which(colnames(test)=="fs")
    NET_unique=function(vec) {
      ws=as.numeric(as.character(vec[i]))
      fs=as.numeric(as.character(vec[j]))
      a=(as.numeric(vec[1])-ws):(as.numeric(as.character(vec[1]))+ws)
      b=(as.numeric(as.character(vec[1]))-ws-fs):(as.numeric(as.character(vec[1]))-ws-1)
      c=(as.numeric(as.character(vec[1]))+ws+fs):(as.numeric(as.character(vec[1]))+ws+1)
      if (vec[3]=="+"){
        nb_w=mean(M$P[a[a>0]])
        nb_f1=mean(M$P[b[b>0]])
        nb_f2=mean(M$P[c[c>0]])
        nb_f=mean(c(nb_f1,nb_f2))
        net=nb_f/nb_w
      }else {
        nb_w=mean(M$N[a[a>0]])
        nb_f1=mean(M$N[b[b>0]])
        nb_f2=mean(M$N[c[c>0]])
        nb_f=mean(c(nb_f1,nb_f2))
        net=nb_f/nb_w
      }
      
      return(net)
    }
    
    t=apply(test,1,NET_unique)
    test$ws=NULL
    test$fs=NULL
    test$net=t
    colnames(test)[which(colnames(test)=="net")]=paste("net_w",window_size,"_f",flanks_size,sep="",collapse=NULL)
    return(test)
  } # Start enrichment measurment
  agreg4<-function(tot2,ff){
    colnames(tot2)[1]="5end"
    colnames(ff)[1]="5end"
    tot2=tot2[order(tot2$'str',tot2$'5end',decreasing=FALSE),]
    p=0
    s=1
    while (s!=p){
      i=1
      s=nrow(tot2)
      while (i <=(nrow(tot2)-1)){
        if ((tot2$'5end'[i]+tot2$PU[i])>=(tot2$'5end'[i+1]-tot2$PU[i+1])& tot2$str[i]==tot2$str[i+1]){
          a=max(tot2$'5end'[i+1]+tot2$PU[i+1],tot2$'5end'[i]+tot2$PU[i])-min(tot2$'5end'[i]-tot2$PU[i],tot2$'5end'[i+1]-tot2$PU[i+1])
          b=mean(c(max(tot2$'5end'[i+1]+tot2$PU[i+1],tot2$'5end'[i]+tot2$PU[i]),min(tot2$'5end'[i]-tot2$PU[i],tot2$'5end'[i+1]-tot2$PU[i+1])))
          tot2$PU[i]=a/2
          tot2$'5end'[i]=b
          tot2=tot2[-(i+1),]
        }
        i=i+1
      }
      p=nrow(tot2)
    }
    
    sumfreq=function(tot2,ff){
      END=as.numeric(as.character(tot2[1]))
      PU=as.numeric(as.character(tot2[5]))
      brin=as.character(tot2[3])
      tot2[2]=(sum(ff$sum_freq[which(ff$`5end`%in%((END-PU):(END+PU)) & ff$str==brin)]))
      return(tot2)
    }
    
    tot2=(data.frame(t(apply(tot2,1,sumfreq,ff=ff))))
    tot2$`X5end`=as.numeric(as.character(tot2$`X5end`))
    tot2$sum_freq=as.numeric(as.character(tot2$sum_freq))
    tot2$PU=as.numeric(as.character(tot2$PU))
    
    return(tot2)
  }
  benizalon2<-function(data,genome_size,str=c("+","-")) { # Make a coverage matrix of the detected TSS
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
  
  detection_demarrages=function(d,wmax=7,flanks=c(3,5),seuil_net=0.1,genome_size){
    wmin=1 # Initial value
    int=seq(wmin,wmax,2)
    res=data.frame(window=int,count=0)
    
    for (i in int){
      
      cat(paste("\n","Start detection with window = ",i,sep="",collapse=NULL))
      
      for (j in flanks){
        dd=NET5(d,i,j,genome_size)
        if(j==flanks[1]){ind=(1:nrow(dd))}
        f=which(dd[,4]<=seuil_net)
        ind=ind[ind %in% f]
      }
      
      net=dd[ind,]
      net$PU=(i-1)/2
      
      if(i == int[1]){
        # window = 1, no aggregation
        fin=net[,c(1:3,5)]
        res$count[res$window==i]=nrow(fin)
      }else {
        cat("\n \t Deletion of overlapping detections")
        Mp=benizalon2(fin,genome_size,str="+")
        Mn=benizalon2(fin,genome_size,str="-")
        
        
        fungarde=function(n,M){
          ri=as.numeric(as.character(n[1]))+as.numeric(as.character(n[5]))
          le=as.numeric(as.character(n[1]))-as.numeric(as.character(n[5]))
          int=le:ri
          return(ifelse(sum(M[int])>0,"NON","OUI"))
        }
        
        np=net[net$str=="+",]
        if(nrow(np)>0){np$garde=apply(X = np,MARGIN = 1,FUN = fungarde, M=Mp)}
        
        nn=net[net$str=="-",]
        if(nrow(nn)>0){nn$garde=apply(X = nn,MARGIN = 1,FUN = fungarde, M=Mn)}
        
        n=rbind(np[np$garde=="OUI",1:5],nn[nn$garde=="OUI",1:5])
        
        if (nrow(n)>0){
          cat("\n \t Agregation of new starts")
          net3=agreg4(n,d)
          colnames(net3)[1]="5end"
          
          fin=rbind(fin,net3[,c(1:3,5)])
          res$count[res$window==i]=nrow(net3)
        } else {
          cat("\n \t No more starts detected with this window")
          res$count[res$window==i]=0}
      }
      
      
    }
    cat ("\n")
    print(res)
    fin=fin[,c(1,4,3,2)]
    colnames(fin)=c("Position","Positional.Uncertainty","str","freq")
    fin=fin[order(fin$Position,decreasing = F),]
    return(fin)
  }
  annot_apply_ARN_ptt=function(data,ref,genome_size){
    options(digits=20)
    
    class<-function(data){ #Final annotation function
      data$Class=as.character(data$Class)
      nb=data.frame(Class=data$Class,P=0,I=0,troisU=0,cinqU=0,Chev=0,Ad=0,Ai=0,Av=0,AA=0,O=0,div=0)
      cla<-function(nb){
        if(gregexpr("P", nb[1])[[1]][1]>=0){nb[2]=1}
        if(gregexpr("I", nb[1])[[1]][1]>=0){nb[3]=1}
        if(gregexpr("3U", nb[1])[[1]][1]>=0){nb[4]=1}
        if(gregexpr("5U", nb[1])[[1]][1]>=0){nb[5]=1}
        if(gregexpr("Chev", nb[1])[[1]][1]>=0){nb[6]=1}
        if(gregexpr("Ad", nb[1])[[1]][1]>=0){nb[7]=1}
        if(gregexpr("Ai", nb[1])[[1]][1]>=0){nb[8]=1}
        if(gregexpr("Av", nb[1])[[1]][1]>=0){nb[9]=1}
        if(gregexpr("AA", nb[1])[[1]][1]>=0){nb[10]=1}
        if(gregexpr("O", nb[1])[[1]][1]>=0){nb[11]=1}
        if(gregexpr("Div", nb[1])[[1]][1]>=0){nb[12]=1}
        return(nb)
      }
      nb=data.frame(t(apply(nb,1,cla)))
      nb$P=as.numeric(as.character(nb$P))
      nb$I=as.numeric(as.character(nb$I))
      nb$troisU=as.numeric(as.character(nb$troisU))
      nb$Ad=as.numeric(as.character(nb$Ad))
      nb$Ai=as.numeric(as.character(nb$Ai))
      nb$AA=as.numeric(as.character(nb$AA))
      nb$O=as.numeric(as.character(nb$O))
      nb$div=as.numeric(as.character(nb$div))
      nb$Av=as.numeric(as.character(nb$Av))
      nb$Chev=as.numeric(as.character(nb$Chev))
      nb$cinqU=as.numeric(as.character(nb$cinqU))
      
      nb$Class=0
      
      then<-function(nb){
        if(nb[2]!=0){nb[1]="P"}
        if(nb[3]!=0){
          if(nb[1]==0){nb[1]="I"}else{nb[1]=paste(nb[1],"I",sep=";")}
        }
        if(nb[4]!=0){
          if(nb[1]==0){nb[1]="3U"}else{nb[1]=paste(nb[1],"3U",sep=";")}
        }
        if(nb[5]!=0){
          if(nb[1]==0){nb[1]="5U"}else{nb[1]=paste(nb[1],"5U",sep=";")}
        }
        if(nb[6]!=0){
          if(nb[1]==0){nb[1]="Chev"}else{nb[1]=paste(nb[1],"Chev",sep=";")}
        }
        if(nb[7]!=0){
          if(nb[1]==0){nb[1]="3A"}else{nb[1]=paste(nb[1],"3A",sep=";")}
        }
        if(nb[8]!=0){
          if(nb[1]==0){nb[1]="Ai"}else{nb[1]=paste(nb[1],"Ai",sep=";")}
        }
        if(nb[9]!=0){
          if(nb[1]==0){nb[1]="5A"}else{nb[1]=paste(nb[1],"5A",sep=";")}
        }
        if(nb[10]!=0){
          if(nb[1]==0){nb[1]="Achev"}else{nb[1]=paste(nb[1],"Achev",sep=";")}
        }
        if(nb[12]!=0){
          if(nb[1]==0){nb[1]="Div"}else{nb[1]=paste(nb[1],"Div",sep=";")}
        }
        
        if(nb[11]!=0){nb[1]="O"}
        return(nb)
      }
      nb=data.frame(t(apply(nb,1,then)))  
      data$Class=nb$Class
      return(data)
    }
    
    data$str=as.character(data$str)
    data$Position=as.numeric(as.character(data$Position))
    data$Positional.Uncertainty=as.numeric(as.character(data$Positional.Uncertainty))
    
    ref=data.frame(colsplit(ref$Location,"\\..",c("Left","Right")),ref[,-1])
    ref=ref[,c(5,1,2,5,6,3,8,7,10)]
    
    colnames(ref)=c("Feature.Type","Left.End.ASAP","Right.End.ASAP","product","Gene.Symbol.ASAP","Strand","ID.ASAP","locus.tag","description.ASAP")
    ref$Gene.Symbol.ASAP=as.character(ref$Gene.Symbol.ASAP)
    
    
    data$Class=0
    data$Comment=0
    
    refn=ref[which(ref$Strand=="-"),] #Annontation on minus strand
    refp=ref[which(ref$Strand=="+"),]#Annontation on plus strand
    
    datap=data[which(data$str=="+"),]
    datan=data[which(data$str=="-"),]
    
    refn=refn[order(refn$Left.End.ASAP,decreasing=F),]# Ordering
    refp=refp[order(refp$Left.End.ASAP,decreasing=F),]
    
    datap=datap[order(datap$Position),]
    datan=datan[order(datan$Position),]
    refn=rbind(c("CDD",0,0,"Start_Genome","Start_Genome","Start_Genome","Start_Genome","Start_Genome","Start_Genome"),refn)
    refp=rbind(c("CDD",0,0,"Start_Genome","Start_Genome","Start_Genome","Start_Genome","Start_Genome","Start_Genome"),refp)
    
    refn=rbind(refn,c("CDD",genome_size,genome_size,"End_Genome","End_Genome","End_Genome","End_Genome","End_Genome","End_Genome"))
    refp=rbind(refp,c("CDD",genome_size,genome_size,"End_Genome","End_Genome","End_Genome","End_Genome","End_Genome","End_Genome"))
    
    #Annotation minus strand
    if(nrow(datan)>=1){
      annot_moins=function(datan){
        datan[1]=as.numeric(as.character(datan[1]))
        datan[6]=as.numeric(as.character(datan[6]))
        p=0# no match
        if (length(which(as.numeric(refp$Right.End.ASAP)<(as.numeric(datan[1])-as.numeric(datan[6])-10000)))>=1){
          j=max(which(as.numeric(refp$Right.End.ASAP)<(as.numeric(datan[1])-as.numeric(datan[6])-10000)))
        }else {j=1}
        if (length(which(as.numeric(refn$Right.End.ASAP)<(as.numeric(datan[1])-as.numeric(datan[6])-10000)))>=1) {
          k=max(which(as.numeric(refn$Right.End.ASAP)<(as.numeric(datan[1])-as.numeric(datan[6])-10000)))
        }else {k=1}
        arret_boucle=0
        while((as.numeric(refp[j,2])-(as.numeric(datan[1])+as.numeric(datan[6])))<=500 & arret_boucle==0){
          if (((as.numeric(datan[1])-as.numeric(datan[6])))>as.numeric(refp[j,2]) &((as.numeric(datan[1])-as.numeric(datan[6]))<as.numeric(refp[j,3])) & ((as.numeric(datan[1])+as.numeric(datan[6]))>as.numeric(refp[j,3]))){
            if(datan[6]==0){a=paste("antisense to gene ",refp[j,5]," (",(as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refp[j,3])),"nt downstream)")} else{a=paste("antisense to gene ",refp[j,5]," (between",(as.numeric(datan[1])+as.numeric(datan[6]))-as.numeric(refp[j,3]),"and",(as.numeric(datan[1])-as.numeric(datan[6]))-as.numeric(refp[j,3]),"nt downstream)")}
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else {datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Ad",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="Ad"}else {datan[7]=paste(datan[7],"Ad",sep=";")}
          }else if (((as.numeric(datan[1])+as.numeric(datan[6]))<=as.numeric(refp[j,3])) & ((as.numeric(datan[1])-as.numeric(datan[6]))>=as.numeric(refp[j,2])) ){
            a=paste("antisense to gene ",refp[j,5])
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Ai",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="Ai"}else {datan[7]=paste(datan[7],"Ai",sep=";")}
          }else if (((as.numeric(datan[1])+as.numeric(datan[6]))<=as.numeric(refp[j,2])) & (as.numeric(refp[j,2])-as.numeric(datan[1])-as.numeric(datan[6]))<=300 ){
            if (datan[6]==0){a=paste((as.numeric(refp[j,2])-as.numeric(datan[1])),"nt diverging with gene ",refp[j,5])}else {a=paste("Between ",(as.numeric(refp[j,2])-as.numeric(datan[1])-as.numeric(datan[6])),"nt and ",(as.numeric(refp[j,2])-as.numeric(datan[1])+as.numeric(datan[6])),"nt diverging with gene ",refp[j,5])}
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Div",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="Div"}else {datan[7]=paste(datan[7],"Div",sep=";")}
          }else if (((as.numeric(datan[1])+as.numeric(datan[6])))>as.numeric(refp[j,2]) &((as.numeric(datan[1])+as.numeric(datan[6]))<as.numeric(refp[j,3])) & ((as.numeric(datan[1])-as.numeric(datan[6]))<as.numeric(refp[j,2]))){
            if(datan[6]==0){a=paste("antisense to gene ",refp[j,5]," (",(as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refp[j,3])),"nt downstream)")} else{a=paste("antisense to gene ",refp[j,5]," (between",(as.numeric(datan[1])+as.numeric(datan[6]))-as.numeric(refp[j,3]),"and",(as.numeric(datan[1])-as.numeric(datan[6]))-as.numeric(refp[j,3]),"nt downstream)")}
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Av",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="Av"}else {datan[7]=paste(datan[7],"Av",sep=";")} 
          }else if (as.numeric(refp[j,2])>=(as.numeric(datan[1])-as.numeric(datan[6])) & (as.numeric(refp[j,3])<=(as.numeric(datan[1])+as.numeric(datan[6])))){
            if(datan[6]==0){a=paste("Ai overlaping gene ",refp[j,5])} else{a=paste("Ai overlaping gene ",refp[j,5])}
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else {datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("AA",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="AA"}else {datan[7]=paste(datan[7],"AA",sep=";")}
          }
          j=j+1
          if(j>nrow(refp)){
            arret_boucle=1
            j=j-1
          }
        }
        arret_boucle=0
        while((as.numeric(datan[1])+as.numeric(datan[6]))>=as.numeric(refn[k,2])-300 & arret_boucle==0){
          if ((as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refn[k,3]))<=250 & (as.numeric(datan[1])-as.numeric(datan[6])>=as.numeric(refn[k,3])) ){
            if (datan[6]==0){a=paste((as.numeric(datan[1])+as.numeric(datan[6]))-as.numeric(refn[k,3]),"nt upstream",refn[k,5])}else{a=paste("Between",(as.numeric(datan[1])-as.numeric(datan[6])-as.numeric(refn[k,3])),"and",(as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refn[k,3])),"nt upstream ",refn[k,5])}
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("P",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="P"}else {datan[7]=paste(datan[7],"P",sep=";")}
          }else if (((as.numeric(refn[k,3]))<(as.numeric(datan[1])+as.numeric(datan[6]))) & ((as.numeric(refn[k,3]))>(as.numeric(datan[1])-as.numeric(datan[6]))) & ((as.numeric(refn[k,2]))<(as.numeric(datan[1]) -as.numeric(datan[6]))) ){
            a=paste("Between",(as.numeric(datan[1])-as.numeric(datan[6])-as.numeric(refn[k,3])),"and",(as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refn[k,3])),"nt upstream ",refn[k,5])
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("5U",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="5U"}else {datan[7]=paste(datan[7],"5U",sep=";")}
          }else if (((as.numeric(refn[k,3]))-(as.numeric(datan[1])+as.numeric(datan[6])))>=0 &((as.numeric(refn[k,2]))-(as.numeric(datan[1])-as.numeric(datan[6])))<=0 ){
            a=paste("Whithin ",refn[k,5],", ",(as.numeric(refn[k,3])-(as.numeric(datan[1])+as.numeric(datan[6]))),"nt after start and ",((as.numeric(refn[k,2]))-(as.numeric(datan[1])-as.numeric(datan[6]))),"nt before end")
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("I",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="I"}else {datan[7]=paste(datan[7],"I",sep=";")}
          } else if (((as.numeric(refn[k,2]))<(as.numeric(datan[1])+as.numeric(datan[6]))) & ((as.numeric(refn[k,2]))>(as.numeric(datan[1])-as.numeric(datan[6]))) & ((as.numeric(refn[k,3]))>(as.numeric(datan[1]) + as.numeric(datan[6]))) ){
            a=paste("3UTR of ",refn[k,5],", ",((as.numeric(refn[k,2]))-(as.numeric(datan[1])+as.numeric(datan[6]))),"nt before the end and",((as.numeric(refn[k,2]))-(as.numeric(datan[1])-as.numeric(datan[6]))),"after")
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("3U",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="3U"}else {datan[7]=paste(datan[7],"3U",sep=";")}
          }else if ((as.numeric(refn[k,3])<=(as.numeric(datan[1])+as.numeric(datan[6]))) & (as.numeric(refn[k,2])>=(as.numeric(datan[1])-as.numeric(datan[6]))) ){
            a=paste("Overlaping gene ",refn[k,5])
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Chev",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="Chev"}else {datan[7]=paste(datan[7],"Chev",sep=";")}
          }
          k=k+1
          if(k>nrow(refn)){
            arret_boucle=1
            k=k-1
          }
        }  
        if (as.character(datan[7])=="0"){
          datan[7]="O"
          a1=as.numeric(as.character(datan[1]))-as.numeric(as.character(datan[6]))
          a2=as.numeric(as.character(datan[1]))+as.numeric(as.character(datan[6]))
          a=max(which(as.numeric(refn$Right.End.ASAP)<(a1)))
          b=min(which(as.numeric(refn$Left.End.ASAP)>(a2)))
          c1=a1-as.numeric(refn[a,3])
          c2=a2-as.numeric(refn[a,3])
          d1=as.numeric(refn[b,2])-(a2)
          d2=as.numeric(refn[b,2])-(a1)
          if(datan[6]!=0){e=paste("Between ",d1,"nt & ",d2,"nt after",refn[b,5],"and between",c1,"nt & ",c2,"nt before",refn[a,5])
          }else e=paste(d1,"nt after",refn[b,5],"and ",c1,"nt before",refn[a,5])
          datan[8]=e
        }
        if(j<=1){j=j+1}
        if(k<=1){k=k+1}
        return(datan)
      }
      datan=data.frame(t(apply(datan,1,annot_moins)))
      
    }
    
    # Annotation plus strand
    if(nrow(datap)>=1){
      annot_plus=function(datap){
        datap[1]=as.numeric(as.character(datap[1]))
        datap[6]=as.numeric(as.character(datap[6]))
        p=0#no match
        if (length(which(as.numeric(refn$Right.End.ASAP)<(as.numeric(datap[1])-as.numeric(datap[6])-10000)))>=1){
          j=max(which(as.numeric(refn$Right.End.ASAP)<(as.numeric(datap[1])-as.numeric(datap[6])-10000)))
        }else {j=1}
        if (length(which(as.numeric(refp$Right.End.ASAP)<(as.numeric(datap[1])-as.numeric(datap[6])-10000)))>=1){
          k=max(which(as.numeric(refp$Right.End.ASAP)<(as.numeric(datap[1])-as.numeric(datap[6])-10000)))
        }else {k=1}
        arret_boucle=0
        while((as.numeric(datap[1])+as.numeric(datap[6]))>=(as.numeric(refn[j,2])) & arret_boucle==0){
          if (as.numeric(refn[j,2])>(as.numeric(datap[1])-as.numeric(datap[6])) & (as.numeric(refn[j,2])<(as.numeric(datap[1])+as.numeric(datap[6]))) & as.numeric(refn[j,3])>(as.numeric(datap[1])+as.numeric(datap[6]))){
            if(datap[6]==0){a=paste("antisense to gene ",refn[j,5]," (",as.numeric(refn[j,2])-(as.numeric(datap[1])-as.numeric(datap[6])),"nt downstream)")} else{a=paste("antisense to gene ",refn[j,5]," (between",as.numeric(refn[j,2])-(as.numeric(datap[1])+as.numeric(datap[6])),"and",as.numeric(refn[j,2])-(as.numeric(datap[1])-as.numeric(datap[6])),"nt downstream)")}
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else {datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Ad",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="Ad"}else {datap[7]=paste(datap[7],"Ad",sep=";")}
          }else if ((as.numeric(refn[j,3])>=(as.numeric(datap[1])+as.numeric(datap[6]))) & (as.numeric(refn[j,2])<=(as.numeric(datap[1])-as.numeric(datap[6])))){
            a=paste("antisense to gene ",refn[j,5])
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Ai",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="Ai"}else {datap[7]=paste(datap[7],"Ai",sep=";")}
          }else if (((as.numeric(datap[1])-as.numeric(datap[6]))>=as.numeric(refn[j,3])) &((as.numeric(datap[1])-as.numeric(datap[6]))-as.numeric(refn[j,3]))<=300){
            if (datap[6]==0){a=paste((as.numeric(datap[1])-as.numeric(refn[j,3])),"nt diverging with gene ",refn[j,5])}else {a=paste("Between ",(as.numeric(datap[1])-as.numeric(refn[j,3])-as.numeric(datap[6])),"nt and ",(as.numeric(datap[1])-as.numeric(refn[j,3])+as.numeric(datap[6])),"nt diverging with gene ",refn[j,5])}
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Div",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="Div"}else {datap[7]=paste(datap[7],"Div",sep=";")}
          }else if (as.numeric(refn[j,3])>(as.numeric(datap[1])-as.numeric(datap[6])) & (as.numeric(refn[j,3])<(as.numeric(datap[1])+as.numeric(datap[6]))) & (as.numeric(refn[j,2])<(as.numeric(datap[1])-as.numeric(datap[6])))){
            if(datap[6]==0){a=paste("antisense to gene ",refn[j,5]," (",as.numeric(refn[j,3])-(as.numeric(datap[1])-as.numeric(datap[6])),"nt downstream)")} else{a=paste("antisense to gene ",refn[j,5]," (between",as.numeric(refn[j,3])-(as.numeric(datap[1])+as.numeric(datap[6])),"and",as.numeric(refn[j,3])-(as.numeric(datap[1])-as.numeric(datap[6])),"nt downstream)")}
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else {datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Av",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="Av"}else {datap[7]=paste(datap[7],"Av",sep=";")}
          }else if (as.numeric(refn[j,2])>=(as.numeric(datap[1])-as.numeric(datap[6])) & (as.numeric(refn[j,3])<=(as.numeric(datap[1])+as.numeric(datap[6])))){
            if(datap[6]==0){a=paste("Ai overlaping gene ",refn[j,5])} else{a=paste("Ai overlaping gene ",refn[j,5])}
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else {datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("AA",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="AA"}else {datap[7]=paste(datap[7],"AA",sep=";")}
          }
          j=j+1
          if(j>nrow(refn)){
            arret_boucle=1
            j=j-1
          }
        }
        arret_boucle=0
        while((as.numeric(datap[1])-as.numeric(datap[6]))>=((as.numeric(refp[k,2]))-10000) & arret_boucle==0){
          if ((as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6])))<=250 & (as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6])))>=0 & (as.numeric(refp[k,2])>=(as.numeric(datap[1])+as.numeric(datap[6]))) ){
            if (datap[6]==0){a=paste((as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"nt upstream",refp[k,5])}else{a=paste("Between",(as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"and",(as.numeric(refp[k,2])-(as.numeric(datap[1])+as.numeric(datap[6]))),"nt upstream ",refp[k,5])}
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("P",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="P"}else {datap[7]=paste(datap[7],"P",sep=";")}
          }else if ((as.numeric(refp[k,2])>(as.numeric(datap[1])-as.numeric(datap[6]))) & (as.numeric(refp[k,2])<(as.numeric(datap[1])+as.numeric(datap[6]))) & (as.numeric(refp[k,3])>(as.numeric(datap[1])+as.numeric(datap[6]))) ){
            if (datap[6]==0){a=paste((as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"nt upstream",refp[k,5])}else{a=paste("Between",(as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"and",(as.numeric(refp[k,2])-(as.numeric(datap[1])+as.numeric(datap[6]))),"nt upstream ",refp[k,5])}
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("5U",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="5U"}else {datap[7]=paste(datap[7],"5U",sep=";")}
          }else if ((as.numeric(refp[k,2])<=(as.numeric(datap[1])-as.numeric(datap[6]))) & (as.numeric(refp[k,3])>=(as.numeric(datap[1])+as.numeric(datap[6])))){
            a=paste("Whithin ",refp[k,5],", ",((as.numeric(refp[k,2]))-(as.numeric(datap[1])-as.numeric(datap[6]))),"nt after start and ",(as.numeric(refp[k,3])-(as.numeric(datap[1])+as.numeric(datap[6]))),"nt before end")
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("I",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="I"}else {datap[7]=paste(datap[7],"I",sep=";")}
          }else if ((as.numeric(refp[k,3])<(as.numeric(datap[1])+as.numeric(datap[6]))) & (as.numeric(refp[k,3])>(as.numeric(datap[1])-as.numeric(datap[6]))) & (as.numeric(refp[k,2])<(as.numeric(datap[1])-as.numeric(datap[6]))) ){
            a=paste("3UTR of ",refp[k,5],", ",((as.numeric(refp[k,3]))-(as.numeric(datap[1])+as.numeric(datap[6]))),"nt before the end and ",((as.numeric(refp[k,3]))-(as.numeric(datap[1])-as.numeric(datap[6]))),"after")
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("3U",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="3U"}else {datap[7]=paste(datap[7],"3U",sep=";")}
          }else if ((as.numeric(refp[k,3])<=(as.numeric(datap[1])+as.numeric(datap[6]))) & (as.numeric(refp[k,2])>=(as.numeric(datap[1])-as.numeric(datap[6]))) ){
            a=paste("Overlaping gene ",refp[k,5])
            if (datap[8]==a) {}else if(datap[8]==0){datap[8]=a}else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Chev",datap[7])[[1]][1]
            if(cla>=0){}else if (datap[7]==0){datap[7]="Chev"}else {datap[7]=paste(datap[7],"Chev",sep=";")}
          }
          k=k+1
          if(k>nrow(refp)){
            arret_boucle=1
            k=k-1
          }
        }  
        if (as.character(datap[7])=="0"){
          datap[7]="O"
          a1=as.numeric(as.character(datap[1]))-as.numeric(as.character(datap[6]))
          a2=as.numeric(as.character(datap[1]))+as.numeric(as.character(datap[6]))
          a=max(which(as.numeric(refp$Right.End.ASAP)<(a1)))
          b=min(which(as.numeric(refp$Left.End.ASAP)>(a2)))
          c1=a1-as.numeric(refp[a,3])
          c2=a2-as.numeric(refp[a,3])
          d1=as.numeric(refp[b,2])-(a2)
          d2=as.numeric(refp[b,2])-(a1)
          if(datap[6]!=0){e=paste("Between ",c1,"nt & ",c2,"nt after",refp[a,5],"and between",d1,"nt & ",d2,"nt before",refp[b,5])
          }else e=paste(c1,"nt after",refp[a,5],"and ",d1,"nt before",refp[b,5])
          datap[8]=e
        }
        if(j<=1){j=j+1}
        if(k<=1){k=k+1}
        
        return(datap) 
      }
      datap=data.frame(t(apply(datap,1,annot_plus)))
    }
    data=rbind(datap,datan)
    data$Position=as.numeric(as.character(data$Position))
    if (length(data$X5end)>0){data$X5end=as.numeric(as.character(data$X5end))}
    if(length(data$mate5end)!=0){data$mate5end=as.numeric(as.character(data$mate5end))}
    # data$freq=as.numeric(as.character(data$freq))
    data$Positional.Uncertainty=as.numeric(as.character(data$Positional.Uncertainty))
    data=data[order(data$Position),]
    data=class(data)
    
    return(data)
  }#End function annot
  
  ### BAM import ----
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
    
    FLAG=c(81,83,97,99,145,147,161,163)
    bam=bam[which(bam$flag%in%FLAG),]
    rm(FLAG)
    
    bam$qwidth=apply(bam,1,qw) ### Apply with snowfall crash
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
    
  
  ### Start detection ----
  flanks=c(1:min_dist)
  dem=detection_demarrages(d,wmax=(2*wmax)+1,flanks=flanks,seuil_net=enrichment,genome_size=genome_size)
  
  dem=dem[which(dem$freq>min_read_number),]
  
  ### annotation ----
  if (!is.null(ptt_file)){
    dem$C1=NA
    dem$C2=NA
    dem=dem[,c(1,3:6,2)]

    dem=annot_apply_ARN_ptt(dem,ptt_file,genome_size)
    
    dem$C1=NULL
    dem$C2=NULL
  }
  
  ###----
  print("DONE")
  return(dem)
}
