NET5=function(test,window_size,flanks_size,genome_size){
    # Calculate enrichment value with following parameters <---flank---[window]---flank--->  
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
  }
