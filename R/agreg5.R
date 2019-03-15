agreg5<-function(tot2,ff){ 
    # Merge overlaping 5'ends
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
          b=mean(c(max(tot2$'5end'[i+1]+tot2$PU[i+1],tot2$'5end'[i]+tot2$PU[i]),
                   min(tot2$'5end'[i]-tot2$PU[i],tot2$'5end'[i+1]-tot2$PU[i+1])))
          tot2$PU[i]=a/2
          tot2$'5end'[i]=b
          tot2=tot2[-(i+1),]
        }
        i=i+1
      }
      p=nrow(tot2)
    }
    
    sumfreq=function(tot2,ff){
      # Calculate read starts of agregated 5'ends 
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
