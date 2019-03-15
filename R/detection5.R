  detection5=function(d,wmax=7,flanks=c(3,5),seuil_net=0.1,genome_size){ 
    # 5'end detection function
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
        Mp=cov5(fin,genome_size,str="+")
        Mn=cov5(fin,genome_size,str="-")
        
        
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
          net3=agreg5(n,d)
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
