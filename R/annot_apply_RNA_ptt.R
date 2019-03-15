annot_apply_RNA_ptt=function(data,ref,genome_size){ 
    # Annotation function
    # Input "ref" = ptt table
    # Input "data" = 6 columns table
    # 1st column = "Position"; 6th column = "Positional.Uncertainty"; One column (2nd, 3rd, 4th ou 5th) has to be "str" (for strand)
    # Return a 8 columns table : The input table, a "Class" column (short annotation) and a "Comment" column (detailed annotation) 
    options(digits=20)
    
    class<-function(data){ #Final annotation function to reorganized classes
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
    
    colnames(ref)=c("Feature.Type","Left.End.ASAP","Right.End.ASAP","product",
                    "Gene.Symbol.ASAP","Strand","ID.ASAP","locus.tag","description.ASAP")
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
          if (((as.numeric(datan[1])-as.numeric(datan[6])))>as.numeric(refp[j,2]) &
              ((as.numeric(datan[1])-as.numeric(datan[6]))<as.numeric(refp[j,3])) & 
              ((as.numeric(datan[1])+as.numeric(datan[6]))>as.numeric(refp[j,3]))){
            if(datan[6]==0){
              a=paste("antisense to gene ",refp[j,5]," (",(as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refp[j,3])),"nt downstream)")} 
            else{a=paste("antisense to gene ",refp[j,5]," (between",
                         (as.numeric(datan[1])+as.numeric(datan[6]))-as.numeric(refp[j,3]),"and",
                         (as.numeric(datan[1])-as.numeric(datan[6]))-as.numeric(refp[j,3]),"nt downstream)")}
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else {datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Ad",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="Ad"}else {datan[7]=paste(datan[7],"Ad",sep=";")}
          }else if (((as.numeric(datan[1])+as.numeric(datan[6]))<=as.numeric(refp[j,3])) & 
                    ((as.numeric(datan[1])-as.numeric(datan[6]))>=as.numeric(refp[j,2])) ){
            a=paste("antisense to gene ",refp[j,5])
            if (datan[8]==a) {}else if(datan[8]==0){datan[8]=a}else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Ai",datan[7])[[1]][1]
            if(cla>=0){}else if (datan[7]==0){datan[7]="Ai"}else {datan[7]=paste(datan[7],"Ai",sep=";")}
          }else if (((as.numeric(datan[1])+as.numeric(datan[6]))<=as.numeric(refp[j,2])) & 
                    (as.numeric(refp[j,2])-as.numeric(datan[1])-as.numeric(datan[6]))<=300 ){
            if (datan[6]==0){
              a=paste((as.numeric(refp[j,2])-as.numeric(datan[1])),"nt diverging with gene ",refp[j,5])
            }else {a=paste("Between ",(as.numeric(refp[j,2])-as.numeric(datan[1])-as.numeric(datan[6])),"nt and ",
                           (as.numeric(refp[j,2])-as.numeric(datan[1])+as.numeric(datan[6])),"nt diverging with gene ",refp[j,5])}
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Div",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="Div"
            }else {datan[7]=paste(datan[7],"Div",sep=";")}
          }else if (((as.numeric(datan[1])+as.numeric(datan[6])))>as.numeric(refp[j,2]) &
                    ((as.numeric(datan[1])+as.numeric(datan[6]))<as.numeric(refp[j,3])) & 
                    ((as.numeric(datan[1])-as.numeric(datan[6]))<as.numeric(refp[j,2]))){
            if(datan[6]==0){
              a=paste("antisense to gene ",refp[j,5]," (",(as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refp[j,3])),"nt downstream)")
            }else{
              a=paste("antisense to gene ",refp[j,5]," (between",(as.numeric(datan[1])+as.numeric(datan[6]))-as.numeric(refp[j,3]),
                      "and",(as.numeric(datan[1])-as.numeric(datan[6]))-as.numeric(refp[j,3]),"nt downstream)")
            }
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Av",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="Av"
            }else {datan[7]=paste(datan[7],"Av",sep=";")} 
          }else if (as.numeric(refp[j,2])>=(as.numeric(datan[1])-as.numeric(datan[6])) & 
                    (as.numeric(refp[j,3])<=(as.numeric(datan[1])+as.numeric(datan[6])))){
            if(datan[6]==0){
              a=paste("Ai overlaping gene ",refp[j,5])
            }else{a=paste("Ai overlaping gene ",refp[j,5])}
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else {datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("AA",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="AA"
            }else {datan[7]=paste(datan[7],"AA",sep=";")}
          }
          j=j+1
          if(j>nrow(refp)){
            arret_boucle=1
            j=j-1
          }
        }
        arret_boucle=0
        while((as.numeric(datan[1])+as.numeric(datan[6]))>=as.numeric(refn[k,2])-300 & arret_boucle==0){
          if ((as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refn[k,3]))<=250 & 
              (as.numeric(datan[1])-as.numeric(datan[6])>=as.numeric(refn[k,3])) ){
            if (datan[6]==0){
              a=paste((as.numeric(datan[1])+as.numeric(datan[6]))-as.numeric(refn[k,3]),"nt upstream",refn[k,5])
            }else{a=paste("Between",(as.numeric(datan[1])-as.numeric(datan[6])-as.numeric(refn[k,3])),"and",
                          (as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refn[k,3])),"nt upstream ",refn[k,5])}
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("P",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="P"
            }else {datan[7]=paste(datan[7],"P",sep=";")}
          }else if (((as.numeric(refn[k,3]))<(as.numeric(datan[1])+as.numeric(datan[6]))) & 
                    ((as.numeric(refn[k,3]))>(as.numeric(datan[1])-as.numeric(datan[6]))) & 
                    ((as.numeric(refn[k,2]))<(as.numeric(datan[1]) -as.numeric(datan[6]))) ){
            a=paste("Between",(as.numeric(datan[1])-as.numeric(datan[6])-as.numeric(refn[k,3])),"and",
                    (as.numeric(datan[1])+as.numeric(datan[6])-as.numeric(refn[k,3])),"nt upstream ",refn[k,5])
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("5U",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="5U"
            }else {datan[7]=paste(datan[7],"5U",sep=";")}
          }else if (((as.numeric(refn[k,3]))-(as.numeric(datan[1])+as.numeric(datan[6])))>=0 & 
                    ((as.numeric(refn[k,2]))-(as.numeric(datan[1])-as.numeric(datan[6])))<=0 ){
            a=paste("Whithin ",refn[k,5],", ",(as.numeric(refn[k,3])-(as.numeric(datan[1])+as.numeric(datan[6]))),
                    "nt after start and ",((as.numeric(refn[k,2]))-(as.numeric(datan[1])-as.numeric(datan[6]))),"nt before end")
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("I",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="I"
            }else {datan[7]=paste(datan[7],"I",sep=";")}
          } else if (((as.numeric(refn[k,2]))<(as.numeric(datan[1])+as.numeric(datan[6]))) & 
                     ((as.numeric(refn[k,2]))>(as.numeric(datan[1])-as.numeric(datan[6]))) & 
                     ((as.numeric(refn[k,3]))>(as.numeric(datan[1]) + as.numeric(datan[6]))) ){
            a=paste("3UTR of ",refn[k,5],", ",((as.numeric(refn[k,2]))-(as.numeric(datan[1])+as.numeric(datan[6]))),
                    "nt before the end and",((as.numeric(refn[k,2]))-(as.numeric(datan[1])-as.numeric(datan[6]))),"after")
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("3U",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="3U"
            }else {datan[7]=paste(datan[7],"3U",sep=";")}
          }else if ((as.numeric(refn[k,3])<=(as.numeric(datan[1])+as.numeric(datan[6]))) & 
                    (as.numeric(refn[k,2])>=(as.numeric(datan[1])-as.numeric(datan[6]))) ){
            a=paste("Overlaping gene ",refn[k,5])
            if (datan[8]==a) {
            }else if(datan[8]==0){
              datan[8]=a
            }else{datan[8]=paste(datan[8],a,sep=";")}
            cla=gregexpr("Chev",datan[7])[[1]][1]
            if(cla>=0){
            }else if (datan[7]==0){
              datan[7]="Chev"
            }else {datan[7]=paste(datan[7],"Chev",sep=";")}
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
          if (as.numeric(refn[j,2])>(as.numeric(datap[1])-as.numeric(datap[6])) & 
              (as.numeric(refn[j,2])<(as.numeric(datap[1])+as.numeric(datap[6]))) & 
              as.numeric(refn[j,3])>(as.numeric(datap[1])+as.numeric(datap[6]))){
            if(datap[6]==0){
              a=paste("antisense to gene ",refn[j,5]," (",as.numeric(refn[j,2])-(as.numeric(datap[1])-as.numeric(datap[6])),
                      "nt downstream)")
            } else{a=paste("antisense to gene ",refn[j,5]," (between",
                           as.numeric(refn[j,2])-(as.numeric(datap[1])+as.numeric(datap[6])),"and",
                           as.numeric(refn[j,2])-(as.numeric(datap[1])-as.numeric(datap[6])),"nt downstream)")}
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else {datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Ad",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="Ad"
            }else {datap[7]=paste(datap[7],"Ad",sep=";")}
          }else if ((as.numeric(refn[j,3])>=(as.numeric(datap[1])+as.numeric(datap[6]))) & 
                    (as.numeric(refn[j,2])<=(as.numeric(datap[1])-as.numeric(datap[6])))){
            a=paste("antisense to gene ",refn[j,5])
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Ai",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="Ai"
            }else {datap[7]=paste(datap[7],"Ai",sep=";")}
          }else if (((as.numeric(datap[1])-as.numeric(datap[6]))>=as.numeric(refn[j,3])) &
                    ((as.numeric(datap[1])-as.numeric(datap[6]))-as.numeric(refn[j,3]))<=300){
            if (datap[6]==0){
              a=paste((as.numeric(datap[1])-as.numeric(refn[j,3])),"nt diverging with gene ",refn[j,5])
            }else {a=paste("Between ",(as.numeric(datap[1])-as.numeric(refn[j,3])-as.numeric(datap[6])),
                           "nt and ",(as.numeric(datap[1])-as.numeric(refn[j,3])+as.numeric(datap[6])),
                           "nt diverging with gene ",refn[j,5])}
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Div",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="Div"
            }else {datap[7]=paste(datap[7],"Div",sep=";")}
          }else if (as.numeric(refn[j,3])>(as.numeric(datap[1])-as.numeric(datap[6])) & 
                    (as.numeric(refn[j,3])<(as.numeric(datap[1])+as.numeric(datap[6]))) & 
                    (as.numeric(refn[j,2])<(as.numeric(datap[1])-as.numeric(datap[6])))){
            if(datap[6]==0){
              a=paste("antisense to gene ",refn[j,5]," (",as.numeric(refn[j,3])-(as.numeric(datap[1])-as.numeric(datap[6])),"nt downstream)")
            }else{a=paste("antisense to gene ",refn[j,5]," (between",as.numeric(refn[j,3])-(as.numeric(datap[1])+as.numeric(datap[6])),
                          "and",as.numeric(refn[j,3])-(as.numeric(datap[1])-as.numeric(datap[6])),"nt downstream)")}
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else {datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Av",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="Av"
            }else {datap[7]=paste(datap[7],"Av",sep=";")}
          }else if (as.numeric(refn[j,2])>=(as.numeric(datap[1])-as.numeric(datap[6])) & 
                    (as.numeric(refn[j,3])<=(as.numeric(datap[1])+as.numeric(datap[6])))){
            if(datap[6]==0){
              a=paste("Ai overlaping gene ",refn[j,5])
            }else{a=paste("Ai overlaping gene ",refn[j,5])}
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else {datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("AA",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="AA"
            }else {datap[7]=paste(datap[7],"AA",sep=";")}
          }
          j=j+1
          if(j>nrow(refn)){
            arret_boucle=1
            j=j-1
          }
        }
        arret_boucle=0
        while((as.numeric(datap[1])-as.numeric(datap[6]))>=((as.numeric(refp[k,2]))-10000) & arret_boucle==0){
          if ((as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6])))<=250 & 
              (as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6])))>=0 & 
              (as.numeric(refp[k,2])>=(as.numeric(datap[1])+as.numeric(datap[6]))) ){
            if (datap[6]==0){
              a=paste((as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"nt upstream",refp[k,5])
            }else{a=paste("Between",(as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"and",
                          (as.numeric(refp[k,2])-(as.numeric(datap[1])+as.numeric(datap[6]))),"nt upstream ",refp[k,5])}
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("P",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="P"
            }else {datap[7]=paste(datap[7],"P",sep=";")}
          }else if ((as.numeric(refp[k,2])>(as.numeric(datap[1])-as.numeric(datap[6]))) & 
                    (as.numeric(refp[k,2])<(as.numeric(datap[1])+as.numeric(datap[6]))) & 
                    (as.numeric(refp[k,3])>(as.numeric(datap[1])+as.numeric(datap[6]))) ){
            if (datap[6]==0){
              a=paste((as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"nt upstream",refp[k,5])
            }else{a=paste("Between",(as.numeric(refp[k,2])-(as.numeric(datap[1])-as.numeric(datap[6]))),"and",
                          (as.numeric(refp[k,2])-(as.numeric(datap[1])+as.numeric(datap[6]))),"nt upstream ",refp[k,5])}
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("5U",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="5U"
            }else {datap[7]=paste(datap[7],"5U",sep=";")}
          }else if ((as.numeric(refp[k,2])<=(as.numeric(datap[1])-as.numeric(datap[6]))) & 
                    (as.numeric(refp[k,3])>=(as.numeric(datap[1])+as.numeric(datap[6])))){
            a=paste("Whithin ",refp[k,5],", ",((as.numeric(refp[k,2]))-(as.numeric(datap[1])-as.numeric(datap[6]))),
                    "nt after start and ",(as.numeric(refp[k,3])-(as.numeric(datap[1])+as.numeric(datap[6]))),"nt before end")
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("I",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="I"
            }else {datap[7]=paste(datap[7],"I",sep=";")}
          }else if ((as.numeric(refp[k,3])<(as.numeric(datap[1])+as.numeric(datap[6]))) & 
                    (as.numeric(refp[k,3])>(as.numeric(datap[1])-as.numeric(datap[6]))) & 
                    (as.numeric(refp[k,2])<(as.numeric(datap[1])-as.numeric(datap[6]))) ){
            a=paste("3UTR of ",refp[k,5],", ",((as.numeric(refp[k,3]))-(as.numeric(datap[1])+as.numeric(datap[6]))),
                    "nt before the end and ",((as.numeric(refp[k,3]))-(as.numeric(datap[1])-as.numeric(datap[6]))),"after")
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("3U",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="3U"
            }else {datap[7]=paste(datap[7],"3U",sep=";")}
          }else if ((as.numeric(refp[k,3])<=(as.numeric(datap[1])+as.numeric(datap[6]))) & 
                    (as.numeric(refp[k,2])>=(as.numeric(datap[1])-as.numeric(datap[6]))) ){
            a=paste("Overlaping gene ",refp[k,5])
            if (datap[8]==a) {
            }else if(datap[8]==0){
              datap[8]=a
            }else{datap[8]=paste(datap[8],a,sep=";")}
            cla=gregexpr("Chev",datap[7])[[1]][1]
            if(cla>=0){
            }else if (datap[7]==0){
              datap[7]="Chev"
            }else {datap[7]=paste(datap[7],"Chev",sep=";")}
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
  }#End annotation function
