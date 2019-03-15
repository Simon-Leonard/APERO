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


  max_dem_int=function(f,dem){ ### Find the biggest start between the actual start et the actual end
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
