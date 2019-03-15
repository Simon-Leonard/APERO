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
