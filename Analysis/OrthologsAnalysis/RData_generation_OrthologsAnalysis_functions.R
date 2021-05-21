




compute.tau <- function(exp){
	if(max(exp)==0){
	  return(NA)
	}
	n=length(exp)
	newexp=exp/max(exp)
	tau=sum(1-newexp)/(n-1)
	return(tau)
}


max.col <- function(exp){
	if(max(exp)==0){
	  return(NA)
	}
	maxcollist <- paste(names(exp[which(exp==max(exp))]), sep=":")
	return(maxcollist)
}

p.value.fromM8D <- function(d){
	pchisq(d, df=1, lower.tail=F)
}







