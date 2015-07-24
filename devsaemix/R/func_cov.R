subset.SaemixData<-function (x, subset, ...) {
    if (missing(subset)) 
        return(x)
    else {
        e <- substitute(subset)
	xdat<-x["data"]
        r <- eval(e, xdat, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    x1<-x
    x1["data"]<-x["data"][r,,drop=FALSE]
    if(length(x1["yorig"])>0) x1["yorig"]<-x["yorig"][r]
    id<-x1["data"][,x1["name.group"]]
    x1["N"]<-length(unique(id))
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-c(nind.obs[match(unique(id),names(nind.obs))])
    x1["nind.obs"]<-nind.obs
    x1["ntot.obs"]<-length(id)
    x1["data"]$index<-rep(1:x1["N"],times=nind.obs)
    return(x1)
}
