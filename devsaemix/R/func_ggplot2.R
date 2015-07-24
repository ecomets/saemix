
#######################	   Default plot options (list)	 ########################

saemix.plot.setoptions<-function(saemixObject) {
# setting default plot options
  plot.opt<-list(
# Options for plot types
    ilist=c(1:saemixObject["data"]["N"]),
    level=0:1,
    smooth=FALSE,
    line.smooth="s",
    indiv.par="map",			# type of individual parameters
    which.par="all",			# which parameters to plot 
    which.cov="all",			# which covariates to plot 
    which.resplot=c("res.vs.x","res.vs.pred","dist.qqplot","dist.hist"), # which type of residual plots
    which.pres=c("wres","npde"),	# which population weighted residuals
    which.poppred=c("ppred"),		# which population predictions to use (ypred=E(f(theta_i)), ppred=f(population parameters))
    indiv.histo=FALSE,			# whether to include an histogram of estimated individual parameters
    cov.value=rep(NA,length(saemixObject["model"]["name.cov"])),
# General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    interactive=FALSE,			# whether the user should be prompted before computing predictions or performing simulations for VPC, npde and wres
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main="",				# title
    xlab="",
    ylab="",
    col="black",
    pch=20,
    lty=1,
    lwd=1,
    xlim=c(),
    ylim=c(),
    xlog=FALSE,
    ylog=FALSE,
    type="b",
    cex=1,
    cex.axis=1,
    cex.lab=1,
    cex.main=1,
    obs.pch=20,
    pred.pch=20,
    obs.col="black",
    obs.lty=1,
    obs.lwd=1,
    obs.pch=20,
    ipred.col="black",
    ipred.lty=2,
    ipred.lwd=1,
    ipred.pch=20,
    ppred.col="black",
    ppred.lty=3,
    ppred.lwd=1,
    ppred.pch=20,
    pcol="black",
    lcol="black",
    fillcol="lightblue1",
    ablinecol="DarkRed",
    ablinelty=2,
    ablinelwd=2,
# 
    range=3,
    col.fillmed="pink",
    col.fillpi="slategray1",
    col.lmed="indianred4",
    col.lpi="slategray4",
    col.pobs="steelblue4",
    col.lobs="steelblue4",
    lty.lmed=2,
    lty.lpi=2,
    lwd.lmed=2,
    lwd.lpi=1,
    lwd.lobs=2,
    lty.lobs=1,
# Options for VPC plot
    vpc.method="equal",			# method (one of "equal"=same nb of points in each interval, "width"=equally spaced intervals (on the log-scale if xlog=TRUE), "user"=user-defined breaks, "optimal"=Marc's optimal binning algorithm); for "user", the breaks must be specified in vpc.breaks (otherwise defaults back to "equal"), while for the other methods the number of bins must be specified in vpc.bin
    vpc.breaks=NULL,			# user-defined breaks
    vpc.bin=10,				# nb of bins
    vpc.beta=0.2,			# value of beta used to compute the variance-based criterion (Jopt,beta(I)) in the clustering algorithm
    vpc.lambda=0.3,			# value of lambda used in the penalised criterion to select the number of bins (if vpc.bin=NULL)
    vpc.interval=0.95,
    vpc.pi=TRUE,
    vpc.obs=TRUE)
    
     if(is.null(plot.opt$name.X)) {
        if(length(saemixObject["data"]["name.X"])>0) plot.opt$name.X<-saemixObject["data"]["name.X"] else plot.opt$name.X<-saemixObject["data"]["name.predictors"][1]
    }
    plot.opt$xlab<-paste(plot.opt$name.X," (",saemixObject["data"]["units"]$x,")", sep="")
     if(length(saemixObject["data"]["name.response"])>0)
    plot.opt$ylab<-paste(saemixObject["data"]["name.response"]," (", saemixObject["data"]["units"]$y,")",sep="")
   return(plot.opt)
}

#################    Function to supersede default plot options	 ##################

replace.plot.options<-function(plot.opt,...) {
  args1<-match.call(expand.dots=TRUE)
  if(length(args1)>2) {
# General arguments: col, pch
    i1<-match("col",names(args1))
    if(!is.na(i1)) {
      plot.opt$col<-eval(args1[[i1]])
      plot.opt$obs.col<-eval(args1[[i1]])
      plot.opt$ipred.col<-eval(args1[[i1]])
      plot.opt$ppred.col<-eval(args1[[i1]])
      plot.opt$pcol<-eval(args1[[i1]])
      plot.opt$lcol<-eval(args1[[i1]])
    }
    i1<-match("pch",names(args1))
    if(!is.na(i1)) {
      plot.opt$pch<-eval(args1[[i1]])
      plot.opt$obs.pch<-eval(args1[[i1]])
      plot.opt$ipred.pch<-eval(args1[[i1]])
      plot.opt$ppred.pch<-eval(args1[[i1]])
    }
# Other arguments
    for(i in 3:length(args1)) {
      if(match(names(args1)[i],names(plot.opt),nomatch=0)>0)    
#    plot.opt[[names(args1)[i]]]<-args1[[i]] else {
    plot.opt[[names(args1)[i]]]<-eval(args1[[i]]) else {
      if(names(args1)[i]!="plot.type") cat("Argument",names(args1)[i],"not available, check spelling.\n")
    }
   }
  }
  return(plot.opt)
}

####################################################################################