####################################################################################
####		Defining class containing data, model and options		####
####################################################################################

###############################
# definition
setClass(Class="SaemixObject",
  representation=representation(
    data="SaemixData",		# Data
    model="SaemixModel",	# Model
    results="SaemixRes",	# Fit results
    rep.data="SaemixRepData",	# Data replicates during algorithm (nb chains)
    sim.data="SaemixSimData", 	# Simulated data
    options="list",		# Options and parameters for algorithm
    prefs="list"		# Options for graphs
  ),
  validity=function(object){
#    cat ("--- Checking SaemixObject object ---\n")
    validObject(object@data)
    validObject(object@model)
    return(TRUE)
  }
)

###############################
# Initialize
setMethod(
  f="initialize",
  signature="SaemixObject",
  definition= function (.Object,data,model,options=list()){
#    cat ("--- initialising SaemixObject --- \n")
    if(model["error.model"]=='exponential') {
      yobs<-data["data"][,data["name.response"]]
      y<-log(cutoff(yobs))
      data["yorig"]<-yobs
      data["data"][,data["name.response"]]<-y
    }
    .Object@data<-data
# Adjusting number of covariates
    if(dim(model@covariate.model)[1]>length(data@name.covariates)) {
      cat("The number of covariates in model (",dim(model@covariate.model)[1],") is larger than the number of covariates in the dataset (",length(data@name.covariates),"), keeping only the first",length(data@name.covariates),":",data@name.covariates,".\n")
      model@covariate.model<-model@covariate.model[1:length(data@name.covariates),]
    }
    if(dim(model@covariate.model)[1]<length(data@name.covariates) & dim(model@covariate.model)[1]>0) {
    	cat("The number of covariates in model (",dim(model@covariate.model)[1],") is smaller than the number of covariates in the dataset (",length(data@name.covariates),"), assuming no covariate-parameter relationship for the remaining covariates; please check covariates:",data@name.covariates,".\n")
    	l1<-rep(0,dim(model@covariate.model)[2])
    	n1<-length(data@name.covariates)-dim(model@covariate.model)[1]
    	model@covariate.model<-rbind(model@covariate.model, matrix(rep(l1,n1),nrow=n1))
    }
# setting the names of the fixed effects
    if(dim(model@covariate.model)[1]>0) {
      nam.with.cov<-rep("",length(model@covariate.model))
      row1<-matrix(rep(model@name.modpar,length(data@name.covariates)), ncol=length(model@name.modpar),byrow=TRUE)
      col1<-matrix(rep(data@name.covariates,length(model@name.modpar)), ncol=length(model@name.modpar))
      idcov<-which(model@covariate.model==1)
      nam.with.cov[idcov]<-paste("beta_",col1[idcov],"(",row1[idcov],")",sep="")
      nam1<-rbind(model@name.modpar,matrix(nam.with.cov, ncol=length(model@name.modpar)))
      nam1<-c(nam1)
      model@name.fixed<-nam1[nam1!=""]
    } else model@name.fixed<-model@name.modpar
    i1.omega2<-model@indx.omega
    model@name.random<-paste("omega2",model@name.modpar[model@indx.omega], sep=".")
    .Object@model<-model
    .Object@model@betaest.model<-matrix(c(rep(1,.Object@model@nb.parameters), c(t(.Object@model@covariate.model))),ncol=.Object@model@nb.parameters,byrow=TRUE)
    
# Covariates
    .Object@model@name.cov<-.Object@data@name.covariates
    if(length(.Object@model@name.cov)>0 & sum(.Object@model@covariate.model)>0) {
      try(rownames(.Object@model@covariate.model)<-.Object@model@name.cov)
      try(rownames(.Object@model@betaest.model)[2:(1+ length(.Object@model@name.cov))]<-.Object@model@name.cov)
    }
    ucov <- rownames(.Object@model@covariate.model)[ rowSums(.Object@model@covariate.model)>0]
    if(length(ucov)>0) {
      for(icov in ucov) {
      	xdat<-subset(.Object@data@data,is.na(get(icov)))
      	if(dim(xdat)[1]>0) {
      		imis<-unique(xdat[,.Object@data@name.group])
      		cat("Missing values for covariate",as.character(icov),"for which a parameter-covariate relationship is estimated: removing subject(s)",imis,"from the dataset.\n")
      	}
      	.Object@data<-subset(.Object@data,!is.na(get(icov)))
      }
    }
# Initialising options
    opt<-saemixControl()
    if(length(options)>0) {
    for(i in names(options)) opt[i]<-options[i]
      if(!opt$fix.seed) {
      rm(.Random.seed)
      runif(1)
      opt$seed<-.Random.seed[5]
      }
    if(is.null(options$nb.chains) & data@N>0) opt$nb.chains<-ceiling(50/data@N)
    if(data@N>0 && data@N<50 & opt$nb.chains<ceiling(50/data@N)) {
      cat("The number of subjects is small, increasing the number of chains to", ceiling(50/data@N),"to improve convergence\n")
      opt$nb.chains<-ceiling(50/data@N)
    }
    if(opt$ipar.lmcmc<2) {
      opt$ipar.lmcmc<-2
      cat("Value of L_MCMC too small, setting it to 2 (computation of the conditional means and variances of the individual parameters)\n")
    }
    }
    opt$nbiter.sa<-round(opt$nbiter.saemix[1]/2)
    opt$nbiter.tot<-sum(opt$nbiter.saemix)
    .Object@options<-opt
# Options for plots
    .Object@prefs<-saemix.plot.setoptions(.Object)
# Object validation
    validObject(.Object)
    return (.Object )
  }
)

###########################	Default options		#############################

#' List of options for running the algorithm SAEM
#' 
#' List containing the variables relative to the optimisation algorithm. All
#' these elements are optional and will be set to default values when running
#' the algorithm if they are not specified by the user.
#' 
#' All the variables are optional and will be set to their default value when
#' running \code{\link{saemix}}.
#' 
#' The function \code{\link{saemix}} returns an object with an element options
#' containing the options used for the algorithm, with defaults set for
#' elements which have not been specified by the user.
#' 
#' These elements are used in subsequent functions and are not meant to be used
#' directly.
#' 
#' @param algorithms a vector of 1s specifying which algorithms are to be run.
#' Defaults to c(1,1,1) (respectively estimation of the individual parameters
#' (MAP estimates), estimation of the Fisher Information Matrix and
#' log-likelihood by linearisation, and estimation of the log-likelihood by
#' importance sampling)
#' @param nbiter.saemix nb of iterations in each step (a vector containing 2
#' elements)
#' @param nb.chains nb of chains to be run in parallel in the MCMC algorithm.
#' Defaults to 1.
#' @param nbiter.burn nb of iterations for burning
#' @param nbiter.mcmc nb of iterations in each kernel during the MCMC step
#' @param proba.mcmc probability of acceptance
#' @param stepsize.rw stepsize for kernels q2 and q3. Defaults to 0.4
#' @param rw.init initial variance parameters for kernels. Defaults to 0.5
#' @param alpha.sa parameter controlling cooling in the Simulated Annealing
#' algorithm. Defaults to 0.97
#' @param fix.seed TRUE (default) to use a fixed seed for the random number
#' generator. When FALSE, the random number generator is initialised using a
#' new seed, created from the current time.  Hence, different sessions started
#' at (sufficiently) different times will give different simulation results.
#' The seed is stored in the element seed of the options list.
#' @param seed seed for the random number generator. Defaults to 123456
#' @param nmc.is nb of samples used when computing the likelihood through
#' importance sampling
#' @param nu.is number of degrees of freedom of the Student distribution used
#' for the estimation of the log-likelihood by Importance Sampling. Defaults to
#' 4
#' @param print.is when TRUE, a plot of the likelihood as a function of the
#' number of MCMC samples when computing the likelihood through importance
#' sampling is produced and updated every 500 samples. Defaults to FALSE
#' @param nbdisplay nb of iterations after which to display progress
#' @param displayProgress when TRUE, the convergence plots are plotted after
#' every nbdisplay iteration, and a dot is written in the terminal window to
#' indicate progress. When FALSE, plots are not shown and the algorithm runs
#' silently. Defaults to TRUE
#' @param nnodes.gq number of nodes to use for the Gaussian quadrature when
#' computing the likelihood with this method (defaults to 12)
#' @param nsd.gq span (in SD) over which to integrate when computing the
#' likelihood by Gaussian quadrature. Defaults to 4 (eg 4 times the SD)
#' @param maxim.maxiter Maximum number of iterations to use when maximising the
#' fixed effects in the algorithm. Defaults to 100
#' @param nb.sim number of simulations to perform to produce the VPC plots or
#' compute npde. Defaults to 1000
#' @param nb.simpred number of simulations used to compute mean predictions
#' (ypred element), taken as a random sample within the nb.sim simulations used
#' for npde
#' @param ipar.lmcmc number of iterations required to assume convergence for
#' the conditional estimates. Defaults to 50
#' @param ipar.rmcmc confidence interval for the conditional mean and variance.
#' Defaults to 0.95
#' @param print whether the results of the fit should be printed out. Defaults
#' to TRUE
#' @param save whether the results of the fit should be saved to a file.
#' Defaults to TRUE
#' @param save.graphs whether diagnostic graphs and individual graphs should be
#' saved to files. Defaults to TRUE
#' @param directory the directory in which to save the results. Defaults to
#' "newdir" in the current directory
#' @param warnings whether warnings should be output during the fit. Defaults
#' to FALSE
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}}, \code{\link{saemix}}
#' @references Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear
#' mixed effects models. Computational Statistics and Data Analysis 49, 4
#' (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords models
#' @examples
#' 
#' 
#' # All default options
#' saemix.options<-saemixControl()
#' 
#' # All default options, changing seed
#' saemix.options<-saemixControl(seed=632545)
#' 
#' 
#' @export saemixControl
saemixControl<-function(map=TRUE,fim=TRUE,ll.is=TRUE,ll.gq=FALSE,nbiter.saemix=c(300,100), nb.chains=1,fix.seed=TRUE,seed=23456,nmc.is=5000,nu.is=4, print.is=FALSE,nbdisplay=100,displayProgress=TRUE,nbiter.burn=5, nbiter.mcmc=c(2,2,2),proba.mcmc=0.4,stepsize.rw=0.4,rw.init=0.5,alpha.sa=0.97,  nnodes.gq=12,nsd.gq=4,maxim.maxiter=100,nb.sim=1000,nb.simpred=100, ipar.lmcmc=50,ipar.rmcmc=0.05, print=TRUE, save=TRUE, save.graphs=TRUE,directory="newdir",warnings=FALSE) {
  if(fix.seed) seed<-seed else {
    rm(.Random.seed)
    runif(1)
    seed<-.Random.seed[5]
  }
  if(ipar.lmcmc<2) {
    ipar.lmcmc<-2
    cat("Value of L_MCMC too small, setting it to 2 (computation of the conditional means and variances of the individual parameters)\n")
  }
  list(map=map,fim=fim,ll.is=ll.is,ll.gq=ll.gq,nbiter.saemix=nbiter.saemix, nbiter.burn=nbiter.burn,nb.chains=nb.chains,fix.seed=fix.seed,seed=seed, nmc.is=nmc.is,nu.is=nu.is,print.is=print.is, nbdisplay=nbdisplay,displayProgress=displayProgress,print=print,save=save, save.graphs=save.graphs,directory=directory,warnings=warnings, nbiter.mcmc=nbiter.mcmc,proba.mcmc=proba.mcmc,stepsize.rw=stepsize.rw, rw.init=rw.init,alpha.sa=alpha.sa,nnodes.gq=nnodes.gq,nsd.gq=nsd.gq, maxim.maxiter=maxim.maxiter,nb.sim=nb.sim,nb.simpred=nb.simpred,
ipar.lmcmc=ipar.lmcmc,ipar.rmcmc=ipar.rmcmc)
}

####################################################################################
####			saemixObject class - accesseur				####
####################################################################################

##' Get/set methods for SaemixObject object
##' 
##' Access slots of a SaemixObject using the object["slot"] format
##' 
##' @keywords methods
##' @exportMethod [
##' @exportMethod [<-

# Getteur
setMethod(
  f ="[",
  signature = "SaemixObject" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "data"={return(x@data)},
    "model"={return(x@model)},
    "results"={return(x@results)},
    "rep.data"={return(x@rep.data)},
    "sim.data"={return(x@sim.data)},
    "options"={return(x@options)},
    "prefs"={return(x@prefs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixObject" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "data"={x@data<-value},
    "model"={x@model<-value},
    "results"={x@results<-value},
    "rep.data"={x@rep.data<-value},
    "sim.data"={x@sim.data<-value},
    "options"={x@options<-value},
    "prefs"={x@prefs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)


####################################################################################
####				Summary method for SaemixObject			####
####################################################################################
setMethod("summary","SaemixObject",
  function(object, print=TRUE, ...) {
    if(length(object@results@fixed.effects)==0) {
      cat("Object of class SaemixObject, no fit performed yet.\n")
      return()
    }
    digits<-2;nsmall<-2
    if(print) {
    cat("----------------------------------------------------\n")
    cat("-----------------  Fixed effects  ------------------\n")
    cat("----------------------------------------------------\n")
    }
    if(length(object@results@se.fixed)==0) {
      tab<-data.frame(c(object@results@name.fixed, object@results@name.res[object@results@indx.res]), c(object@results@fixed.effects,object@results@respar[object@results@indx.res]))
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-data.frame(c(object@results@name.fixed, object@results@name.res[object@results@indx.res]), c(object@results@fixed.effects,object@results@respar[object@results@indx.res]),c(object@results@se.fixed,object@results@se.respar[object@results@indx.res]), stringsAsFactors=FALSE)
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
      if(length(object@results@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[object@results@indx.cov]<-1-normcdf(abs(wstat[object@results@indx.cov]))
      tab<-cbind(tab,"p-value"=pval,stringsAsFactors=FALSE)
      }
    }
    tab.fix<-tab
    for(i in 2:dim(tab)[2]) {
     xcol<-as.double(as.character(tab[,i]))
     idx<-which(!is.na(xcol) & xcol!="-")
     tab[idx,i]<-format(xcol[idx],digits=digits,nsmall=nsmall)
    }
    if(print) {
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("-----------  Variance of random effects  -----------\n")
    cat("----------------------------------------------------\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n")
    }
    if(length(object@results@se.omega)==0) {
      tab<-data.frame(object@results@name.random, diag(object@results@omega)[object@results@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-data.frame(object@results@name.random, diag(object@results@omega)[object@results@indx.omega], object@results@se.omega[object@results@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
    }
    tab.random<-tab
    for(i in 2:dim(tab)[2]) 
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
    if(print) {
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("------  Correlation matrix of random effects  ------\n")
    cat("----------------------------------------------------\n")
    }
    tab<-cov2cor(object@results@omega[object@results@indx.omega, object@results@indx.omega,drop=FALSE])
    tab.corr<-tab
    for(i in 1:dim(tab)[2])
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits,nsmall=nsmall)
    try(colnames(tab)<-rownames(tab)<-object@results@name.random)
    if(print) print(tab,quote=FALSE)
    l1<-rep(NA,3)
    tab.ll<-data.frame(Method=c("Linearisation","Importance Sampling","Gaussian Quadrature"),"-2xLL"=l1,AIC=l1,BIC=l1,BIC.h=l1)
    if(length(object@results@ll.lin)>0 | length(object@results@ll.is)>0 | length(object@results@ll.gq)>0) {
    	if(print) {
    cat("----------------------------------------------------\n")
    cat("---------------  Statistical criteria  -------------\n")
    cat("----------------------------------------------------\n")
    	}
    if(length(object@results@ll.lin)>0) {
    	if(print) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL  =",(-2*object@results@ll.lin),"\n")
    cat("      AIC   =",object@results@aic.lin,"\n")
    cat("      BIC   =",object@results@bic.lin,"\n")
    cat("      BIC.h =",object@results@bic.h.lin,"\n")
    	}
    tab.ll[1,2:5]<-c((-2*object@results@ll.lin),object@results@aic.lin, object@results@bic.lin, object@results@bic.h.lin)
    }
    if(length(object@results@ll.is)>0) {
    	if(print) {
    cat("\nLikelihood computed by importance sampling\n")
    cat("      -2LL  =",(-2*object@results@ll.is),"\n")
    cat("      AIC   =",object@results@aic.is,"\n")
    cat("      BIC   =",object@results@bic.is,"\n")
    cat("      BIC.h =",object@results@bic.h.is,"\n")
    	}
    tab.ll[2,2:5]<-c((-2*object@results@ll.is),object@results@aic.is, object@results@bic.is, object@results@bic.h.is)
    }  
    if(length(object@results@ll.gq)>0) {
    	if(print) {
    cat("\nLikelihood computed by Gaussian quadrature\n")
    cat("      -2LL  =",(-2*object@results@ll.gq),"\n")
    cat("      AIC   =",object@results@aic.gq,"\n")
    cat("      BIC   =",object@results@bic.gq,"\n")
    cat("      BIC.h =",object@results@bic.h.gq,"\n")
    	}
    tab.ll[3,2:5]<-c((-2*object@results@ll.gq),object@results@aic.gq, object@results@bic.gq, object@results@bic.h.gq)
    }
    if(print) cat("----------------------------------------------------\n")
    }
    tab<-data.frame(Id=unique(object@data@data[,object@data@name.group]), object@results@cond.mean.psi)
    try(colnames(tab)[-c(1)]<-object@model@name.modpar,silent=TRUE)
    npar<-length(object@results@name.fixed)
    coef<-list(fixed=tab.fix[1:npar,2],random=list(map.psi=object@results@map.psi, cond.mean.psi=tab))
    sigma<-tab.fix[-c(1:npar),2]

    res<-list(fixed.effects=tab.fix,sigma=sigma,random.effects=tab.random, correlation.matrix=tab.corr,logLik=tab.ll,coefficients=coef)
    if(length(object@results@fim)>0) res$FIM<-object@results@fim
    res$data<-list(N=object@data@N,nobs=list(ntot=object@data@ntot.obs, nind=object@data@nind.obs),data=object@data@data)
    if(length(object@results@ypred)>0 | length(object@results@ipred)>0  | length(object@results@ppred)>0 | length(object@results@icpred)>0) {
      res$fitted<-list(population=list(pop.param=object@results@ppred, pop.mean=object@results@ypred),individual=list(map.ipred=object@results@ipred, cond.ipred=object@results@icpred))
    }
     if(length(object@results@wres)>0 | length(object@results@iwres)>0  | length(object@results@icwres)>0 | length(object@results@pd)>0) {
      res$residuals<-list(population=list(wres=object@results@wres), individual=list(map.iwres=object@results@iwres,cond.iwres=object@results@icwres, pd=object@results@pd, npde=object@results@npde))
    }
   
    invisible(res)
 }
)

####################################################################################
####			Print and show methods for SaemixObject			####
####################################################################################

setMethod("print","SaemixObject",
  function(x,nlines=10,...) {
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    print(x@data,nlines=nlines)
    cat("-----------------------------------\n")
    cat("----          Model            ----\n")
    cat("-----------------------------------\n")
    print(x@model)
    cat("-----------------------------------\n")
    cat("----    Key algorithm options  ----\n")
    cat("-----------------------------------\n")
    if(x@options$map) cat("    Estimation of individual parameters (MAP)\n")
    if(x@options$fim) cat("    Estimation of standard errors and linearised log-likelihood\n")
    if(x@options$ll.is) cat("    Estimation of log-likelihood by importance sampling\n")
    if(x@options$ll.gq) cat("    Estimation of log-likelihood by gaussian quadrature\n")
    if(as.integer(x@options$map+x@options$fim+x@options$ll.is+x@options$ll.gq)==0) cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),x@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of iterations: ",st1,"\n")
    cat("    Number of chains: ",x@options$nb.chains,"\n")
    cat("    Seed: ",x@options$seed,"\n")
    if(x@options$ll.is) cat("    Number of MCMC iterations for IS: ",x@options$nmc.is,"\n")
    cat("    Simulations:\n")
    cat("        nb of simulated datasets used for npde: ",x@options$nb.sim,"\n")
    cat("        nb of simulated datasets used for VPC: ",x@options$nb.simpred,"\n")
    cat("    Input/output\n")
    cat("        save the results to a file: ",x@options$save,"\n")
    cat("        save the graphs to files: ",x@options$save.graphs,"\n")
    if(x@options$save | x@options$save.graphs) cat("        directory where results should be saved: ",x@options$directory,"\n")
    cat("----------------------------------------------------\n")
    cat("----                  Results                   ----\n")
    print(x@results)
  }
)

setMethod("show","SaemixObject",
  function(object) {
#    cat("Object of class SaemixObject\n")
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------------\n")
    cat("----         Data and Model          ----\n")
    cat("-----------------------------------------\n")
#    show(object@data)
    cat("Data\n")
    cat("    Dataset",object@data@name.data,"\n")
    st1<-paste(object@data@name.response," ~ ",paste(object@data@name.predictors, collapse=" + ")," | ", object@data@name.group,sep="")
    cat("    Longitudinal data:",st1,"\n\n")
#    show(object@model)
    cat("Model:\n")
    if(length(object@model@description)>0 && nchar(object@model@description)>0) {
      cat("   ",object@model@description,"\n")}
    fix1<-ifelse(object@model@fixed.estim==1,""," [fixed]")
    cat("    ",object@model@nb.parameters,"parameters:", paste(object@model@name.modpar,fix1,sep=""),"\n")
    cat("     error model:",object@model@error.model,"\n")
    if(dim(object@model@covariate.model)[1]>0) {
      cat("     covariate model:\n")
      print(object@model@covariate.model) 
    } else cat("     No covariate\n")
    cat("\n")
    cat("Key options\n")
    if(object@options$map) cat("    Estimation of individual parameters (MAP)\n")
    if(object@options$fim) cat("    Estimation of standard errors and linearised log-likelihood\n")
    if(object@options$ll.is) cat("    Estimation of log-likelihood by importance sampling\n")
    if(object@options$ll.gq) cat("    Estimation of log-likelihood by gaussian quadrature\n")
    if(as.integer(object@options$map+object@options$fim+object@options$ll.is+object@options$ll.gq)==0) cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),object@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of iterations: ",st1,"\n")
    cat("    Number of chains: ",object@options$nb.chains,"\n")
    cat("    Seed: ",object@options$seed,"\n")
    if(object@options$ll.is) cat("    Number of MCMC iterations for IS: ",object@options$nmc.is,"\n")
    cat("    Input/output\n")
    if(object@options$save)
      cat("        save the results to a file: ",object@options$save,"\n") else 
      cat("        results not saved\n")
    if(object@options$save.graphs)
      cat("        save the graphs to files: ",object@options$save.graphs,"\n") else 
      cat("        no graphs\n")
    if(object@options$save | object@options$save.graphs) cat("        directory where results are saved: ",object@options$directory,"\n")
    if(FALSE) {
    if(length(object@rep.data)>0) irep<-1 else irep<-0
    if(length(object@sim.data)>0) isim<-1 else isim<-0
    if(irep>0 | isim>0) {
      cat("-----------------------------------\n")
      cat("----      Other components     ----\n")
      cat("-----------------------------------\n")
      if(irep>0) cat("    Replicated data on ",object@options$nb.chains," chains\n")
      if(isim>0) cat("    Simulated data, ",object@options$nb.sim," simulations\n")
    }
    cat("-----------------------------------\n")
  }
  if(length(object@results@fixed.effects)>0) {
    cat("----------------------------------------------------\n")
    cat("----                  Results                   ----\n")
    show(object@results)
  }}
)

# Could be print, with only head of data
setMethod("showall","SaemixObject",
  function(object) {
#    cat("Object of class SaemixObject\n")
    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    showall(object@data)
    cat("-----------------------------------\n")
    cat("----          Model            ----\n")
    cat("-----------------------------------\n")
    show(object@model)
    cat("-----------------------------------\n")
    cat("----      Algorithm options    ----\n")
    cat("-----------------------------------\n")
    if(object@options$map) cat("    Estimation of individual parameters (MAP)\n")
    if(object@options$fim) cat("    Estimation of standard errors and linearised log-likelihood\n")
    if(object@options$ll.is) cat("    Estimation of log-likelihood by importance sampling\n")
    if(object@options$ll.gq) cat("    Estimation of log-likelihood by gaussian quadrature\n")
    if(as.integer(object@options$map+object@options$fim+object@options$ll.is+object@options$ll.gq)==0) cat("    Algorithms: estimation only\n")
    st1<-paste(c("K1=","K2="),object@options$nbiter.saemix,sep="",collapse=", ")
    cat("    Number of chains: ",object@options$nb.chains,"\n")
    cat("    Number of iterations: ",st1,"\n")
    cat("        nb of iterations for SA: ",object@options$nbiter.sa,"\n")
    cat("        nb of burning iterations: ",object@options$nbiter.burn,"\n")
    cat("    Seed:\n")
    cat("        setting a random seed: ",!object@options$fix.seed,"\n")
    cat("        seed for the random number generator: ",object@options$seed,"\n")
    if(object@options$ll.is) {
    	cat("    Estimation of LL by Importance Sampling:\n")
    	cat("        number of MCMC samples: ",object@options$nmc.is,"\n")
    	cat("        nu for IS: ",object@options$nu.is,"\n")
    	cat("        produce plots during the estimation of LL by IS: ",object@options$print.is,"\n")
    }
    cat("    Input/output\n")
    cat("        display progress during the estimation process: ",object@options$displayProgress,"\n")
    cat("        nb of iterations after which to display progress: ",object@options$nbdisplay,"\n")
    cat("        print out the results after the fit: ",object@options$print,"\n")
    cat("        save the results to a file: ",object@options$save,"\n")
    cat("        save the graphs to files: ",object@options$save.graphs,"\n")
    if(object@options$save | object@options$save.graphs) cat("        directory where results should be saved: ",object@options$directory,"\n")
    cat("        whether warnings should be output during the fit: ",object@options$warnings,"\n")
    cat("    SAEM algorithm\n")
    cat("        number of MCMC iterations for each kernel: ",object@options$nbiter.mcmc,"\n")
    cat("        probability of acceptance: ",object@options$proba.mcmc,"\n")
    cat("        : ",object@options$stepsize.rw,"\n")
    cat("        : ",object@options$rw.init,"\n")
    cat("        : ",object@options$alpha.sa,"\n")
    cat("        maximum nb of iterations for estimation of fixed effects: ",object@options$maxim.maxiter,"\n")
    if(object@options$ll.gq) {
    	cat("    Estimation of LL by Gaussian Quadrature:\n")
    	cat("        number of nodes: ",object@options$nnodes.gq,"\n")
    	cat("        width of integral: ",object@options$nsd.gq,"\n")
    }
    cat("    Simulations:\n")
    cat("        nb of simulated datasets used for npde: ",object@options$nb.sim,"\n")
    cat("        nb of simulated datasets used for VPC: ",object@options$nb.simpred,"\n")
    cat("    Estimation of individual parameters\n")
    cat("        nb of iterations: ",object@options$ipar.lmcmc,"\n")
    cat("        : ",object@options$ipar.rhomcmc,"\n")
    cat("        size of confidence interval: ",object@options$ipar.rmcmc,"\n")
    cat("-----------------------------------\n")
    if(length(object@results@fixed.effects)>0) {
      cat("----------------------------------------------------\n")
      cat("----                  Results                   ----\n")
      show(object@results)
    }
  }
)

####################################################################################
####			SaemixObject class - method to predict			####
####################################################################################

setMethod(f="predict",
  signature="SaemixObject",
  def=function(object,newdata=NULL,type=c("ipred", "ypred", "ppred", "icpred"), se.fit=FALSE, ...) {
    type<-match.arg(type)
#    se.fit<-match.arg(se.fit) # doesn't work with logical type, change
#    if(se.fit) cat("Currently predict() does not handle argument se.fit=TRUE.\n")
    if(missing(newdata)) {
      xpred<-fitted(object,type,...)
      if(length(xpred)==0) {
        namObj<-deparse(substitute(object))
        opred<-saemix.predict(object)
        assign(namObj,opred,envir=parent.frame()) # update object invisibly with predictions
        xpred<-fitted(object,type,...)
      }
    } else {
# Ignore type - when newdata is given, predictions correspond to the population predictions using the final estimates
      id<-newdata["data"][,newdata["name.group"]]
      if(length(newdata["name.covariates"])==0) tab<-data.frame(id=id) else
        tab<-data.frame(id=id,newdata["data"][, newdata["name.covariates",drop=FALSE]])
      # Mcovariates: extract columns needed in the model
      temp2<-unique(tab)
      temp<-tab[!duplicated(id),,drop=FALSE]
      if(dim(temp)[1]!=dim(temp2)[1]) cat("Some covariates have time-varying values; only the first is taken into account in the current version of the algorithm.\n")
      #temp<-temp[order(temp[,1]),]
      
      if(length(newdata["name.covariates"])>0) {
        Mcovariates<-data.frame(id=rep(1,newdata["N"]),temp[,2:dim(temp)[2]])} else {
          Mcovariates<-data.frame(id=rep(1,newdata["N"]))
        }
      j.cov<-which(rowSums(object["model"]["betaest.model"])>0)
      names(j.cov)[1]<-names(Mcovariates)[1]
      Mcovariates<-Mcovariates[,names(j.cov),drop=FALSE] # eliminate all the unused covariates
      for(icol in dim(Mcovariates)[2])
        if(is.factor(Mcovariates[,icol])) Mcovariates[,icol]<-as.numeric(Mcovariates[,icol])-1
      # COV: design matrix
      COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
      for(j in 1:-object["model"]["nb.parameters"]) {
        jcov<-which(object["model"]["betaest.model"][,j]==1)
        aj<-as.matrix(Mcovariates[,jcov])
        COV<-cbind(COV,aj)
      }
      # population parameters given by phi=Ci*mu (COV %*% Lcovariates) and psi=h(phi)
      pop.phi<-COV %*% object["results"]["MCOV"]
      pop.psi<-transphi(pop.phi,object["model"]["transform.par"])
      structural.model<-object["model"]["model"]
      chdat<-new(Class="SaemixRepData",data=newdata, nb.chains=1)
      IdM<-chdat["dataM"]$IdM
      XM<-chdat["dataM"][,newdata["name.predictors"],drop=FALSE]
      # Model predictions with pop.psi
      fpred<-structural.model(pop.psi, IdM, XM)
      if(object["model"]["error.model"]=="exponential")
        fpred<-log(cutoff(fpred))
      xpred<-fpred
    }
    xpred
  }
)

saemix.predict<-function(object) {
  if(length(object["results"]["map.psi"])==0)
    object<-map.saemix(object)
  # en principe n'arrive jamais car on les calcule pdt le fit...
  if(length(object["results"]["cond.mean.phi"])==0)
    object<-conddist.saemix(object)
  saemix.res<-object["results"]
  xind<-object["data"]["data"][,object["data"]["name.predictors"],drop=FALSE]
  if(object["model"]["error.model"]=="exponential") yobs<-object["data"]["yorig"] else yobs<-object["data"]["data"][,object["data"]["name.response"]]
  index<-object["data"]["data"][,"index"]
  # Individual predictions
  ipred<-object["model"]["model"](saemix.res["map.psi"][, 2:dim(saemix.res["map.psi"])[2]],index,xind)
  ires<-object["data"]["data"][,object["data"]["name.response"]]-ipred
  psiM<-transphi(saemix.res["cond.mean.phi"],object["model"]["transform.par"])
  icond.pred<-object["model"]["model"](psiM,index,xind)
  saemix.res["ipred"]<-ipred
  saemix.res["icpred"]<-icond.pred
  # Individual weighted residuals
  pres<-saemix.res["respar"]
  gpred<-error(ipred,pres)
  iwres<-(ipred-yobs)/gpred
  gpred<-error(icond.pred,pres)
  icwres<-(icond.pred-yobs)/gpred
  saemix.res["iwres"]<-iwres
  saemix.res["icwres"]<-icwres
  # Population predictions using the population parameters [ f(mu) ]
  psiM<-transphi(saemix.res["mean.phi"],object["model"]["transform.par"])
  ppred<-object["model"]["model"](psiM,index,xind)
  saemix.res["ppred"]<-ppred
  if(length(saemix.res["predictions"])==0) 
    saemix.res["predictions"]<-data.frame(ppred=ppred,ipred=ipred,icpred=icond.pred,ires=ires,iwres=iwres,icwres=icwres) else {
      saemix.res["predictions"]$ppred<-ppred
      saemix.res["predictions"]$ipred<-ipred
      saemix.res["predictions"]$icpred<-icond.pred
      saemix.res["predictions"]$ires<-ires
      saemix.res["predictions"]$iwres<-iwres
      saemix.res["predictions"]$icwres<-icwres
    }
  # Population weighted residuals: needs the individual variance-covariance matrix => use compute.sres to estimate these by simulations
  object["results"]<-saemix.res
  return(object)
}

####################################################################################
####			SaemixObject class - method to plot			####
####################################################################################

setMethod(f="plot",
  signature="SaemixObject",
  def=function(x,y,...) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("plot.type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"reduced"
#    cat("plot.type=",plot.type,"\n")
    if(plot.type[1]=="reduced") plot.type<-c("data","convergence","likelihood", "observations.vs.predictions")
    if(plot.type[1]=="full") plot.type<-c("data","convergence","likelihood", "observations.vs.predictions","residuals.scatter","residuals.distribution","vpc")
    
    pltyp<-c("data","convergence","likelihood","individual.fit", "population.fit", "both.fit","observations.vs.predictions","residuals.scatter", "residuals.distribution","vpc","npde","random.effects","marginal.distribution", "correlations","parameters.vs.covariates","randeff.vs.covariates")
    ifnd<-pmatch(plot.type,pltyp)
    if(sum(is.na(ifnd))>0) {
      cat("The following plot types were not found or are ambiguous:", plot.type[is.na(ifnd)],"\n")
    }
    ifnd<-ifnd[!is.na(ifnd)]
    if(length(ifnd)==0) return()
    plot.type<-pltyp[ifnd]
    interactive<-x["prefs"]$interactive
    id.pred<-match(plot.type,c("observations.vs.predictions","individual.fit", "residuals.scatter","residuals.distribution"))
    if(x@prefs$which.poppred=="ppred") id.pred<-c(id.pred,match(plot.type, c("population.fit", "both.fit")))
    id.map<-match(plot.type,c("randeff.vs.covariates","parameters.vs.covariates"))
    id.sim<-match(plot.type, c("vpc"))
    id.res<-match(plot.type, c("npde","residuals.scatter", "residuals.distribution"))
    if(x@prefs$which.poppred=="ypred") id.sim<-c(id.sim,match(plot.type, c("population.fit", "both.fit")))
    id.pred<-id.pred[!is.na(id.pred)]
    id.sim<-id.sim[!is.na(id.sim)]
    id.map<-id.map[!is.na(id.map)]
    id.res<-id.res[!is.na(id.res)]
    namObj<-deparse(substitute(x))
#    cat(namObj,"\n")
    if(length(id.pred)>0) {
      if(length(x["results"]["ipred"])==0 | length(x["results"]["iwres"])) {
        boolpred<-TRUE
        if(interactive) {
      cok<-readline(prompt="Computations will be performed to obtain model predictions, proceed ? (y/Y) [default=yes] ")
        if(!cok %in% c("y","Y","yes","")) return()
        }
        if(boolpred) {
          x<-saemix.predict(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.sim)>0) {
      if(length(x["sim.data"]["N"])==0 || x["sim.data"]["nsim"]==0) {
        boolpred<-TRUE
        if(interactive) {
        cok<-readline(prompt="Simulations will be performed. This might take a while, proceed ? (y/Y) [default=yes] ")
        if(!cok %in% c("y","Y","yes","")) return()
        } else {
        	cat("Performing simulations under the model.\n")
        }
        if(boolpred) {
          x<-simul.saemix(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.res)>0) {
      if(length(x["results"]["npde"])==0) {
        boolpred<-TRUE
        if(interactive) {
        cok<-readline(prompt="Simulations will be performed to obtain residuals, VPC and npde. This might take a while, proceed ? (y/Y) [default=yes] ")
        if(!cok %in% c("y","Y","yes","")) return()
        }
        if(boolpred) {
          x<-compute.sres(x)
          assign(namObj,x,envir=parent.frame())
        }
      }
    }
    if(length(id.map)>0) {
      if(length(x["results"]["map.eta"])==0) {
        cat("Computing ETA estimates and adding them to fitted object.\n")
	x<-compute.eta.map(x)
        assign(namObj,x,envir=parent.frame())
      }
    }
    for(ipl in plot.type) {
      switch (EXPR=ipl,
    "data"={
       cat("Plotting the data\n")
       saemix.plot.data(x,...)
    },
    "convergence"={
       cat("Plotting convergence plots\n")
       saemix.plot.convergence(x,...)
    },
    "likelihood"={  
       cat("Plotting the likelihood\n")
       saemix.plot.llis(x,...)
    },
    "observations.vs.predictions"={
       cat("Plotting observations versus predictions\n")      
       saemix.plot.obsvspred(x,...)
    },
    "individual.fit"={
      cat("Plotting individual fits\n")
      saemix.plot.fits(x,...)
    },
    "population.fit"={
      cat("Plotting fits obtained with population predictions\n")
      saemix.plot.fits(x,level=0,...)
    },
    "both.fit"={
      cat("Plotting the fits overlaying individual and population predictions\n")
      saemix.plot.fits(x,level=c(0,1),...)
    },
    "residuals.scatter"={
      cat("Plotting scatterplots of residuals\n")
      saemix.plot.scatterresiduals(x,...)
    },
    "residuals.distribution"={
      cat("Plotting the distribution of residuals\n")
      saemix.plot.distribresiduals(x,...)
    },
    "random.effects"={
      saemix.plot.randeff(x,...)
    },
    "correlations"={
      if(length(x@model@indx.omega>1)) saemix.plot.correlations(x,...)
    },
    "parameters.vs.covariates"={
      if(length(x@data@name.covariates)==0) {
        cat("No covariates in the dataset\n")
        return()
      } else saemix.plot.parcov(x,...)
    },
    "randeff.vs.covariates"={
      if(length(x@data@name.covariates)==0) {
        cat("No covariates in the dataset\n")
        return()
      } else saemix.plot.randeffcov(x,...)
    },
    "marginal.distribution"={
      saemix.plot.distpsi(x,...)
    },
    "vpc"={
      cat("Plotting VPC\n")
      saemix.plot.vpc(x,...)
    },
    "npde"={
      cat("Plotting npde\n")
      saemix.plot.npde(x,...)
    },
    cat("Plot ",ipl," not implemented yet\n")
     )
   }
  }
)

####################################################################################
####			Likelihood and tests		####
####################################################################################

##' Extract likelihood from a saemixObject resulting from a call to saemix
##'
##' The likelihood in saemix can be computed by one of three methods: linearisation (linearisation of the model), importance sampling (stochastic integration) and gaussian quadrature (numerical integration). The linearised likelihood is obtained as a byproduct of the computation of the Fisher Information Matrix (argument FIM=TRUE in the options given to the saemix function).
##' If no method argument is given, this function will attempt to extract the likelihood computed by importance sampling (method="is"), unless the object contains the likelihood computed by linearisation, in which case the function will extract this component instead.
##' If the requested likelihood is not present in the object, it will be computed and aded to the object before returning.
##'
##' @aliases logLik
##'
##' @param object name of an SaemixObject object
##' @param ... additional arguments; for SaemixObject objects, an argument 'method' can be given as a character string, one of c("is","lin","gq"), to select one of the available approximations to the log-likelihood (is: Importance Sampling; lin: linearisation and gq: Gaussian Quadrature). See documentation for details
##' 
##' @references Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.
##' 
##' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
#' @author Emmanuelle Comets \email{emmanuelle.comets@@inserm.fr}
#' @author Audrey Lavenu
#' @author Marc Lavielle.
#' @seealso \code{\link{AIC}},\code{\link{BIC}}, \code{\link{saemixControl}}, \code{\link{saemix}}
#' @exportMethod logLik

# Extract likelihood and number of parameters
setMethod("logLik",signature="SaemixObject",
  function(object,...) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("method",names(args1))
    if(!is.na(i1)) {
    	str1<-as.character(args1[[i1]])
    	if(str1 %in% c("is","lin","gq")) method<-str1
    } else {
      if(length(object@results@ll.is)!=0 | length(object@results@ll.lin)==0) method<-"is" else method<-"lin"
    }
# Compute the requested LL if not present
   namObj<-deparse(substitute(object))
    if(method=="is" & length(object@results@ll.is)==0) {
        object<-llis.saemix(object)
        assign(namObj,object,envir=parent.frame())
    }
    if(method=="gq" & length(object@results@ll.gq)==0) {
        object<-llgq.saemix(object)
        assign(namObj,object,envir=parent.frame())
    }
    if(method=="lin" & length(object@results@ll.lin)==0) {
        object<-fim.saemix(object)
        assign(namObj,object,envir=parent.frame())
    }
# OR: return if desired LL has not been computed
#     if(method=="gq" & length(object@results@ll.gq)==0) {
#     	cat("The log-likelihood by Gaussian Quadrature has not yet been computed.\n")
#     	invisible(NULL)
#     }
#     if(method=="is" & length(object@results@ll.is)==0) {
#     	cat("The log-likelihood by Importance Sampling has not yet been computed.\n")
#     	invisible(NULL)
#     }
#     if(method=="lin" & length(object@results@ll.lin)==0) {
#     	cat("The log-likelihood by linearisation has not yet been computed.\n")
#     	invisible(NULL)
#     }
    val<-switch(method,is=object@results@ll.is,lin=object@results@ll.lin, gq=object@results@ll.gq)
    res<-paste(format(val,digits=2,nsmall=2)," (df=",object@results@npar.est,")", sep="")
    res<-matrix(res,dimnames=list(paste("LL by",method)," "))
    print(res)
    attr(val, "nall") <- object@data@N
    attr(val, "nobs") <- object@data@ntot.obs
    attr(val, "df") <- object@results@npar.est
    class(val) <- "logLik"
    invisible(val)
  }
)
# Extract AIC
setMethod("AIC",signature="SaemixObject",
          function(object,...) {
            args1<-match.call(expand.dots=TRUE)
            i1<-match("method",names(args1))
            if(!is.na(i1)) {
              str1<-as.character(args1[[i1]])
              if(str1 %in% c("is","lin","gq")) method<-str1
            } else {
              if(length(object@results@ll.is)!=0 | length(object@results@ll.lin)==0) method<-"is" else method<-"lin"
            }
            # Compute the requested LL if not present
            namObj<-deparse(substitute(object))
            if(method=="is" & length(object@results@ll.is)==0) {
              object<-llis.saemix(object)
              assign(namObj,object,envir=parent.frame())
            }
            if(method=="gq" & length(object@results@ll.gq)==0) {
              object<-llgq.saemix(object)
              assign(namObj,object,envir=parent.frame())
            }
            if(method=="lin" & length(object@results@ll.lin)==0) {
              object<-fim.saemix(object)
              assign(namObj,object,envir=parent.frame())
            }
            val<-switch(method,is=object@results@aic.is,lin=object@results@aic.lin, gq=object@results@aic.gq)
            val
          }
)

# Extract likelihood and number of parameters
setMethod("BIC","SaemixObject",
          function(object,...) {
            args1<-match.call(expand.dots=TRUE)
            i1<-match("method",names(args1))
            if(!is.na(i1)) {
              str1<-as.character(args1[[i1]])
              if(str1 %in% c("is","lin","gq")) method<-str1
            } else {
              if(length(object@results@ll.is)!=0 | length(object@results@ll.lin)==0) method<-"is" else method<-"lin"
            }
            # Compute the requested LL if not present
            namObj<-deparse(substitute(object))
            if(method=="is" & length(object@results@ll.is)==0) {
              object<-llis.saemix(object)
              assign(namObj,object,envir=parent.frame())
            }
            if(method=="gq" & length(object@results@ll.gq)==0) {
              object<-llgq.saemix(object)
              assign(namObj,object,envir=parent.frame())
            }
            if(method=="lin" & length(object@results@ll.lin)==0) {
              object<-fim.saemix(object)
              assign(namObj,object,envir=parent.frame())
            }
            val<-switch(method,is=object@results@bic.is,lin=object@results@bic.lin, gq=object@results@bic.gq)
            val
          }
)

# Log-likelihood ratio test, given 2 models
# ECO TODO

####################################################################################
####			saemixObject class - accesseurs parametres		####
####################################################################################

# Extract individual parameter estimates (psi_i)
setMethod("psi","SaemixObject",
  function(object,indiv.par) {
    if(missing(indiv.par)) indiv.par<-"map"
    namObj<-deparse(substitute(object))
    if(indiv.par=="map") { # mode
      if(length(object@results@map.psi)==0) {
        object<-map.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      psi<-object@results@map.psi[,-c(1)]
    } else { # conditional means
      if(length(object@results@cond.mean.psi)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      psi<-object@results@cond.mean.psi
    }
    return(psi)
  }
)

# Extract individual parameter estimates on non-transformed scale (phi_i)
setMethod("phi","SaemixObject",
  function(object,indiv.par) {
    if(missing(indiv.par)) indiv.par<-"map"
    namObj<-deparse(substitute(object))
    if(indiv.par=="map") { # mode
      if(length(object@results@map.phi)==0) {
        object<-map.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      phi<-object@results@map.phi[,-c(1)]
    } else { # conditional means
      if(length(object@results@cond.mean.phi)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      phi<-object@results@cond.mean.phi
    }
    return(phi)
  }
)

# Extract individual estimates of random effects (eta_i)
setMethod("eta","SaemixObject",
  function(object,indiv.par) {
    if(missing(indiv.par)) indiv.par<-"map"
    namObj<-deparse(substitute(object))
#    cat("Nom objet",namObj,"\n")
    if(indiv.par=="map") { # mode
      if(length(object@results@map.eta)==0) {
        object<-compute.eta.map(object)
        assign(namObj,object,envir=parent.frame())
      }
      eta<-object@results@map.eta[,-c(1)]
    } else { # conditional means
      if(length(object@results@cond.mean.eta)==0) {
        object<-conddist.saemix(object)
        assign(namObj,object,envir=parent.frame())
      }
      eta<-object@results@cond.mean.eta
    }
    return(eta)
  }
)

# Extract coefficients
setMethod("coef","SaemixObject",
  function(object) {
#    if(missing(level)) level<-1
#    if(missing(indiv.par)) indiv.par<-"map"
    pfix<-object@results@fixed.effects[object@results@indx.fix]
    names(pfix)<-object@results@name.fixed[object@results@indx.fix]
#c(object@results@fixed.effects,object@name.res[object@indx.res])
#    names(pfix)<-c(object@results@name.fixed,object@name.res[object@indx.res])    
    pop.phi<-object@results@mean.phi
    pop.psi<-transphi(pop.phi,object@model@transform.par)
    ind.psi<-list(map=object@results@map.psi[,-c(1)], cond=object@results@cond.mean.psi)
    ind.phi<-list(map=object@results@map.phi[,-c(1)], cond=object@results@cond.mean.phi)
    eta<-list(map=object@results@map.eta[,-c(1)],cond=object@results@cond.mean.eta)
    colnames(pop.phi)<-colnames(pop.psi)<-colnames(ind.phi$map)<-names(pfix)
    colnames(ind.psi$map)<-colnames(ind.phi$cond)<-colnames(ind.psi$cond)<-names(pfix)
    coef<-list(fixed=pfix,population=list(psi=pop.psi,phi=pop.phi), individual=list(psi=ind.psi,phi=ind.phi,eta=eta))
    return(coef)
  }
)

# Extract residuals and fitted
setMethod("resid","SaemixObject",
          function (object, type = c("ires", "wres", "npde", "pd", "iwres", "icwres"), ...) 
          {
            type <- match.arg(type)
            resid(object@results, type=type, ...)
          }
)

setMethod("fitted","SaemixObject",
          function (object, type = c("ipred", "ypred", "ppred", "icpred"), ...) 
          {
            type <- match.arg(type)
            fitted(object@results, type=type, ...)
          }
)

####################################################################################
