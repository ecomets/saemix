#######################	Conditional means - estimates of PSI_i ########################
conddist.saemix<-function(saemixObject,nsamp=1,max.iter=NULL,...) {
# Estimate conditional means and estimates for the individual parameters PSI_i using the MCMC algorithm
# nsamp= number of MCMC samples
# kmax= max nb of iterations
# returns an array 
  N<-saemixObject["data"]["N"]
  nb.parameters<-saemixObject["model"]["nb.parameters"]
  if(is.null(max.iter)) kmax<-sum(saemixObject["options"]$nbiter.saemix)*2 else kmax<-max.iter
# using several Markov chains
  chdat<-new(Class="SaemixRepData",data=saemixObject["data"], nb.chains=nsamp)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,saemixObject["data"]["name.predictors"],drop=FALSE]
  io<-matrix(data=0,nrow=N,ncol=max(saemixObject["data"]["nind.obs"]))
  for(i in 1:N)
    io[i,1:saemixObject["data"]["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),nsamp))
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])

  ind.eta<-saemixObject["model"]["indx.omega"]
  nb.etas<-length(ind.eta)
  omega.eta<-saemixObject["results"]["omega"][ind.eta,ind.eta]
  diag.omega<-mydiag(saemixObject["results"]["omega"])
  domega2<-do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))* saemixObject["options"]$rw.ini), nb.etas))
  VK<-rep(c(1:nb.etas),2)
  phi<-array(data=0,dim=c(N,nb.parameters, nsamp))

  domega<-cutoff(mydiag(saemixObject["results"]["omega"][ind.eta,ind.eta]), .Machine$double.eps)
  omega.eta<-saemixObject["results"]["omega"][ind.eta,ind.eta]
  omega.eta<-omega.eta-mydiag(mydiag(saemixObject["results"]["omega"][ind.eta, ind.eta]))+mydiag(domega)
  chol.omega<-chol(omega.eta)
  pres<-saemixObject["results"]["respar"]
  
# Preparing plots
  if(saemixObject["options"]$displayProgress) {
    plot.opt<-saemixObject["prefs"]
    plot.opt$xlab<-"Iteration"
    plot.opt<-replace.plot.options(plot.opt,...)
    change.ylab<-FALSE
    if(plot.opt$ylab!=saemixObject["prefs"]$ylab & length(plot.opt$which.par)==1) change.ylab<-TRUE
    change.main<-FALSE
    if(plot.opt$main!=saemixObject["prefs"]$main & length(plot.opt$which.par)==1) change.main<-TRUE
    np<-nb.etas
    if(length(plot.opt$mfrow)==0) {
      n1<-round(sqrt(np))
      n2<-ceiling(np/n1)
      if(n1>5 | n2>5) {
        n1<-3
        n2<-4
#      cat("Changing the plot layout\n")
      }
      plot.opt$mfrow<-c(n1,n2)
    }
  }
# Simulation MCMC
# initialisation a phiM=estimation des parametres individuels
  mean.phiM<-do.call(rbind,rep(list(saemixObject["results"]["mean.phi"]),nsamp))
  phiM<-do.call(rbind,rep(list(saemixObject["results"]["cond.mean.phi"]),nsamp))
  etaM<-phiM[,ind.eta]-mean.phiM[,ind.eta]  
  psiM<-transphi(phiM,saemixObject["model"]["transform.par"])
  fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
  if(saemixObject["model"]["error.model"]=="exponential")
     fpred<-log(cutoff(fpred))
  gpred<-error(fpred,pres)
  DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
  U.y<-colSums(DYF)
  phiMc<-phiM
  
  econd<-sdcond<-array(0,dim=c(nb.parameters,N*kmax,nsamp)) # parametres individuels
  ebar<-sdbar<-array(0,dim=c(nb.parameters,kmax, nsamp)) # moyenne des parametres individuels
  cat("Estimating the conditional mean and variance of the distribution of individual parameters\n")
  k<-1
  while(k<=kmax) { # Set a maximum nb of iterations

  for(u in 1:saemixObject["options"]$nbiter.mcmc[1]) { # 1er noyau
    etaMc<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
    phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
    psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
    fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
    if(saemixObject["model"]["error.model"]=="exponential")
      fpred<-log(cutoff(fpred))
    gpred<-error(fpred,pres)
    DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
    Uc.y<-colSums(DYF)
    deltau<-Uc.y-U.y
    ind<-which(deltau<(-1)*log(runif(NM)))
    etaM[ind,]<-etaMc[ind,]
    U.y[ind]<-Uc.y[ind]
  }

  U.eta<-0.5*rowSums(etaM*(etaM%*%solve(omega.eta)))

# Second stage

  if(saemixObject["options"]$nbiter.mcmc[2]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-1
    for (u in 1:saemixObject["options"]$nbiter.mcmc[2]) {
     for(vk2 in 1:nb.etas) {
       etaMc<-etaM
       etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2],nrow=1) # 2e noyau ? ou 1er noyau+permutation?
       phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
       psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
       fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
       if(saemixObject["model"]["error.model"]=="exponential")
         fpred<-log(cutoff(fpred))
       gpred<-error(fpred,pres)
       DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)**2+log(gpred)
       Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
       Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
       deltu<-Uc.y-U.y+Uc.eta-U.eta
       ind<-which(deltu<(-1)*log(runif(NM)))
       etaM[ind,]<-etaMc[ind,]
       U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
       U.eta[ind]<-Uc.eta[ind]
       nbc2[vk2]<-nbc2[vk2]+length(ind)
       nt2[vk2]<-nt2[vk2]+NM
     }
    }
    domega2[,nrs2]<-domega2[,nrs2]*(1+saemixObject["options"]$stepsize.rw* (nbc2/nt2-saemixObject["options"]$proba.mcmc))
  }

  if(saemixObject["options"]$nbiter.mcmc[3]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-k%%(nb.etas-1)+2
#    if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
    for (u in 1:saemixObject["options"]$nbiter.mcmc[3]) {
      if(nrs2<nb.etas) {
        vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
        nb.iter2<-nb.etas
      } else {
        vk<-0:(nb.etas-1)
        nb.iter2<-1
      }
      for(k2 in 1:nb.iter2) {
        vk2<-VK[k2+vk]
        etaMc<-etaM
        etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2])
        phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
        psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
        fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
        if(saemixObject["model"]["error.model"]=="exponential")
         fpred<-log(cutoff(fpred))
        gpred<-error(fpred,pres)
        DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)**2+log(gpred)
        Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
        Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
        deltu<-Uc.y-U.y+Uc.eta-U.eta
        ind<-which(deltu<(-log(runif(NM))))
        etaM[ind,]<-etaMc[ind,]
        U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
        U.eta[ind]<-Uc.eta[ind]
        nbc2[vk2]<-nbc2[vk2]+length(ind)
        nt2[vk2]<-nt2[vk2]+NM
      }
    }
    domega2[,nrs2]<-domega2[,nrs2]*(1+saemixObject["options"]$stepsize.rw* (nbc2/nt2-saemixObject["options"]$proba.mcmc))
  }
  
  phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
    if(k==1) {
      eik<-array(t(phiM),dim=c(nb.parameters,N,nsamp))
      varik<-0*eik
      ebar[,k,]<-apply(eik,c(1,3),mean)
      ij<-(k-1)*N
      econd[,(ij+1):(ij+N),]<-eik
      sdcond[,(ij+1):(ij+N),]<-varik
    } else {
      eik1<-eik
      eik<-eik*(k-1)/k+array(t(phiM),dim=c(nb.parameters,N,nsamp))/k
      varik<-varik*(k-1)/k+(array(t(phiM),dim=c(nb.parameters,N,nsamp))**2)/k+ (eik1**2)*(k-1)/k-(eik**2)
      sdik<-sqrt(varik)
      ij<-(k-1)*N
      econd[,(ij+1):(ij+N),]<-eik
      sdcond[,(ij+1):(ij+N),]<-sdik
      ebar[,k,]<-apply(eik,c(1,3),mean)
      sdbar[,k,]<-apply(sdik,c(1,3),mean)
      if(k>=saemixObject["options"]$ipar.lmcmc) {
      ibeg<-max(1,k-saemixObject["options"]$ipar.lmcmc)
      ekmax<-(ebar[,k,]+abs(ebar[,k,])* saemixObject["options"]$ipar.rmcmc)
      ekmin<-(ebar[,k,]-abs(ebar[,k,])*saemixObject["options"]$ipar.rmcmc)
        if(nsamp==1) {
          vec<-apply(ebar[,ibeg:k,],1,max)
          vec2<-apply(ebar[,ibeg:k,],1,min)
          iek<-sum(vec>ekmax)+sum(vec2<ekmin)
          vec<-apply(sdbar[,ibeg:k,],1,max)
          vec2<-apply(sdbar[,ibeg:k,],1,min)
          sdek<-sum(vec>(sdbar[,k,]+abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))+ sum(vec2<(sdbar[,k,]-abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))
        } else {
        vec<-apply(ebar[,ibeg:k,],c(1,3),max)
	vec2<-apply(ebar[,ibeg:k,],c(1,3),min)
        iek<-sum(vec>ekmax)+sum(vec2<ekmin)
	if(ibeg==1) ibeg<-2
        vec<-apply(sdbar[,ibeg:k,],c(1,3),max)
        vec2<-apply(sdbar[,ibeg:k,],c(1,3),min)
        sdek<-sum(vec>(sdbar[,k,]+abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))+ sum(vec2<(sdbar[,k,]-abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))
        }
        if(is.na(sdek)) sdek<-0
        if(iek==0 & sdek==0) break
      }
    }
    if((k%%50)==50) {
      cat(".")
      if(saemixObject["options"]$displayProgress) {
#        par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
        par(mfrow=plot.opt$mfrow)
        limy<-cbind(apply(cbind(apply(ebar[,1:k,],1,min),ekmin),1,min), apply(cbind(apply(ebar[,1:k,],1,max),ekmax),1,max))
        for(ipar in ind.eta) {
        laby<-saemixObject["model"]["name.modpar"][ipar]
	plot(1:k,ebar[ipar,1:k,1],type="n",xlab=plot.opt$xlab, ylab=laby,ylim=limy[ipar,],main=plot.opt$main)
        for(isamp in 1:nsamp) 
	  lines(1:k,ebar[ipar,1:k,isamp],col=plot.opt$col,lty=plot.opt$lty, lwd=plot.opt$lwd)
        abline(h=apply(ebar[,k,],1,mean), col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
        abline(h=apply(ekmax,1,mean)[ipar], col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
        abline(h=apply(ekmin,1,mean)[ipar], col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
        }
      }
    }
    k<-k+1
  }
  cat("\n")
  if(k>=kmax)
    cat("Computing the empirical conditional mean and variance: maximum number of iterations reached without meeting convergence criterion (max.iter=",kmax,")\n") else cat("Convergence achieved in",k,"iterations\n")
  eta.cond<-matrix(0,nrow=dim(etaM)[1],ncol=nb.parameters)
  eta.cond[,ind.eta]<-etaM
  eta.cond<-array(t(eta.cond),dim=c(nb.parameters,N,nsamp))
  eta.cond<-t(apply(eta.cond,c(1,2),mean))
  cond.shrinkage<-100*(1-apply(eta.cond,2, var)/mydiag(saemixObject["results"]["omega"]))
  names(cond.shrinkage)<-paste("Sh.",names(cond.shrinkage),".%",sep="")
  resh.eik<-resh.varik<-array(0,dim=c(nrow=N,ncol=nb.parameters,nsamp))
  psi.samp<-resh.eik
  for(isamp in 1:nsamp) {
    resh.eik[,,isamp]<-t(eik[,,isamp])
    resh.varik[,,isamp]<-t(varik[,,isamp])
    psi.samp[,,isamp]<-transphi(resh.eik[,,isamp],saemixObject["model"]["transform.par"])
  }

  saemixObject["results"]["cond.mean.phi"]<-apply(resh.eik,c(1,2),mean)
  saemixObject["results"]["cond.var.phi"]<-apply(resh.varik,c(1,2),mean)
  saemixObject["results"]["cond.mean.eta"]<-eta.cond
  saemixObject["results"]["cond.shrinkage"]<-cond.shrinkage
  saemixObject["results"]["phi.samp"]<-resh.eik
  saemixObject["results"]["phi.samp.var"]<-resh.varik
  saemixObject["results"]["psi.samp"]<-psi.samp

# ECO TODO: verifier ce qui se passe quand un parametre n'est pas estime => faut-il prendre la moyenne des valeurs ou grace a chol.omega qui vaut 0, ca marche ? Verifier d'ailleurs si on a chol.omega quand il manque un parametre (normalement oui)

  return(saemixObject)
}

###########################  Individual MAP estimates 	#############################

map.saemix<-function(saemixObject) {
# Compute the MAP estimates of the individual parameters PSI_i
  i1.omega2<-saemixObject["model"]["indx.omega"]
  iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
  id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
  yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  id.list<-unique(id)
  phi.map<-saemixObject["results"]["phi"]
  
  cat("Estimating the individual parameters, please wait a few moments...\n")
  for(i in 1:saemixObject["data"]["N"]) {
    cat(".")
    isuj<-id.list[i]
    xi<-xind[id==isuj,,drop=FALSE]
#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
    yi<-yobs[id==isuj]
    idi<-rep(1,length(yi))
    mean.phi1<-saemixObject["results"]["mean.phi"][i,i1.omega2]
    phii<-saemixObject["results"]["phi"][i,]
    phi1<-phii[i1.omega2]
    phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"])
    phi.map[i,i1.omega2]<-phi1.opti$par
  }
  cat("\n")
  map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
  map.psi<-data.frame(id=id.list,map.psi)
  map.phi<-data.frame(id=id.list,phi.map)
  colnames(map.psi)<-c(saemixObject["data"]["name.group"], saemixObject["model"]["name.modpar"])
  saemixObject["results"]["map.psi"]<-map.psi
  saemixObject["results"]["map.phi"]<-map.phi
  return(saemixObject)
}

compute.eta.map<-function(saemixObject) {
# Compute individual estimates of the MAP random effects from the MAP estimates of the parameters
# returns the parameters (psi), newly computed if needs be, the corresponding random effects, and the associated shrinkage
  if(length(saemixObject["results"]["map.psi"])) {
      saemixObject<-map.saemix(saemixObject)
  }
  psi<-saemixObject["results"]["map.psi"][,-c(1)]
  phi<-transpsi(as.matrix(psi),saemixObject["model"]["transform.par"])

# Computing COV again here (no need to include it in results)  
#  COV<-matrix(nrow=dim(saemix.model["Mcovariates"])[1],ncol=0)
  COV<-matrix(nrow=saemixObject["data"]["N"],ncol=0)
  for(j in 1:saemixObject["model"]["nb.parameters"]) {
    jcov<-which(saemixObject["model"]["betaest.model"][,j]==1)
    aj<-as.matrix(saemixObject["model"]["Mcovariates"][,jcov])
    COV<-cbind(COV,aj)
  }
  eta<-phi-COV%*%saemixObject["results"]["MCOV"] 
  shrinkage<-100*(1-apply(eta,2,var)/mydiag(saemixObject["results"]["omega"]))
  names(shrinkage)<-paste("Sh.",names(shrinkage),".%",sep="")
  colnames(eta)<-paste("ETA(",colnames(eta),")",sep="")
  eta<-cbind(id=saemixObject["results"]["map.psi"][,1],eta)
  
  saemixObject["results"]["map.eta"]<-eta
  saemixObject["results"]["map.shrinkage"]<-shrinkage
  
  return(saemixObject)
}

###########################  Fisher Information Matrix 	#############################

fim.saemix<-function(saemixObject) {
# Estimate the Fisher Information Matrix and the s.e. of the estimated parameters  

	saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]

  covariance.model<-0*saemix.model["covariance.model"]
  diag(covariance.model)<-mydiag(saemix.model["covariance.model"])
  omega<-0*saemix.res["omega"]
  diag(omega)<-mydiag(saemix.res["omega"])
  hat.phi<-saemix.res["cond.mean.phi"]
  nphi<-dim(hat.phi)[2]
  dphi<-cutoff(abs(colMeans(hat.phi))*1e-4,1e-10)
  coefphi<-c(0,-1,1)

  F<-array(data=0,dim=c(saemix.data["ntot.obs"],nphi,length(coefphi)))
  gs<-matrix(0,saemix.data["ntot.obs"],4)

  for (l in 1:length(coefphi)) {
    for (j in 1:nphi) {
        phi<-hat.phi
        phi[,j]<-phi[,j]+coefphi[l]*dphi[j]
        psi<-transphi(phi,saemix.model["transform.par"])
        f <- saemix.model["model"](psi, saemix.data["data"][,"index"],xind)
      	if(saemix.model["error.model"]=='exponential') 
		f<-log(cutoff(f))        
        F[,j,l]<-f
    	}
  }

  ind.covariates<-which(saemix.model["betaest.model"]>0)
  f0<-F[,1,1]
  g0<-cutoff(saemix.res["respar"][1]+saemix.res["respar"][2]*abs(f0))
  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 #gradient of f 
  z<-matrix(0,saemix.data["ntot.obs"],1)
  j2<-0

	for (i in 1:saemix.data["N"]) {
    j1<-j2+1
    j2<-j2+saemix.data["nind.obs"][i]
    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF[j1:j2,,drop=FALSE]%*%hat.phi[i,]
  }

# ECO ici modifie car role de covariate.estim pas clair
# covariate.estim=si un parametre (et ses covariables associees) sont estimees ou non
  covariate.estim<-matrix(rep(saemix.model["fixed.estim"], dim(saemix.model["betaest.model"])[1]),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))*saemix.model["betaest.model"]
  j<-which(saemix.model["betaest.model"]>0)
  ind.fixed.est<-(covariate.estim[j]>0)

# hw=waitbar(1,'Estimating the population parameters (SAEM). Wait...')

  ll.lin<- -0.5*saemix.data["ntot.obs"]*log(2*pi)
  Fmu<-0
  FO<-0
  j2<-0
	for (i in 1:saemix.data["N"]) {
#waitbar(i/N,hw)
    ni<-saemix.data["nind.obs"][i]
    j1<-j2+1
    j2<-j2+ni
    yi<-yobs[j1:j2]
    DFi<-DF[j1:j2,,drop=FALSE]
    f0i<-f0[j1:j2]
    g0i<-g0[j1:j2]
    zi<-z[j1:j2]
    Ai<-kronecker(diag(nphi),as.matrix(saemix.model["Mcovariates"][i,]))
    Ai<-Ai[,ind.covariates,drop=FALSE]
    DFAi<-DFi%*%Ai
    Gi<-DFi%*%omega%*%t(DFi) + mydiag(g0i^2,nrow=ni)  #variance of zi
  
    Gi<-round(Gi*1e10)/1e10
    VD<-try(eigen(Gi))
    if(class(VD)=="try-error") {
      cat("Unable to compute the FIM by linearisation.\n")
      return(saemixObject)
    }
    D<-Re(VD$values)
    V<-Re(VD$vectors)
    IGi<-V%*%mydiag(1/D,nrow=length(D))%*%t(V)
    Dzi<-zi-DFAi%*%saemix.res["betas"]
    
    if (sum(ind.fixed.est)>0) {
        DFAiest<-DFAi[,ind.fixed.est,drop=FALSE]
        Fmu<-Fmu-t(DFAiest)%*%IGi%*%DFAiest
    }
    
    ###########################################################
    OP<-NULL
    for (k in 1:nphi) {
        for (l in 1:nphi) {
            if (covariance.model[k,l]==1) {
                OPkl<-DFi[,k,drop=FALSE]%*%t(DFi[,l,drop=FALSE])
                OP<-cbind(OP,c(OPkl))
        	}
        }
    }
    if (sum(saemix.res["indx.res"]==1)>0) {
        SIi<-2*g0i
        dSIi<-mydiag(SIi)
        OP<-cbind(OP,c(dSIi))
    	}
    if (sum(saemix.res["indx.res"]==2)>0) {
        SIi<-2*f0i%*%g0i
        dSIi<-mydiag(SIi)
	  OP<-cbind(OP,c(dSIi))
    	}
    kl<-0
    FG<-matrix(0,dim(OP)[2],ni*ni)
    for (k in 1:ni)
	{
        for (l in 1:ni)
	{
            FGkl<- -IGi[,k]%*%t(IGi[l,])/2
            kl<-kl+1
            FG[,kl]<-t(OP)%*%c(FGkl)
        }
    }
    FO<-FO+FG%*%OP
    
    ###############################################################
    ll.lin <- ll.lin - 0.5*log(det(Gi)) - 0.5*t(Dzi)%*%IGi%*%Dzi 
  }
#partie precedente verifiee pas a pas avec matlab

  if(saemix.model["error.model"]=='exponential')
    ll.lin<-ll.lin-sum(yobs)

  if (sum(ind.fixed.est)>0) {
    Mparam<-matrix(0,dim(saemix.model["betaest.model"])[1], dim(saemix.model["betaest.model"])[2])
    Mparam[1,]<-saemix.model["transform.par"]
    Mtp<-Mparam[saemix.model["betaest.model"]>0]    
    Mtp<-Mtp[ind.fixed.est]
    dbetas <- dtransphi(saemix.res["betas"][ind.fixed.est],Mtp)
    Mupth<-mydiag(1/dbetas,nrow=length(dbetas))
    Fth<-t(Mupth)%*%Fmu%*%Mupth
    Cth<-try(solve(-Fth))
    if(class(Cth)=="try-error") {
      cat("Error computing the Fisher Information Matrix: singular system.\n")
      Cth<-NA*Fth
    }
  } else {
    Cth<-NULL
  }

  fim<-rbind(cbind(Fth,matrix(0,dim(Fth)[1],dim(FO)[2])), cbind(matrix(0,dim(FO)[1],dim(Fth)[2]),FO)) 

  sTHest<-sqrt(mydiag(Cth))
#sTH<-matrix(0,1,length(saemix.res["betas"]))
  sTH<-rep(0,length(saemix.res["betas"]))
  sTH[ind.fixed.est]<-sTHest
  se.fixed<-sTH

  CO<-try(solve(-FO))
    if(class(CO)=="try-error") {
      CO<-NA*FO
      cat("Error computing the Fisher Information Matrix: singular system.\n")
  }
  sO<-sqrt(mydiag(CO))
  nb.omega2<-length(saemix.model["indx.omega"])
  se.omega<-matrix(0,nphi,1)
  se.omega[saemix.model["indx.omega"]]<-sO[1:nb.omega2]
  se.res<-matrix(0,2,1)
  se.res[saemix.res["indx.res"]]<-sO[(nb.omega2+1):length(sO)]    

# ECO TODO : pourquoi negatif ??? FIM = -fim calculee ici ?
  saemix.res["se.fixed"]<-se.fixed
  saemix.res["se.omega"]<-c(se.omega)
  saemix.res["se.respar"]<-c(se.res)
  saemix.res["ll.lin"]<-c(ll.lin )
  saemix.res["fim"]<-fim
  saemix.res["aic.lin"]<-(-2)*saemix.res["ll.lin"]+ 2*saemix.res["npar.est"]
  saemix.res["bic.lin"]<-(-2)*saemix.res["ll.lin"]+ log(saemix.data["N"])*saemix.res["npar.est"]

##################################
#delete(hw)
  saemixObject["results"]<-saemix.res
  return(saemixObject)
}

#######################	Model simulations and residuals ########################

simul.saemix<-function(saemixObject,nsim=saemixObject["options"]$nb.sim, predictions=TRUE,res.var=TRUE,uncertainty=FALSE) {
# Simulate individual parameters from the population distribution
# predictions: if TRUE, use the parameters to predict observations
# res.var: if TRUE, add residual error to the predictions to obtain simulated data
# uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]

  N<-saemix.data["N"]
  ind.eta<-saemix.model["indx.omega"]
  nb.etas<-length(ind.eta)
  NM <- N*nsim  
  domega<-cutoff(mydiag(saemix.res["omega"][ind.eta, ind.eta]),.Machine$double.eps)
  omega.eta<-saemix.res["omega"][ind.eta,ind.eta]
  omega.eta<-omega.eta-mydiag(mydiag(saemix.res["omega"][ind.eta,ind.eta]))+mydiag(domega)
  chol.omega<-chol(omega.eta)

  phiM<-mean.phiM<-do.call(rbind,rep(list(saemix.res["mean.phi"]),nsim))
  etaM<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
  phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
  psiM<-transphi(phiM,saemix.model["transform.par"])

  if(predictions) {
    index<-rep(1:N,times=saemix.data["nind.obs"])
    IdM<-kronecker(c(0:(nsim-1)),rep(N,saemix.data["ntot.obs"]))+rep(index,nsim)
    XM<-do.call(rbind,rep(list(xind),nsim))
    pres<-saemix.res["respar"]
    sim.pred<-sim.data<-NULL
      fpred<-saemix.model["model"](psiM, IdM, XM)
      sim.pred<-fpred
      if(res.var) {
        if(saemix.model["error.model"]=="exponential")
          fpred<-log(cutoff(fpred))
        gpred<-error(fpred,pres)
        eps<-rnorm(length(fpred))
        sim.data<-fpred+gpred*eps
      }
  } else {
    sim.pred<-sim.data<-IdM<-c()
  }
  sim.psi<-data.frame(id=rep(unique(saemix.data["data"][, saemix.data["name.group"]]),nsim),psiM)
  colnames(sim.psi)<-c(saemix.data["name.group"],saemix.model["name.modpar"])
  datasim<-data.frame(idsim=rep(index,nsim),irep=rep(1:nsim, each=saemix.data["ntot.obs"]),ypred=sim.pred,ysim=sim.data)
  ysim<-new(Class="SaemixSimData",saemix.data,datasim)
  ysim["sim.psi"]<-sim.psi
  saemixObject["sim.data"]<-ysim
  
  return(saemixObject)
}

compute.sres<-function(saemixObject) {
# Compute standardised residuals (WRES, npd and npde) using simulations
# saemix.options$nb.sim simulated datasets used to compute npd, npde, and VPC
# saemix.options$nb.simpred simulated datasets used to compute ypred and WRES ? for the moment saemix.options$nb.sim used for both
  nsim<-saemixObject["options"]$nb.sim
  if(length(saemixObject["sim.data"]["N"])==0 || saemixObject["sim.data"]["nsim"]!=nsim) {
	  cat("Simulating data using nsim =",nsim,"simulated datasets\n")
	  saemixObject<-simul.saemix(saemixObject,nsim)
  }
# ECO TODO: maybe here be more clever and use simulations if available (adding some if not enough, truncating if too much ?)  
  ysim<-saemixObject["sim.data"]["datasim"]$ysim
  idsim<-saemixObject["sim.data"]["datasim"]$idsim
  idy<-saemixObject["data"]["data"][,"index"]
  yobsall<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  ypredall<-pd<-npde<-wres<-c()
#  pde<-c()
  cat("Computing WRES and npde ")
  for(isuj in 1:saemixObject["data"]["N"]) {
    if(isuj%%10==1) cat(".")
    ysimi<-matrix(ysim[idsim==isuj],ncol=nsim)
#    ysimi.pred<-ysimi[,1:saemixObject["options"]$nb.simpred]
    yobs<-yobsall[idy==isuj]
    tcomp<-apply(cbind(ysimi,yobs),2,"<",yobs)
    if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
    pdsuj<-rowMeans(tcomp)
#      pdsuj[pdsuj==0]<-1/nsim
#      pdsuj[pdsuj==1]<-1-1/nsim
    pdsuj<-pdsuj+0.5/(nsim+1)
# pdsuj=0 pour yobs<min(tabobs$yobs)
# pdsuj=1 : jamais
# dc dist va de 0 ? 1-1/(nsim+1)
# dc dist+0.5/(nsim+1) va de 0.5/(nsim+1) ? 1-0.5/(nsim+1)
# est-ce que ?a recentre ma distribution ? ECO TODO CHECK
    pd<-c(pd,pdsuj)
    ypred<-rowMeans(ysimi)
    ypredall<-c(ypredall,ypred)
    xerr<-0
    if(length(yobs)==1) {
      npde<-c(npde,qnorm(pdsuj))
      wres<-c(wres,(yobs-ypred)/sd(t(ysimi)))
    } else {
# Decorrelation
    vi<-cov(t(ysimi))
    xmat<-try(chol(vi))
    if(is.numeric(xmat)) {
      sqvi<-try(solve(xmat))
      if(!is.numeric(sqvi)) 
        xerr<-2
      } else xerr<-1
    if(xerr==0) {
    #decorrelation of the simulations
      decsim<-t(sqvi)%*%(ysimi-ypred)
      decobs<-t(sqvi)%*%(yobs-ypred)
      wres<-c(wres,decobs)
#    ydsim<-c(decsim)
#    ydobs<-decobs
    #Computing the pde
      tcomp<-apply(cbind(decsim,decobs),2,"<",decobs)
      if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
      pdesuj<-rowMeans(tcomp)
      pdesuj<-pdesuj+0.5/(nsim+1)
#      pde<-c(pde,pdesuj)
      npde<-c(npde,qnorm(pdesuj))
    } else {
      npde<-c(npde,rep(NA,length(yobs)))
      wres<-c(wres,rep(NA,length(yobs)))
    }
    }
  }
  cat("\n")
  saemixObject["results"]["npde"]<-npde
  saemixObject["results"]["wres"]<-wres
  saemixObject["results"]["ypred"]<-ypredall # = [ E_i(f(theta_i)) ]
  saemixObject["results"]["pd"]<-pd
  if(length(saemixObject["results"]["predictions"])==0) saemixObject["results"]["predictions"]<-data.frame(ypred=ypredall,wres=wres,pd=pd,npde=npde) else {
  	saemixObject["results"]["predictions"]$ypred<-ypredall
  	saemixObject["results"]["predictions"]$wres<-wres
  	saemixObject["results"]["predictions"]$npde<-npde
  	saemixObject["results"]["predictions"]$pd<-pd
  }
#  return(list(ypred=ypredall,pd=pd,npd=npd,wres=wres,sim.data=ysim, sim.pred=x$sim.pred))
  return(saemixObject)
}

###########################	Computational fcts	#############################
# Redefining diag function, too many problems with the R version
mydiag <- function (x = 1, nrow, ncol) {
	if (is.matrix(x)) {
		if (nargs() > 1L) 
			stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
		if ((m <- min(dim(x))) == 0L) 
			return(vector(typeof(x), 0L))
		y <- c(x)[1L + 0L:(m - 1L) * (dim(x)[1L] + 1L)]
		nms <- dimnames(x)
		if (is.list(nms) && !any(sapply(nms, is.null)) && identical((nm <- nms[[1L]][seq_len(m)]), 
																																nms[[2L]][seq_len(m)])) 
			names(y) <- nm
		return(y)
	}
	if (is.array(x) && length(dim(x)) != 1L) 
		stop("'x' is an array, but not 1D.")
	if (missing(x)) 
		n <- nrow
	else n <- length(x)
	if (!missing(nrow)) 
		n <- nrow
	if (missing(ncol)) 
		ncol <- n
	p <- ncol
	y <- array(0, c(n, p))
	if ((m <- min(n, p)) > 0L) 
		y[1L + 0L:(m - 1L) * (n + 1L)] <- x
	y
}

cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}
cutoff.max<-function(x) max(x,.Machine$double.xmin)
cutoff.eps<-function(x) max(x,.Machine$double.eps)
cutoff.res<-function(x,ares,bres) max(ares+bres*abs(x),.Machine$double.xmin)

# Inverse of the normal cumulative distribution fct: using erfcinv from ?pnorm
norminv<-function(x,mu=0,sigma=1)  mu-sigma*qnorm(x,lower.tail=FALSE)

# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

error<-function(f,ab,type="combined") {
  g<-cutoff(ab[1]+ab[2]*abs(f))
  return(g)
}

ssq<-function(ab,y,f) { # Sum of squares
	g<-abs(ab[1]+ab[2]*f)
	e<-sum(((y-f)/g)**2+2*log(g))
	return(e)
}

transpsi<-function(psi,tr) {
  phi<-psi
#  if(is.null(dim(psi))) phi<-as.matrix(t(phi),nrow=1)
# ECO TODO: pourquoi ce test ??? Dans le cas ou psi est un vecteur ?
  i1<-which(tr==1) # log-normal
  phi[,i1]<-log(phi[,i1])
  i2<-which(tr==2) # probit
  phi[,i2]<-norminv(phi[,i2])
  i3<-which(tr==3) # logit
  phi[,i3]<-log(phi[,i3]/(1-phi[,i3]))
  if(is.null(dim(psi))) phi<-c(phi)
  return(phi)
}

transphi<-function(phi,tr) {
  psi<-phi
#  if(is.null(dim(psi))) psi<-as.matrix(t(psi),nrow=1)
  i1<-which(tr==1) # log-normal
  psi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-normcdf(psi[,i2])
  i3<-which(tr==3) # logit
  psi[,i3]<-1/(1+exp(-psi[,i3]))
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

derivphi<-function(phi,tr) {
# Fonction calculant la derivee de h pour tracer la distribution des parametres
  psi<-phi # identite
  i1<-which(tr==1) # log-normal
  psi[,i1]<-1/exp(phi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-1/(sqrt(2*pi))*exp(-(phi[,i2]**2)/2)
  i3<-which(tr==3) # logit
  psi[,i3]<-2+exp(phi[,i3])+exp(-phi[,i3])
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

dtransphi<-function(phi,tr) {
  psi<-phi
  if(is.null(dim(phi))) {
     dpsi<-as.matrix(t(rep(1,length(phi))))
     psi<-as.matrix(t(phi),nrow=1)
  } else 
    dpsi<-matrix(1,dim(phi)[1],dim(phi)[2])
  i1<-which(tr==1) # log-normal
  dpsi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  dpsi[,i2]<-1/dnorm(qnorm(dpsi[,i2]))   # derivee de la fonction probit, dqnorm <- function(p) 1/dnorm(qnorm(p))
  i3<-which(tr==3) # logit
  dpsi[,i3]<-1/(2+exp(-psi[,i3])+exp(psi[,i3]))
  if(is.null(dim(phi))) dpsi<-c(dpsi)
  return(dpsi)
}

compute.Uy<-function(b0,phiM,pres,args,Dargs,DYF) {
# Attention, DYF variable locale non modifiee en dehors
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nchains))
  psiM<-transphi(phiM,Dargs$transform.par)
  fpred<-Dargs$structural.model(psiM,Dargs$IdM,Dargs$XM)
  if(Dargs$error.model=="exponential")
     fpred<-log(cutoff(fpred))
  gpred<-error(fpred,pres)
  DYF[args$ind.ioM]<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
  U<-sum(DYF)
  return(U)
}

conditional.distribution<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model,pres,err) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
  if(err=="exponential")
    fi<-log(cutoff(fi))
  gi<-cutoff((pres[1]+pres[2]*abs(fi)))
  Uy<-sum(0.5*((yi-fi)/gi)**2+log(gi))
  dphi<-phi1-mphi
  Uphi<-0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}

trnd.mlx<-function(v,n,m) {
  r<-rnorm(n*m)*sqrt(v/2/gammarnd.mlx(v/2,n,m))
  return(r=matrix(r,nrow=n,ncol=m))
}

gammarnd.mlx<-function(a,n,m) {
  nm<-n*m
  y0 <- log(a)-1/sqrt(a)
  c <- a - exp(y0)
  b <- ceiling(nm*(1.7 + 0.6*(a<2)))
  y <- log(runif(b))*sign(runif(b)-0.5)/c + log(a)
  f <- a*y-exp(y) - (a*y0 - exp(y0))
  g <- c*(abs((y0-log(a))) - abs(y-log(a)))
  reject <- ((log(runif(b)) + g) > f)
  y<-y[!reject]
  if(length(y)>=nm) x<-exp(y[1:nm]) else 
    x<-c(exp(y),gammarnd.mlx(a,(nm-length(y)),1))
#  x<-matrix(x,nrow=n,ncol=m) # not useful ?
  return(x)
}

tpdf.mlx<-function(x,v) {
# TPDF_MLX  Probability density function for Student's T distribution

    term<-exp(lgamma((v + 1) / 2) - lgamma(v/2))
    return(term/(sqrt(v*pi)*(1+(x**2)/v)**((v+1)/2)))
}

###########################	Functions for npde	#############################

kurtosis<-function (x) 
{
#from Snedecor and Cochran, p 80
    x<-x[!is.na(x)]
    m4<-sum((x - mean(x))^4)
    m2<-sum((x - mean(x))^2)
    kurt<-m4*length(x)/(m2**2)-3
    return(kurtosis=kurt)
}
skewness<-function (x) 
{
#from Snedecor and Cochran, p 79
    x<-x[!is.na(x)]
    m3<-sum((x - mean(x))^3)
    m2<-sum((x - mean(x))^2)
    skew<-m3/(m2*sqrt(m2/length(x)))
    return(skewness=skew)
}
   
testnpde<-function(npde) 
{
    cat("---------------------------------------------\n")
    cat("Distribution of npde:\n")
    sev<-var(npde)*sqrt(2/(length(npde)-1))
    sem<-sd(npde)/sqrt(length(npde))
    cat("           mean=",format(mean(npde),digits=4),"  (SE=",format(sem,digits=2),")\n")
    cat("       variance=",format(var(npde),digits=4),"  (SE=",format(sev,digits=2),")\n")
    cat("       skewness=",format(skewness(npde),digits=4),"\n")
    cat("       kurtosis=",format(kurtosis(npde),digits=4),"\n")
    cat("---------------------------------------------\n\n")
    myres<-rep(0,4)
    y<-wilcox.test(npde)
    myres[1]<-y$p.val
    y<-shapiro.test(npde)
    myres[3]<-y$p.val

    # test de variance pour 1 ?chantillon
    # chi=s2*(n-1)/sigma0 et test de H0={s=sigma0} vs chi2 ? n-1 df
    semp<-sd(npde)
    n1<-length(npde)
    chi<-(semp**2)*(n1-1)
    y<-2*min(pchisq(chi,n1-1),1-pchisq(chi,n1-1))
    myres[2]<-y
    xcal<-3*min(myres[1:3])
    myres[4]<-min(1,xcal)
    names(myres)<-c("  Wilcoxon signed rank test ","  Fisher variance test      ",
    "  SW test of normality      ","Global adjusted p-value     ")
    cat("Statistical tests\n")
    for(i in 1:4) {
      cat(names(myres)[i],": ")
      #if (myres[i]<1) 
      cat(format(myres[i],digits=3)) 
      #else cat(myres[i])
      if(as.numeric(myres[i])<0.1 & as.numeric(myres[i])>=0.05) cat(" .")
      if(as.numeric(myres[i])<0.05) cat(" *")
      if(as.numeric(myres[i])<0.01) cat("*")
      if(as.numeric(myres[i])<0.001) cat("*")
      cat("\n")
    }
      cat("---\n")
      cat("Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 \n")
    cat("---------------------------------------------\n")
    return(myres)
}
