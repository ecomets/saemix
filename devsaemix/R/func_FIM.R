###########################  Fisher Information Matrix & LL by linearisation 	#############################



#' Computes the Fisher Information Matrix by linearisation
#' 
#' Estimate by linearisation the Fisher Information Matrix and the standard
#' error of the estimated parameters.
#' 
#' The inverse of the Fisher Information Matrix provides an estimate of the
#' variance of the estimated parameters theta. This matrix cannot be computed
#' in closed-form for nonlinear mixed-effect models; instead, an approximation
#' is obtained as the Fisher Information Matrix of the Gaussian model deduced
#' from the nonlinear mixed effects model after linearisation of the function f
#' around the conditional expectation of the individual Gaussian parameters.
#' This matrix is a block matrix (no correlations between the estimated fixed
#' effects and the estimated variances).
#' 
#' @param saemixObject an object returned by the \code{\link{saemix}} function
#' @return The function returns an updated version of the object saemix.fit in
#' which the following elements have been added: \describe{
#' \item{se.fixed:}{standard error of fixed effects, obtained as part of the
#' diagonal of the inverse of the Fisher Information Matrix (only when
#' fim.saemix has been run, or when the saemix.options$algorithms[2] is 1)}
#' \item{se.omega:}{standard error of the variance of random effects, obtained
#' as part of the diagonal of the inverse of the Fisher Information Matrix
#' (only when fim.saemix has been run, or when the saemix.options$algorithms[2]
#' is 1)} \item{se.res:}{standard error of the parameters of the residual error
#' model, obtained as part of the diagonal of the inverse of the Fisher
#' Information Matrix (only when fim.saemix has been run, or when the
#' saemix.options$algorithms[2] is 1)} \item{fim:}{Fisher Information Matrix}
#' \item{ll.lin:}{ likelihood calculated by linearisation} }
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}}
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
#' # Running the main algorithm to estimate the population parameters
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), 
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="constant")
#' 
#' saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # Estimating the Fisher Information Matrix using the result of saemix 
#' # & returning the result in the same object
#' # fim.saemix(saemix.fit)
#' 
#' 
#' @export fim.saemix
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
