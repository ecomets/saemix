#####################################################################################
# R code for paper:
# Parameter estimation in nonlinear mixed effect models using saemix, an R implementation of the SAEM algorithm
# Authors: Emmanuelle Comets, Audrey Lavenu, Marc Lavielle
# Submitted to JSS: February 2014
# Revision received: April 2015
# Revised:

# Setting up working directory (other links relative)
# Required: download the files from XXX
# unzip these files to subdirectory

work.dir<-"/home/eco/work/monolix/rversion/current/saemix"
# work.dir<-getwd()
setwd(work.dir)

# Toggle using library or using development code
develVersion<-TRUE

# Toggle estimation for datasets from Plan et al. 2012, scenarios R3A and S3A
estimPlan12<-FALSE

# Toggle save feature
# saving graphs and results feature disabled in all the runs
save.results<-TRUE
save.results<-FALSE
save.dir<-"PaperExamples"
if(save.results & !file_test("-d",save.dir)) dir.create(save.dir)

#####################################################################################
if(develVersion) {
  source("devsaemix/R/global.R")
  
  # SaemixData
  source("devsaemix/R/SaemixData.R")
  # SaemixRes
  source("devsaemix/R/SaemixRes.R")
  # SaemixModel
  source("devsaemix/R/SaemixModel.R")
  # SaemixObject
  source("devsaemix/R/SaemixObject.R")
  
  # Functions
  source("devsaemix/R/func_main.R")
  source("devsaemix/R/func_aux.R")
  source("devsaemix/R/main_initialiseMainAlgo.R")
  source("devsaemix/R/main_estep.R")
  source("devsaemix/R/main_mstep.R")
  source("devsaemix/R/func_distcond.R")
  source("devsaemix/R/func_FIM.R")
  source("devsaemix/R/compute_LL.R")
  source("devsaemix/R/func_simulations.R")
  source("devsaemix/R/func_plots.R")
} else {
  # Pre-requisite: install saemix library version 2.0
  library(saemix)
}

#####################################################################################
#### Example 3.1 - Orange Trees

data(Orange)
head(Orange)

orange<-saemixData(name.data=Orange,name.group="Tree",name.predictor="age", name.response="circumference", units=list(x="Days",y="mm"))

plot(orange)

logistic.model<-function(psi,id,xidep) { 
	age<-xidep[,1]  
	lambda1<-psi[id,1]
	lambda2<-psi[id,2]
	lambda3<-psi[id,3]
	resp<-lambda1/(1+exp(-(age-lambda2)/lambda3))
	return(resp) 
}

# IIV on lambda1 only

orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(200,800,400),ncol=3,byrow=TRUE, dimnames=list(NULL,c("lambda1","lambda2", "lambda3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3))

opt<-list(seed=94352514,save=save.results,save.graphs=save.results,directory=file.path(save.dir,"paperOrange"))

orange.fit<-saemix(orange.model,orange,opt)

#####################################################################################
#### Example 3.2 - Theophylline PK

if(develVersion) {
  theo.saemix<-read.table(file.path(work.dir,"keepsaemix","data","theo.saemix.tab"),header=T)
} else data(theo.saemix)

theo.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
   name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"), 
   name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), 
   name.X="Time")

model1cpt<-function(psi,id,xidep) { 
	  dose<-xidep[,1]
	  tim<-xidep[,2]  
	  ka<-psi[id,1]
	  V<-psi[id,2]
	  CL<-psi[id,3]
	  k<-CL/V
	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
	  return(ypred)
}

###########################
# Model with a covariate effect (effect of Weight on CL)

theo.model<-saemixModel(model=model1cpt,
   description="One-compartment model with first-order absorption", 
   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, 
   c("ka","V","CL"))),transform.par=c(1,1,1), 
   covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE))

opt<-list(save=save.results,save.graphs=save.results,directory=file.path(save.dir,"paperTheoCov"))

theo.fit<-saemix(theo.model,theo.data,opt)
theo.fit<-llgq.saemix(theo.fit)

###########################
# Model without covariate

theo.model.base<-saemixModel(model=model1cpt,
   description="One-compartment model with first-order absorption", 
   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, 
   c("ka","V","CL"))),transform.par=c(1,1,1))

opt<-list(save=save.results,save.graphs=save.results,directory=file.path(save.dir,"paperTheoNoCov"))

theo.base<-saemix(theo.model.base,theo.data,opt)
theo.base<-llgq.saemix(theo.base)

ll1<-theo.base["results"]["ll.lin"]*(-2)
ll2<-theo.fit["results"]["ll.lin"]*(-2)
1-pchisq(ll1-ll2,1)

ll1<-theo.base["results"]["ll.is"]*(-2)
ll2<-theo.fit["results"]["ll.is"]*(-2)
1-pchisq(ll1-ll2,1)

ll1<-theo.base["results"]["ll.gq"]*(-2)
ll2<-theo.fit["results"]["ll.gq"]*(-2)
1-pchisq(ll1-ll2,1)

theo.fit["results"]["ll.is"]

######################################
# Diagnostic plots

# Plotting individual fits with selected options
par(mfrow=c(2,2))
plot(theo.fit,plot.type="individual.fit",new=FALSE,ilist=1:4,smooth=TRUE,ylog=T, pch=1, col="Blue",xlab="Time in hr",ylab="Theophylline concentrations (mg/L)")
# Plots of the observations versus predictions
plot(theo.fit, plot.type="observations.vs.predictions")

# Scatterplots and distribution of residuals
plot(theo.fit, plot.type="npde")

# VPC
plot(theo.fit, plot.type="vpc")

#####################################################################################
#### Example 3.3 - Simulations from Elodie

setwd(file.path(work.dir,'PaperExamples'))
library(xtable)

######################################################
# Estimation for the simulated dataset
### two scenarios (sparse and rich, with Hill model and Hill factor=3)
### four settings for each scenario (2 sets of options and 2 sets of CI)
scenariosPlan<-c("R3A","S3A")
fitsettings<-c("ptrue_default","pwrong_default","ptrue_tuned","pwrong_tuned")

# Paths
namdir.data<-file.path(work.dir,'PaperExamples',"dataPlan12")
dir.fit<-file.path(work.dir,'PaperExamples',"fitPlan12")
if(save.results & !file_test("-d",dir.fit)) dir.create(dir.fit)

# Hill model
modelhill<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (4 columns, E0, Emax, E50, gamma)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  emax<-psi[id,2]
  e50<-psi[id,3]
  gamma<-psi[id,4]
  f<-e0+emax*dose**gamma/(e50**gamma+dose**gamma)
  return(f)
}

# Initial conditions
ptrue<-c(30,500,2,5,0.490,0.245,0.490,0.090,4)
pwrong<-c(60,1000,1,10,0.100,0.01,0.1,0.1,1)
nsim<-100
prefix.data<-"mc-sim-1-"

# Fitting the datasets
# Note - isolate LL by IS and GQ
#Rprof("profil_fitsJSS.out")
if(estimPlan12) 
  source("paperJSS_fitPlan12.R")

#Rprof(NULL)
#system("R CMD Rprof profil_fitsJSS.out")

######################################################
# Reading the results from the fits and saving to a LaTeX table
dir.create(file.path(dir.fit,"results"))

system("rm latex.res")
xtab2<-xtab.test<-NULL
for(scenar in scenariosPlan) {
  # Simulation values
  ci.param<-ptrue
  if(length(grep("P",scenar))>0) {ci.param[9]<-0.01}
  ci.param[3]<-as.integer(substr(scenar,2,2))	
  rhosim<-ci.param[6]/sqrt(ci.param[5]*ci.param[7])
  true.par<-data.frame(nampar=colnames(par.nocov)[3:12],simval=c(ci.param[c(4,1:2)],ci.param[c(3,9,8,5,7,6)],rhosim))
  true.par[5,2]<-sqrt(true.par[5,2])
  
  # Bias and ST on estimated parameters
  res.bias<-NULL
  for(settings in fitsettings) {
    namres.dir<-paste(scenar,"_",settings,sep="")
    namdir<-file.path(dir.fit,namres.dir)
    par.nocov<-read.table(file.path(namdir,"parnocov.res"),header=T)
    se.nocov<-read.table(file.path(namdir,"senocov.res"),header=T)
    rbias<-par.nocov[,3:12]
    for(icol in 3:12) {
      nampar<-colnames(par.nocov)[icol]
      trueval<-true.par[match(nampar,true.par[,1]),2]
      rbias[,(icol-2)]<-100*(par.nocov[,icol]-trueval)/trueval  
    }
    x1<-summary(rbias,digits=3,nsmall=3,trim=T) # Relative Bias
    x2<-apply(rbias,2,function(x) sqrt(mean(x**2))) # Relative RMSE
    x3<-matrix(unlist(strsplit(x1[4,],"Mean   :")),ncol=2,byrow=T)[,2]
    x4<-format(as.double(x3),digits=2,nsmall=2,drop0trailing =T)
    l2<-paste(x4," (",format(as.double(x2),digits=2,nsmall=2,trim=T,drop0trailing =T),")",sep="")
    x3<-matrix(unlist(strsplit(x1[1,],"Min.   :")),ncol=2,byrow=T)[,2]
    x4<-format(as.double(x3),digits=2,nsmall=2,drop0trailing =T,trim=T)
    x3<-matrix(unlist(strsplit(x1[6,],"Max.   :")),ncol=2,byrow=T)[,2]
    x5<-format(as.double(x3),digits=2,nsmall=2,drop0trailing =T,trim=T)
    xtab<-data.frame(mean=l2,range=paste("[",x4,"--",x5,"]",sep=""))  
    if(is.null(res.bias)) res.bias<-xtab else res.bias<-cbind(res.bias,xtab)
  }
  rownames(res.bias)<-colnames(rbias)
  xtab2<-rbind(xtab2,res.bias[,c(0:3)*2+1])
  nam.file<-paste(scenar,"_resbias.txt",sep="")
  write.table(res.bias,file.path(dir.fit,"results",nam.file),quote=F,row.names=F)
  colnames(xtab2)<-c("True","False","True","False")
  rownames(xtab2)<-c("E$_0$","E$_{\\rm max}$","ED$_{50}$","$\\gamma$","$\\sigma$","$\\omega^2_{E_0}$","$\\omega^2_{E_max}$","$\\omega^2_{ED_{50}}$","$cov(E_{\\rm max},ED_{50})$","$\\rho$")
  tit<-ifelse(scenar=="S3A","sparse sampling","rich sampling")
  print(xtable(xtab2, align="rcccc",caption=paste("Relative bias (RRMSE) in the scenario with",tit),label=paste("tab:RBias.",scenar,sep=""), floating = FALSE),file="latex.res",append=TRUE,include.colnames=T,quote=F,sanitize.text.function=function(x){x})
}

#####################################################################################
#### Example 3.4 - Scalability
#\hl{Eco: TODO} simulations en augmentant le nb de sujets et de prélèvements (N=50, 100, 500). Simulations avec des modèles avec plus d'effets aléatoires (ici p=3-4, augmenter à p=6, p=10).


# simulation example based on simulations performed in: /home/eco/work/monolix/rversion/hotline/junePAGE/rsimul.r

simulate.toggle<-TRUE
estimate.toggle<-TRUE

#########################
### Simulation functions


#########################
### Models
# Emax model
modelemax<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, E50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  emax<-psi[id,2]
  e50<-psi[id,3]
  f<-e0+emax*dose/(e50+dose)
  return(f)
}

# Hill model - already loaded in example 3.3
modelhill

# PK models for 1 cpt (2 par), 2 cpt (4 par) and 3 cpt (6 par) models with bolus absorption 
# parameterised as expontenial models in closed form

expone<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, E50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  time<-xidep[,1]
  dose<-xidep[,2]
  A1<-psi[id,1]
  alp1<-psi[id,2]
  f<-dose*A1*exp(-alp1*time)
  return(f)
}

exptwo<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, E50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  time<-xidep[,1]
  dose<-xidep[,2]
  A1<-psi[id,1]
  alp1<-psi[id,2]
  A2<-psi[id,3]
  alp2<-psi[id,4]
  f<-dose*(A1*exp(-alp1*time)+A2*exp(-alp2*time))
  return(f)
}

expthree<-function(psi,id,xidep) {
  # input:
  #   psi : matrix of parameters (3 columns, E0, Emax, E50)
  #   id : vector of indices 
  #   xidep : dependent variables (same nb of rows as length of id)
  # returns:
  #   a vector of predictions of length equal to length of id
  time<-xidep[,1]
  dose<-xidep[,2]
  A1<-psi[id,1]
  alp1<-psi[id,2]
  A2<-psi[id,3]
  alp2<-psi[id,4]
  A3<-psi[id,5]
  alp3<-psi[id,6]
  f<-dose*(A1*exp(-alp1*time)+A2*exp(-alp2*time)+A3*exp(-alp3*time))
  return(f)
}

###########################
# Error model
error.model<-"combined"

# Parameter values - Hill and Emax models
ptrue<-c(5,30,500,2,0.09,0.25,0.25,0.09,4,0.01)
pwrong<-c(10,60,1000,1,0.1,0.1,0.1,0.1,1,0.0625)
pres<-sqrt(ptrue[9:10])

# Parameter values - exponential models
A1<-100
A2<-75
A3<-50
alp1<-1
alp2<-0.4
alp3<-0.2
# alp1<-2
# alp2<-0.7
# alp3<-0.2

# Times for the model with 3 exponentials
times.exp<-c(0,1,3,6,12,24) # 6 times, optimised with PFIM, best one-group protocol
#times.exp<-c(0,1,2,3,4,6,12,24) # 8 times

# x<-seq(0,12,0.5)  
# plot(x,(A1*exp(-alp1*x)+A2*exp(-alp2*x)+A3*exp(-alp3*x)),type="l",lwd=3)
# lines(x,(A1*exp(-alp1*x)+A2*exp(-alp2*x)),col='DarkRed',lwd=2)
# lines(x,(A1*exp(-alp1*x)),col='DarkBlue',lwd=2)

# General settings
nsim<-100

# Directories
base.data.sim<-file.path(work.dir,"PaperExamples","simulations","data")
dir.create(file.path(work.dir,"PaperExamples","simulations"))
dir.create(base.data.sim)
base.data.prof<-file.path(work.dir,"PaperExamples","simulations","profiling")
dir.create(base.data.prof)

#########################
### Scalability - varying the number of subjects and number of samples

name.sim<-"samplesize"

# Settings
Nsuj<-c(30,100,200,500)

fulldoses<-c(0,10,50,75,100,200,300,500,600,800,900,1000)
dose.set1<-which(fulldoses %in% c(0,300,1000))
dose.set2<-which(fulldoses %in% c(0,50,100,300,500,1000))

### Simulate data
nisuj<-1
niobs<-1
for(nisuj in 1:4) {
  for(niobs in 1:3) {
    nsuj<-Nsuj[nisuj]
    xvar<-switch(niobs,fulldoses[dose.set1],fulldoses[dose.set2],fulldoses)

    # create directory to save data
    dir.data.sim<-file.path(base.data.sim,paste(name.sim,"_N",nsuj,"n",length(xvar),sep=""))
    dir.create(dir.data.sim)
    
    if(simulate.toggle) {
      for(irep in 1:nsim) {
        # simulate and save psi
        psi<-NULL
        for(i in 1:4) psi<-cbind(psi,ptrue[i]*exp(rnorm(nsuj,sd=sqrt(ptrue[(i+4)]))))
        psi<-round(psi,2)
        # if sparse design (nobs=3), assign 100/300 to half of the subjects each, all subjects have 0 and 1000
        colnames(psi)<-c("E0","Emax","E50","gamma")
        if(niobs==1) {
          xvar2<-xvar; xvar2[2]<-100
          xvar1<-c(rep(xvar,nsuj/2),rep(xvar2,nsuj/2))
        } else xvar1<-rep(xvar,nsuj)
        xidep<-matrix(data=xvar1,ncol=1)
        id<-rep(c(1:nsuj),each=length(xvar))
        ypred<-modelhill(psi,id,xidep)
        # check
        # plot(xidep[,1],ypred)
        # lines(xvar,ptrue[1]+ptrue[2]*xvar**ptrue[4]/(ptrue[3]**ptrue[4]+xvar**ptrue[4]),col="DarkRed",lwd=2)
        
        # simulate and save datasets
        eps<-rnorm(length(ypred),sd=pres[1])+rnorm(length(ypred),sd=pres[2])*abs(ypred)
        ysim<-ypred+eps
        simdata<-data.frame(id=id,dose=xidep[,1],response=round(ysim,2))
        
        nam.file.dat<-paste("simdata",irep,".dat",sep="")
        nam.file.psi<-paste("simpsi",irep,".dat",sep="")
        write.table(simdata,file.path(dir.data.sim,nam.file.dat),quote=F,row.names=F,na=".")
        write.table(psi,file.path(dir.data.sim,nam.file.psi),quote=F,row.names=F,na=".")
      }
    }
    
    nam.profile<-file.path(base.data.prof,paste("memory_",name.sim,"_N",nsuj,"n",length(xvar),".out",sep=""))
    Rprof(nam.profile)
    if(estimate.toggle) {
      xtab<-NULL
      base.data.fit<-file.path(work.dir,"PaperExamples","simulations","fit")
      hillmodel<-saemixModel(model=modelhill,description="Emax model with sigmoidicity (Hill model)",
              psi0=matrix(c(ptrue[1:4],0,0,0,0), ncol=4,byrow=T, dimnames=list(NULL,c("E0","Emax","EC50","gamma"))), 
              transform.par=c(1,1,1,1), omega.init=diag(ptrue[5:8]), error.init=pres, error.model=error.model)
      saemix.options<-list(displayProgress=FALSE, save.graphs=FALSE, directory=base.data.fit)
      for(irep in 1:nsim) {
        nam.file.dat<-paste("simdata",irep,".dat",sep="")
        saemix.data<-saemixData(name.data=file.path(dir.data.sim,nam.file.dat),header=T,name.group="id",name.predictors="dose",name.response="response",units=list(x="mg",y="-"))
        yfit<-saemix(hillmodel,saemix.data,saemix.options)
        nam.file.fit<-paste("parameters",irep,".res",sep="")
        cmd<-paste("cp ",file.path(base.data.fit,"pop_parameters.txt"),file.path(dir.data.sim,nam.file.fit))
        system(cmd)
        nam.file.fit<-paste("indiv",irep,".res",sep="")
        cmd<-paste("cp ",file.path(base.data.fit,"indiv_parameters.txt"),file.path(dir.data.sim,nam.file.fit))
        system(cmd)
        l1<-c(irep,coef(yfit)$fixed,diag(yfit@results@omega),yfit@results@respar,yfit@results@se.fixed, yfit@results@se.omega,yfit@results@se.respar)
        xtab<-rbind(xtab,l1)
      }
      colnames(xtab)[1]<-"Replication"
      colnames(xtab)[6:9]<-paste("omega",colnames(xtab)[2:5],sep=".")
      colnames(xtab)[10:11]<-c("a","b")
      colnames(xtab)[12:21]<-paste("SE",colnames(xtab)[2:11],sep=".")
      write.table(xtab,file.path(dir.data.sim,"table_parestim.res"),quote=F,row.names=F,na=".")
    }
    Rprof()

    # Memory profiling
    nam.profile.summary<-file.path(base.data.prof,paste("memory_",name.sim,"_N",nsuj,"n",length(xvar),".summary",sep=""))
    system(paste("R CMD Rprof",nam.profile,">",nam.profile.summary))
  }
}

# Bias/RMSE
restab<-NULL
for(nisuj in 1:4) {
  for(niobs in 1:3) {
    nsuj<-Nsuj[nisuj]
    xvar<-switch(niobs,fulldoses[dose.set1],fulldoses[dose.set2],fulldoses)
    dir.data.sim<-file.path(base.data.sim,paste(name.sim,"_N",nsuj,"n",length(xvar),sep=""))
    xtab<-read.table(file.path(dir.data.sim,"table_parestim.res"),header=T,na=".")
    l1<-c(nsuj,length(xvar))
    for(i in 1:8) l1<-c(l1,mean((xtab[,(i+1)]-ptrue[i]))/ptrue[i])
    for(i in 9:10) l1<-c(l1,mean((xtab[,(i+1)]-sqrt(ptrue[i])))/sqrt(ptrue[i]))
    restab<-rbind(restab,l1)
  }
}
colnames(restab)<-c("Nsuj","nsamp",colnames(xtab[2:11]))

# Simulation N=200, n=3 ???
# check distribution of times ?

########################################
# Scalability - Model complexity (1, 2, 3 exponents)

name.sim<-"complexity"

# Settings
nsuj<-100
sd.iiv<-.3 # 40% IIV on all parameters
exp.par<-c(A1,alp1,A2,alp2,A3,alp3)
pres<-c(5,0)
# Same datasets for the 3 models, predictions using the model with 3 exponents

# create directory to save data
#for(i in 1:3) dir.data.sim<-file.path(base.data.sim,paste(name.sim,"_exp",i,sep=""))
dir.data.sim<-file.path(base.data.sim,name.sim)
dir.create(dir.data.sim)

if(simulate.toggle) {
  for(irep in 1:nsim) {
    # simulate and save psi
    psi<-NULL
    for(i in 1:6) psi<-cbind(psi,exp.par[i]*exp(rnorm(nsuj,sd=sd.iiv)))
    psi<-round(psi,2)
    colnames(psi)<-c("A1","alp1","A2","alp2","A3","alp3")
    xidep<-matrix(data=rep(times.exp,nsuj),ncol=1)
    xidep<-cbind(xidep,rep(1,dim(xidep)[1]))
    id<-rep(c(1:nsuj),each=length(times.exp))
    ypred<-expthree(psi,id,xidep)
    # check
    # plot(xidep[,1],ypred)
    # xvar<-xidep[,1]
    # lines(xvar,exp.par[1]*exp(-exp.par[2]*xvar),col="DarkRed",lwd=2)
    
    # simulate and save datasets
    eps<-rnorm(length(ypred),sd=pres[1])
    ysim<-ypred+eps
    simdata<-data.frame(id=id,time=xidep[,1],dose=xidep[,2],response=round(ysim,2))
    
    nam.file.dat<-paste("simdata",irep,".dat",sep="")
    nam.file.psi<-paste("simpsi",irep,".dat",sep="")
    write.table(simdata,file.path(dir.data.sim,nam.file.dat),quote=F,row.names=F,na=".")
    write.table(psi,file.path(dir.data.sim,nam.file.psi),quote=F,row.names=F,na=".")
  }
}

modsmx.exp1<-saemixModel(model=expone,description="Exponential model - 1 exponential",
        psi0=matrix(c(exp.par[1:2],0,0), ncol=2,byrow=T, dimnames=list(NULL,c("A1","alp1"))),
        transform.par=c(1,1)) #,error.model="combined")
modsmx.exp2<-saemixModel(model=exptwo,description="Exponential model - 1 exponential",
        psi0=matrix(c(exp.par[1:4],0,0,0,0), ncol=4,byrow=T, dimnames=list(NULL,c("A1","alp1","A2","alp2"))),
        transform.par=c(1,1,1,1)) #,omega.init=diag(rep(0.09,4)),error.model="combined")
modsmx.exp3<-saemixModel(model=expthree,description="Exponential model - 1 exponential",
        psi0=matrix(c(exp.par[1:6],rep(0,6)), ncol=6,byrow=T, dimnames=list(NULL,c("A1","alp1","A2","alp2","A3","alp3"))),
        transform.par=c(1,1,1,1,1,1)) # ,error.model="combined")    # omega.init=diag(rep(0.09,6)),
base.data.fit<-file.path(work.dir,"PaperExamples","simulations","fit")
# saemix.options<-list(displayProgress=TRUE, save.graphs=TRUE, directory=base.data.fit)
saemix.options<-list(displayProgress=FALSE, save.graphs=FALSE, directory=base.data.fit)

for(imod in 1:3) {
  saemix.model<-switch(imod,modsmx.exp1,modsmx.exp2,modsmx.exp3)
  nam.profile<-file.path(base.data.prof,paste("memory_",name.sim,"_exp",imod,".out",sep=""))
  Rprof(nam.profile)
  if(estimate.toggle) {
    dir.data.fit<-file.path(base.data.fit,paste(name.sim,"_exp",imod,sep=""))
    dir.create(dir.data.fit)
    xtab<-NULL
    for(irep in 1:nsim) {
      nam.file.dat<-paste("simdata",irep,".dat",sep="")
      saemix.data<-saemixData(name.data=file.path(dir.data.sim,nam.file.dat),header=T,name.group="id",name.predictors=c("time","dose"),name.response="response",units=list(x="hr",y="mg/L"))
      yfit<-saemix(saemix.model,saemix.data,saemix.options)
      nam.file.fit<-paste("parameters",irep,".res",sep="")
      cmd<-paste("cp ",file.path(base.data.fit,"pop_parameters.txt"),file.path(dir.data.fit,nam.file.fit))
      system(cmd)
      nam.file.fit<-paste("indiv",irep,".res",sep="")
      cmd<-paste("cp ",file.path(base.data.fit,"indiv_parameters.txt"),file.path(dir.data.fit,nam.file.fit))
      system(cmd)
      l1<-c(irep,coef(yfit)$fixed,diag(yfit@results@omega),yfit@results@respar,yfit@results@se.fixed, yfit@results@se.omega,yfit@results@se.respar)
      xtab<-rbind(xtab,l1)
    }
    colnames(xtab)[1]<-"Replication"
    npar<-2*imod
    colnames(xtab)[(2+npar):(1+2*npar)]<-paste("omega",colnames(xtab)[2:(npar+1)],sep=".")
    colnames(xtab)[(2*npar+2):(2*npar+3)]<-c("a","b")
    colnames(xtab)[(2*npar+4):(4*npar+5)]<-paste("SE",colnames(xtab)[2:(2*npar+3)],sep=".")
    write.table(xtab,file.path(dir.data.fit,"table_parestim.res"),quote=F,row.names=F,na=".")
  }
  Rprof()
}

# Memory profiling
for(imod in 1:3) {
  nam.profile<-file.path(base.data.prof,paste("memory_",name.sim,"_exp",imod,".out",sep=""))
  nam.profile.summary<-file.path(base.data.prof,paste("memory_",name.sim,"_exp",imod,".summary",sep=""))
  system(paste("R CMD Rprof",nam.profile,">",nam.profile.summary))
}
########################################
# Scalability - Number of fixed parameters
# Adding covariate effects

########################################
# Graphical results

name.sim<-"complexity"
tprof<-c()
for(imod in 1:3) {
  nam.profile<-file.path(base.data.prof,paste("memory_",name.sim,"_exp",imod,".out",sep=""))
  nam.profile.summary<-file.path(base.data.prof,paste("memory_",name.sim,"_exp",imod,".summary",sep=""))
  cmd<-paste("grep 'Total run time'",nam.profile.summary)
  x<-system(cmd,intern=T)
  x<-unlist(strsplit(unlist(strsplit(x,"Total run time: "))," seconds."))
  tprof<-c(tprof,as.double(x))
}
name.sim<-"samplesize"
for(nisuj in 1:4) {
  for(niobs in 1:3) {
    nsuj<-Nsuj[nisuj]
    xvar<-switch(niobs,fulldoses[dose.set1],fulldoses[dose.set2],fulldoses)
    nam.profile.summary<-file.path(base.data.prof,paste("memory_",name.sim,"_N",nsuj,"n",length(xvar),".summary",sep=""))
    cmd<-paste("grep 'Total run time'",nam.profile.summary)
    x<-system(cmd,intern=T)
    x<-unlist(strsplit(unlist(strsplit(x,"Total run time: "))," seconds."))
    tprof<-c(tprof,as.double(x))
  }
}
tprof.ss<-matrix(tprof[-c(1:3)],byrow=T,ncol=3)
rownames(tprof.ss)<-Nsuj
colnames(tprof.ss)<-c(3,6,12)

if(save.results) {
  paper.dir<-"/home/eco/xtex/saemix/saemixJSS/reviewJSS/figures"
#  postscript(file.path(paper.dir,"runtimes.eps"),horizontal=T)
  pdf(file.path(paper.dir,"runtimes.pdf"),paper="a4r")
  par(mfrow=c(1,2))
  zecol<-c("DarkRed","DarkBlue","DarkGreen")
  plot(Nsuj,tprof.ss[,1],pch=20,type="n",ylim=c(0,max(tprof.ss)),xlab="Number of subjects", ylab="Runtimes for 100 simulations (s)",lwd=2,col=zecol[1])
  for(i in 1:3) {
    points(Nsuj,tprof.ss[,i],pch=20,col=zecol[i])
    lines(Nsuj,tprof.ss[,i],lwd=2,col=zecol[i])
  }
  legend(30,max(tprof.ss),paste(c(3,6,12),"samples"),col=zecol,pch=20,lwd=2)
  plot(c(2,4,6),tprof[1:3],type="b",xlab="Number of random effectds", ylab="Runtimes for 100 simulations (s)",lwd=2,pch=20)
  lines(c(2,4,6),tprof[1:3])
  dev.off()
}

#####################################################################################

