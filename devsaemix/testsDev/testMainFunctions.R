####################################################################################
####				Tests of the main function			####
####################################################################################
setwd("/home/eco/work/monolix/rversion/current/libcov/covsaemix/testsDev/")

############################
# Defining generic for new methods
source("../R/global.R")

############################
# SaemixData
source("../R/SaemixData.R")

# SaemixRes
source("../R/SaemixRes.R")

# SaemixModel
source("../R/SaemixModel.R")

# SaemixObject
source("../R/SaemixObject.R")

# Functions
source("../R/func_main.R")
source("../R/func_aux.R")
source("../R/func_plots.R")
source("../R/compute_LL.R")
source("../R/main_initialiseMainAlgo.R")
source("../R/main_estep.R")
source("../R/main_mstep.R")

####################################################################################
####					Theophylline				####
####################################################################################

##### Model function

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

dat<-read.table("../data/theo.saemix.tab",header=T)
dat$Weight[dat$Id==3]<-NA

##### Data
xdat<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))

##### Model with combined residual error
ymod<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

sopt<-list(map=FALSE,fim=F,ll.is=FALSE,save.graphs=FALSE,save=FALSE,seed=12345)
myfit<-saemix(ymod,xdat,sopt)

if(FALSE) {
	sopt<-list(save.graphs=TRUE,save=TRUE,directory="/home/eco/work/monolix/rversion/current/libcov/betatest/theoCovDebug")
	myfit<-saemix(ymod,xdat,sopt)
}
sopt<-list(save.graphs=FALSE,save=FALSE,seed=12345)
myfit<-saemix(ymod,xdat,sopt)


combfit<-saemix(ymod,xdat)

#model<-ymod;data<-xdat;control<-sopt

##### Model with exponential error
ymod.exp<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,1.8,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(0,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="exponential")

yexp<-saemix(ymod.exp,xdat)
# yexp<-saemix(ymod.exp,xdat,list(save.graphs=FALSE,save=FALSE))

plot(yexp,plot.type="data")
plot(yexp,plot.type="convergence")
plot(yexp,plot.type="likelihood")
plot(yexp,plot.type="population.fit")
plot(yexp,plot.type="individual.fit")
plot(yexp,plot.type="observations.vs.predictions")
plot(yexp,plot.type="residuals.scatter")
plot(yexp,plot.type="residuals.distribution")
plot(yexp,plot.type="random.effects")
plot(yexp,plot.type="correlations")
plot(yexp,plot.type="parameters.vs.covariates")
plot(yexp,plot.type="randeff.vs.covariates")
plot(yexp,plot.type="marginal.distribution")
plot(yexp,plot.type="npde")
plot(yexp,plot.type="vpc")

dat<-read.table("../data/PD1.saemix.tab",header=T)
vec<-ifelse(dat$gender==1,"W","M")
dat$gender<-as.factor(vec)
xchar<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))

####################################################################################
####					Emax simulations			####
####################################################################################
dat<-read.table("../data/PD1.saemix.tab",header=T)
vec<-ifelse(dat$gender==1,"W","M")
dat$gender<-as.factor(vec)
saemix.data1<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"))

modelemax<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, E0, Emax, EC50)
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

saemix.model<-saemixModel(model=modelemax,description="Emax model", psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined")

saemix.options<-list(nb.chains=1, save=FALSE,save.graphs=FALSE)

saemix.fit1<-saemix(saemix.model,saemix.data1,saemix.options)

plot(saemix.fit1,plot.type="data")
plot(saemix.fit1,plot.type="convergence")
plot(saemix.fit1,plot.type="likelihood")
plot(saemix.fit1,plot.type="population.fit",ilist=1:12)
plot(saemix.fit1,plot.type="individual.fit",ilist=1:12,smooth=T)
plot(saemix.fit1,plot.type="observations.vs.predictions")
plot(saemix.fit1,plot.type="residuals.scatter")
plot(saemix.fit1,plot.type="residuals.distribution")
plot(saemix.fit1,plot.type="random.effects")
plot(saemix.fit1,plot.type="correlations")
plot(saemix.fit1,plot.type="parameters.vs.covariates")
plot(saemix.fit1,plot.type="randeff.vs.covariates")
plot(saemix.fit1,plot.type="marginal.distribution")
plot(saemix.fit1,plot.type="npde")
plot(saemix.fit1,plot.type="vpc")

####################################################################################
