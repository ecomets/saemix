#####################################################################################
# Library
setwd("/home/eco/work/monolix/rversion/current/saemix")

############################
# Defining generic for new methods
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

# Toggle save feature
save.results<-TRUE
save.results<-FALSE
save.dir<-"Demo"
if(save.results & !file_test("-d",save.dir)) dir.create(save.dir)

#####################################################################################
# Theophylline

# Data - changing gender to M/F
theo.saemix<-read.table("devsaemix/data/theo.saemix.tab",header=T,na=".")
theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")


# Doc
# theo.saemix<-read.table("devsaemix/data/theo.saemix.tab",header=T,na=".")
# saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

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
# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))

estonly<-list(seed=39546,save=save.results,save.graphs=save.results,map=F,fim=F,ll.is=F,directory=file.path(save.dir,"theoCov2"))
theo.onlypop<-saemix(saemix.model,saemix.data,estonly)

saemix.options<-list(seed=632545,save=save.results,save.graphs=save.results,directory=file.path(save.dir,"theoNoCov"))

Rprof("Demo/theoCF1cpt.out")
theo.nocov<-saemix(saemix.model,saemix.data,saemix.options)
Rprof(NULL)

plot(saemix.fit,plot.type="individual")

# Model with covariates
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

saemix.options<-list(seed=39546,save=save.results,save.graphs=save.results,directory=file.path(save.dir,"theoCov1"))
cov.fit<-saemix(saemix.model,saemix.data,saemix.options)

estonly<-list(seed=39546,save=save.results,save.graphs=save.results,map=F,fim=F,lls.is=F,directory=file.path(save.dir,"theoCov2"))
cov.fit<-saemix(saemix.model,saemix.data,estonly)
fim.saemix(cov.fit)

# Changing gender to M/F
theo.saemix<-read.table("devsaemix/data/theo.saemix.tab",header=T,na=".")
theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

saemix.options<-list(seed=39546,save=save.results,save.graphs=save.results,directory=file.path(save.dir,"theoCov3"))
cov.fit2<-saemix(saemix.model,saemix.data,saemix.options)

# One missing covariate
theo.saemix<-read.table("devsaemix/data/theo.saemix.tab",header=T,na=".")
theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
theo.saemix$Weight[theo.saemix$Id==3]<-NA
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,1,1),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

saemix.options<-list(seed=39546,save=save.results,save.graphs=save.results,directory=file.path(save.dir,"theoCov4"))
cov.fit3<-saemix(saemix.model,saemix.data,saemix.options)
summary(cov.fit3)
summary(cov.fit3@data) # Check 11 subjects

#####################################################################################
# Theophylline with ODE model
library(deSolve)

# Data - changing gender to M/F
theo.saemix<-read.table("devsaemix/data/theo.saemix.tab",header=T,na=".")
theo.saemix$Sex<-ifelse(theo.saemix$Sex==1,"M","F")
data.theo<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

########################
# ad-hoc implementation
# ODE system for 1st order elimination - setting up deSolve
ode1cpt<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dGut <- (-ka)*Gut
    dCentral <- ka*Gut - CL/V*Central
    return(list(c(dGut, dCentral)))
  })
}

# testing ode1cpt
par1cpt<-c(ka=1.53,   # first-order absorption
         V=31.5,   # volume of central compartment
         CL=4.1)   # clearance

modelODE1cpt<-function(psi,id,xidep) { 
  ypred<-rep(0,length(id))
  yini<-c(Gut=0,Central=0)
  for(isuj in unique(id)) {
    param<-c(ka=psi[isuj,1],V=psi[isuj,2],CL=psi[isuj,3])
    yini[1]<-xidep[match(isuj,id),1] # Dose
    tim<-c(0,xidep[id==isuj,2]) # Need to integrate from time=0
    out<-ode(yini,tim,ode1cpt,param)
    ypred[id==isuj]<-out[-c(1),3]/param[2]
  }
  return(ypred)
}

yini<-c(Gut=320,Central=0)
xtim<-seq(0,25,by=0.2)
out1<-ode(yini,xtim,ode1cpt,par1cpt)

plot(data.theo)
lines(xtim,out1[,3]/par1cpt[2],lwd=3,col="DarkRed")

yini[1]*param[1]/(param[2]*(param[1]-param[3]/param[2]))*(exp(-param[3]/param[2]*tim)-exp(-param[1]*tim))

# Compare saemix model functions for 1cpt and ODE- should get same results
xidep<-saemix.data@data[,c("Dose","Time")]
id<-saemix.data@data[,c("Id")]
psi<-exp(theo.nocov@results@phi)

ypred.ode<-modelODE1cpt(psi,id,xidep)
ypred.cf<-model1cpt(psi,id,xidep)

plot(ypred.ode,ypred.cf)
abline(0,1)
summary(ypred.ode-ypred.cf)
summary((ypred.ode-ypred.cf)/(ypred.ode+ypred.cf)*2)

########################
# saemix run
model.ode1cpt<-saemixModel(model=modelODE1cpt,description="One-compartment model with first-order absorption, ODE",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))
saemix.options<-list(seed=632545,save=save.results,save.graphs=save.results,directory=file.path(save.dir,"theoODE1cpt"))

Rprof("Demo/theoODE1cpt.out")
theoODE.fit<-saemix(model.ode1cpt,data.theo,saemix.options)
Rprof(NULL)

# Comparing run times w/r closed-form solution
system("R CMD Rprof Demo/theoCF1cpt.out")
system("R CMD Rprof Demo/theoODE1cpt.out")

########################
# MM model
# ODE system for MM
odeMM<-function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dGut <- (-ka)*Gut
    dCentral <- ka*Gut - Vm*Central/(Km/Vc+Central)
    return(list(c(dGut, dCentral)))
  })
}

parMM<-c(ka=1.5,   # first-order absorption
         Vc=31.5,   # volume of central compartment
         Vm=15,   # 
         Km=5)   # 

# testing odeMM
yini<-c(Gut=320,Central=0)
outMM<-ode(yini,xtim,odeMM,parMM)
plot(saemix.data)
lines(xtim,out1[,3]/par1cpt[2],lwd=3,col="DarkRed")
lines(xtim,outMM[,3]/parMM[2],lwd=3,col="DarkBlue")

# saemix function
modelMM<-function(psi,id,xidep) { 
  ypred<-rep(0,length(id))
  yini<-c(Gut=0,Central=0)
  for(isuj in unique(id)) {
    parMM<-c(ka=psi[isuj,1],Vc=psi[isuj,2],Vm=psi[isuj,3],Km=psi[isuj,4])
    yini[1]<-xidep[match(isuj,id),1] # Dose
    tim<-c(0,xidep[id==isuj,2]) # Need to integrate from time=0
    out<-ode(yini,tim,odeMM,parMM)
    ypred[id==isuj]<-out[-c(1),3]/parMM[2]
  }
  return(ypred)
}

# Testing modelMM
xidep<-saemix.data@data[,c("Dose","Time")]
id<-saemix.data@data[,c("Id")]
psi<-exp(theo.nocov@results@phi)
psi<-cbind(psi,psi[,3]*4)
psi[,3]<-psi[,2]/1.5

ypred<-modelMM(psi,id,xidep)
plot(saemix.data)
lines(xidep[,2],ypred,col='DarkRed',lwd=2)

# MM model with saemix
saemix.model<-saemixModel(model=modelMM,description="Michaelis-Menten elimination with deSolve",psi0=matrix(c(1.,20,15,4,0,0,0,-0.01),ncol=4,byrow=TRUE, dimnames=list(NULL, c("ka","V","Vm","Km"))),transform.par=c(1,1,1,1),error.model="combined")
theoMM.fit<-saemix(saemix.model,saemix.data,saemix.options)

#####################################################################################
# Simulated PD data

PD1.saemix<-read.table("devsaemix/data/PD1.saemix.tab",header=T)
PD2.saemix<-read.table("devsaemix/data/PD2.saemix.tab",header=T)
saemix.data<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"), units=list(x="mg",y="-",covariates=c("-")))

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
saemix.model<-saemixModel(model=modelemax,description="Emax growth model",psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, c("E0","Emax","EC50"))),transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

saemix.options<-list(nb.chains=1,seed=42,save=save.results,save.graphs=save.results,directory=file.path(save.dir,"PD1cov"))

# Plotting the data
plot(saemix.data,main="Simulated data PD1")

saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
plot(saemix.fit,plot.type="individual",ilist=1:4,smooth=T,cex.lab=0.8,cex.axis=0.8,cex.main=0.8)

# Compare models with and without covariates with LL by Importance Sampling
# SE not computed
model1<-saemixModel(model=modelemax,description="Emax growth model", psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1),covariate.model=matrix(c(0,0,0), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

model2<-saemixModel(model=modelemax,description="Emax growth model",psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

saemix.options<-list(map=F,nb.chains=3,seed=765754, nbiter.saemix=c(500,300),save=save.results,save.graphs=save.results, directory=file.path(save.dir,"PD1nocov2"))
fit1<-saemix(model1,saemix.data,saemix.options)

saemix.options<-list(map=F,nb.chains=3,seed=765754, nbiter.saemix=c(500,300),save=save.results,save.graphs=save.results, directory=file.path(save.dir,"PD1cov2"))
fit2<-saemix(model2,saemix.data,saemix.options)
wstat<-(-2)*(fit1["results"]["ll.is"]-fit2["results"]["ll.is"])

cat("LRT test for covariate effect on EC50: p-value=",1-pchisq(wstat,1),"\n")

# Dataset with only one point for some subjects and one subject with only NA
PD1.saemix<-read.table("devsaemix/data/PD1.saemix.tab",header=T)
ione<-c(75,83,95)
isamp<-c()
for(i in ione) {
	idx<-which(PD1.saemix$subject==i)
	isamp<-c(isamp,idx[!idx%in%sample(idx,1)])
}
PD1.saemix<-PD1.saemix[-isamp,]
PD1.saemix[PD1.saemix$subject==33,"response"]<-NA
saemix.data<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),name.predictors=c("dose"),name.response=c("response"),name.covariates=c("gender"), units=list(x="mg",y="-",covariates=c("-")))
summary(saemix.data)
saemix.options<-list(map=F,nb.chains=3,seed=765754, nbiter.saemix=c(500,300),save=save.results,save.graphs=save.results, directory=file.path(save.dir,"PD1nocov3"))
fit1.bis<-saemix(model1,saemix.data,saemix.options)
summary(fit1.bis@data)

fit1@results@fixed.effects
fit1.bis@results@fixed.effects

#####################################################################################
#### Exponential error model
saemix.data1<-saemixData(name.data="devsaemix/data/PD1.saemix.tab",header=TRUE, name.group=c("subject"),name.predictors=c("dose"),name.response=c("response"), name.covariates=c("gender"), units=list(x="mg",y="-",covariates="-"))

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

saemix.model<-saemixModel(model=modelemax,description="Emax model", psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,1),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="exponential")

saemix.options<-list(nb.chains=1, save=save.results,save.graphs=save.results, directory=file.path(save.dir,"PD1cov3"))

saemix.fit1<-saemix(saemix.model,saemix.data1,saemix.options)

plot(saemix.fit1,plot.type="individual.fit",ilist=1:12)
plot(saemix.fit1,plot.type="population.fit",ilist=1:12)
plot(saemix.fit1,plot.type="both",ilist=1:12)

plot(saemix.fit1,plot.type="vpc")

#####################################################################################
# Weight gain of cows

cow.saemix<-read.table("devsaemix/data/cow.saemix.tab",header=T)
saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("cow"), name.predictors=c("time"),name.response=c("weight"),name.covariates=c("birthyear","twin","birthrank"),units=list(x="days",y="kg",covariates=c("yr","-","-")))

growthcow<-function(psi,id,xidep) {
	# input:
	#   psi : matrix of parameters (3 columns, a, b, k)
	#   id : vector of indices 
	#   xidep : dependent variables (same nb of rows as length of id)
	# returns:
	#   a vector of predictions of length equal to length of id
	x<-xidep[,1]
	a<-psi[id,1]
	b<-psi[id,2]
	k<-psi[id,3]
	f<-a*(1-b*exp(-k*x))
	return(f)
}
saemix.model<-saemixModel(model=growthcow,description="Exponential growth model", psi0=matrix(c(700,0.9,0.02,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

saemix.options<-list(nb.chains=1,nbiter.saemix=c(200,100),  seed=4526,save=save.results,save.graphs=save.results, directory=file.path(save.dir,"cow"))

# Plotting the data
plot(saemix.data,xlab="Time (day)",ylab="Weight of the cow (kg)")

saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

#####################################################################################
# Yield

yield.saemix<-read.table("devsaemix/data/yield.saemix.tab",header=T)
saemix.data<-saemixData(name.data=yield.saemix,header=TRUE,name.group=c("site"),
  name.predictors=c("dose"),name.response=c("yield"),
  name.covariates=c("soil.nitrogen"),units=list(x="kg/ha",y="t/ha",
  covariates=c("kg/ha")))

# Model: linear + plateau
yield.LP<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, ymax, xmax, slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  ymax<-psi[id,1]
  xmax<-psi[id,2]
  slope<-psi[id,3]
  f<-ymax+slope*(x-xmax)
#  cat(length(f),"  ",length(ymax),"\n")
  f[x>xmax]<-ymax[x>xmax]
  return(f)
}
saemix.model<-saemixModel(model=yield.LP,description="Linear plus plateau model", 
  psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,
  c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")

saemix.options<-list(nb.chains=1,seed=666,save=save.results,save.graphs=save.results, directory=file.path(save.dir,"yield1"))

# Plotting the data
plot(saemix.data,xlab="Fertiliser dose (kg/ha)", ylab="Wheat yield (t/ha)")

saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Plotting the data
plot(saemix.data,xlab="Fertiliser dose (kg/ha)", ylab="Wheat yield (t/ha)")

# Comparing the likelihoods obtained by linearisation and importance sampling 
# to the likelihood obtained by Gaussian Quadrature
 saemix.fit<-llgq.saemix(saemix.fit)
{
 cat("LL by Importance sampling, LL_IS=",saemix.fit["results"]["ll.is"],"\n")
 cat("LL by linearisation, LL_lin=",saemix.fit["results"]["ll.lin"],"\n")
 cat("LL by Gaussian Quadrature, LL_GQ=",saemix.fit["results"]["ll.gq"],"\n")
}

# Testing for an effect of covariate soil.nitrogen on Xmax
saemix.model2<-saemixModel(model=yield.LP,description="Linear plus plateau model", 
  psi0=matrix(c(8,100,0.2,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, 
  c("Ymax","Xmax","slope"))),covariate.model=matrix(c(0,1,0),ncol=3,byrow=TRUE), 
  transform.par=c(0,0,0),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, 
  byrow=TRUE),error.model="constant")
saemix.options<-list(nb.chains=1,seed=666,save=save.results,save.graphs=save.results, directory=file.path(save.dir,"yield2"))

 saemix.fit2<-saemix(saemix.model2,saemix.data,saemix.options)
# BIC for the two models
{
  cat("Model without covariate, BIC=",saemix.fit["results"]["bic.is"],"\n")
  cat("Model with covariate, BIC=",saemix.fit2["results"]["bic.is"],"\n")
  pval<-1-pchisq(-2*saemix.fit["results"]["ll.is"]+2*saemix.fit2["results"]["ll.is"],1)
  cat("        LRT: p=",pval,"\n")
}

#####################################################################################
# Logit transformation

# Simulating an Imax model

modelimax<-function(psi,id,xidep) {
	# input:
	#   psi : matrix of parameters (3 columns, E0, Imax, EC50)
	#   id : vector of indices 
	#   xidep : dependent variables (same nb of rows as length of id)
	# returns:
	#   a vector of predictions of length equal to length of id
	dose<-xidep[,1]
	e0<-psi[id,1]
	imax<-psi[id,2]
	e50<-psi[id,3]
	f<-e0*(1-imax*dose/(e50+dose))
	return(f)
}

xdos<-c(0,10,60,90)
nsuj<-50
param<-c(100,0.8,15)
omega<-rep(0.3,3)
sig<-c(4,0.2)

tab<-data.frame(e0=param[1]*exp(rnorm(nsuj,sd=omega[1])),imax=1/(1+exp(-rnorm(nsuj,mean=param[2],sd=omega[1]))),e50=param[3]*exp(rnorm(nsuj,sd=omega[3])))
id<-rep(1:nsuj,each=length(xdos))

x1<-param[2]+rnorm(nsuj,sd=omega[1])
summary(exp(x1)/(1+exp(x1)))

##################################################################################### 
##### Orange tree data - one issue

orange<-saemixData(name.data=Orange,name.group="Tree",name.predictors="age", name.response="circumference",units=list(x="Days",y="mm"))

logistic.model<-function(psi,id,xidep) { 
	age<-xidep[,1]  
	phi1<-psi[id,1]
	phi2<-psi[id,2]
	phi3<-psi[id,3]
	resp<-phi1/(1+exp(-(age-phi2)/phi3))
	return(resp) 
}

# OK: IIV on phi1 only
init.psi<-c(200,800,400)
init.omega<-diag(c(rep(0.5,3)*init.psi)**2)
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=save.results,save.graphs=save.results)

orange.fit<-saemix(orange.model,orange,options)

# FAIL: IIV on phi1 only, phi1 fixed (eg trying to estimate only variance of phi1)
# Probably a stupid example anyway

orange.model2<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(0,1,1),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=save.results,save.graphs=save.results)

# orange.fit2<-saemix(orange.model2,orange,options)

#####################################################################################
