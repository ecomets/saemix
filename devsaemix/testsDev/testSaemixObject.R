####################################################################################
####				Tests for SaemixObject				####
####################################################################################

############################
# Defining generic for new methods
setGeneric(name="read.saemixData",
  def=function(object){standardGeneric("read.saemixData")}
)

setGeneric(name="showall",
  def=function(object){standardGeneric("showall")}
)

setGeneric(name="psi",
  def=function(object,indiv.par){standardGeneric("psi")}
)

setGeneric(name="phi",
  def=function(object,indiv.par){standardGeneric("phi")}
)

setGeneric(name="eta",
  def=function(object,indiv.par){standardGeneric("eta")}
)

############################
# SaemixData
source("../R/SaemixData.R")

# SaemixRes
source("../R/SaemixRes.R")

# SaemixModel
source("../R/SaemixModel.R")

# SaemixObject
source("../R/SaemixObject.R")
source("../R/func_plots.R")

####################################################################################
####				Tests for SaemixObject				####
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

############ Model without covariate, SUCCESS

##### Data
xdat<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"))
print(xdat)

##### Model
ymod<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

##### Object
yobj<-new(Class="SaemixObject",data=xdat,model=ymod)
print(yobj)

############ Model with covariate, SUCCESS

dat<-read.table("../data/theo.saemix.tab",header=T)

##### Data
xdat2<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))
print(xdat)

##### Model
ymod2<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

##### Object

yobj2<-new(Class="SaemixObject",data=xdat2,model=ymod2)
print(yobj2)

############ Model with covariates and missing values, SUCCESS

dat<-read.table("../data/theo.saemix.tab",header=T)
dat$Weight[dat$Id==3]<-NA

##### Data
xdat2<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))
print(xdat)

##### Model
ymod2<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

##### Object

yobj2<-new(Class="SaemixObject",data=xdat2,model=ymod2)
print(yobj2)


############ Model with covariates not in the dataset, SUCCESS

dat<-read.table("../data/theo.saemix.tab",header=T)

##### Data
xdat3<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))
print(xdat)

##### Model
ymod3<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1,1,1,1),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

##### Object

yobj3<-new(Class="SaemixObject",data=xdat3,model=ymod3)
print(yobj3)

############ Model with covariates not in the dataset, SUCCESS

##### Data
xdat4<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))
print(xdat)

##### Model
ymod4<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1,1,1,1),ncol=3,byrow=TRUE), fixed.estim=c(1,0,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

##### Object

yobj4<-new(Class="SaemixObject",data=xdat4,model=ymod4)
print(yobj4)


############ TODO: covariables as factor

dat<-read.table("../data/theo.saemix.tab",header=T)
vec<-ifelse(dat$Sex==1,"F","M")
dat$Sex<-vec
dat$Weight[dat$Id==3]<-NA

##### Data
xdat.cat<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))


####################################################################################
