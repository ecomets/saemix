####################################################################################
####				Tests of the main function			####
####################################################################################

setwd("/home/eco/work/monolix/rversion/libcov/covsaemix/testsDev/")
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
source("../R/func_aux.R")
source("../R/func_plots.R")
source("../R/func_main.R")

####################################################################################
####				Debugging	- orange trees		####
####################################################################################

# Pb reading the data - solved
orange<-saemixData(name.data=Orange,name.group="Tree",name.predictors="age", name.response="circumference",units=list(x="Days",y="mm"))

# Only one eta in the model, pb with cholesky (why ??) - solved
if(FALSE) {
	or1<-as.data.frame(Orange)
	or1[,1]<-as.integer(as.character(or1[,1]))
	orange<-saemixData(name.data=or1,name.group="Tree",name.predictors="age", name.response="circumference",units=list(x="Days",y="mm"))
}

logistic.model<-function(psi,id,xidep) { 
	age<-xidep[,1]  
	phi1<-psi[id,1]
	phi2<-psi[id,2]
	phi3<-psi[id,3]
	resp<-phi1/(1+exp(-(age-phi2)/phi3))
	return(resp) 
}

############################## Model object creation errors
init.psi<-c(200,800,400)
init.omega<-diag(c(rep(0.5,3)*init.psi)**2)

# no IIV - should raise an error
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(0,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

# IIV on phi1 only, phi1 fixed - should raise an error
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(0,1,1),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

############################## Different combinations of IIV/fixed parameters
# Full IIV
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(200,800,400),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))))

options<-list(seed=94352514,save=FALSE,save.graphs=FALSE)

orange.fit<-saemix(orange.model,orange,options)

# IIV on phi1 only
init.psi<-c(200,800,400)
init.omega<-diag(c(rep(0.5,3)*init.psi)**2)
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit<-saemix(orange.model,orange,options)


# IIV on phi3 only
init.psi<-c(200,800,400)
init.omega<-diag(c(rep(0.5,3)*init.psi)**2)
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(0,0,0,0,0,0,0,0,1),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit.test1<-saemix(orange.model,orange,options)

# IIV on phi1 only, and phi2 fixed
orange.model2<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(1,0,1),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit2<-saemix(orange.model2,orange,options)

# IIV on phi3 only, and phi2 fixed
orange.model2<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(1,0,1),covariance.model=matrix(c(0,0,0,0,0,0,0,0,1),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit.test2<-saemix(orange.model2,orange,options)

# IIV on phi1 only, phi2 and phi3 fixed
orange.model3<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(1,0,0),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit3<-saemix(orange.model3,orange,options)


# IIV on phi1 only, and phi1 fixed - fails [normal, no more IIV]
orange.model2<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(0,1,1),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit2<-saemix(orange.model2,orange,options)

# IIV on phi1 & phi2, and phi1 fixed - works
orange.model4<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(0,1,1),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit4<-saemix(orange.model4,orange,options)


# Different CI
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(10,1200,100),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit<-saemix(orange.model,orange,options)

# Different CI, covariance matrix
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(20,1000,100),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=matrix(c(1000,0,0,0,1000,0,0,0,10),ncol=3,byrow=3))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit<-saemix(orange.model,orange,options)

# LN distribution
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(10,1200,100),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),transform.par=c(1,1,1))

options<-list(seed=943514,save=FALSE,save.graphs=FALSE)

orange.fit<-saemix(orange.model,orange,options)

# IIV on 2 parameters
orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(20,1000,100),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,1,0,0,0,0),ncol=3,byrow=3))

options<-list(seed=94352514,save=FALSE,save.graphs=FALSE)

orange.fit<-saemix(orange.model,orange,options)

# Only one parameter to estimate : fails

logistic.par1<-function(psi,id,xidep) { 
	age<-xidep[,1]  
	phi1<-psi[id,1]
	phi2<-724
	phi3<-345
	resp<-phi1/(1+exp(-(age-phi2)/phi3))
	return(resp) 
}
orange.par1<-saemixModel(model=logistic.par1,description="Logistic growth with only 1 paramater", psi0=matrix(c(200),ncol=1,byrow=TRUE,dimnames=list(NULL,c("phi1"))))

options<-list(seed=32594514,save=FALSE,save.graphs=FALSE)

orange.fit<-saemix(orange.par1,orange,options)

########### Debugging code

saemix.model3<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(200,800,400),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3))
saemix.model3@omega.init

saemix.model4<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(200,800,400),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),transform.par=c(0,1,1),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3))
saemix.model4@omega.init
orange.fit<-saemix(saemix.model4,orange,options)

###########################
# predict function

predict(orange.fit)

###########################
########### Debugging code
saemix.data2<-orange
saemix.model2<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(200,800,400),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3))
saemix.model2<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(200,800,400),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=diag(c(4000,4000,4000)))
saemix.options2<-list(seed=943514,save=FALSE,save.graphs=FALSE)

# Running only first step
step1<-TRUE
step2<-FALSE
step3<-FALSE
step4.SA<-FALSE
step4.max<-FALSE
step4.RS<-FALSE
step5<-FALSE
niter<-1

# Running all the algorithm
step1<-step2<-step3<-step4.SA<-step4.max<-step4.RS<-step5<-TRUE
niter<-15

source("code_begin.R")

for(kiter in 1:niter) {
	source("code_iterations.R")	
}

####################################################################################
# Directory issue

system("rm -r mytestdir")

orange.model<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(c(200,800,400),ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))))
options<-list(seed=94352514,save=TRUE,save.graphs=TRUE,directory="mytestdir")

orange.fit<-saemix(orange.model,orange,options)
system("ls -ltr mytestdir/")

orange.fit<-saemix(orange.model,orange,options)
system("ls -ltr mytestdir/")

# Limiting the number of subjects in the plot
plot(orange.fit,plot.type="population.fit",ilist=1:12)
plot(orange.fit,plot.type="individual.fit",ilist=1:12,smooth=T)

####################################################################################
####				Debugging			####
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

##### Data, missing value in covariate
dat<-read.table("../data/theo.saemix.tab",header=T)
dat$Weight[dat$Id==3]<-NA

xdat<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))

##### Model with combined residual error
ymod<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

ymod1<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

sopt<-list(algorithms=c(1,1,1),save.graphs=FALSE,save=FALSE)

sopt<-list(algorithms=c(1,1,1))
sfit<-saemix(ymod,xdat,sopt)
sfit1<-saemix(ymod1,xdat,sopt)

plot(sfit,plot.type="parameters.vs.covariates")

##### Data, covariate with character value
dat<-read.table("../data/theo.saemix.tab",header=T)
dat$Weight[dat$Id==3]<-NA
vec<-dat$Sex
dat$Sex<-ifelse(vec==1,"M","F")

xdat<-saemixData(name.data=dat,header=T, name.group="Id", name.predictors=c("Dose","Time"),name.response="Concentration", name.covariates=c("Sex","Weight"),name.X="Time", units=list(x="",y="-",covariates=c("-","kg")))

##### Model with combined residual error
ymod<-saemixModel(model=model1cpt, description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,1),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

sopt<-list(algorithms=c(1,1,1))
sfit<-saemix(ymod,xdat,sopt)

# Subsets
v1<-subset(xdat,Id<3)
is.factor(v1@data$Weight)
is.factor(v1@data$Sex)
summary(v1@data)

####################################################################################
# Nathalie

setwd("/home/eco/work/monolix/rversion/examples/nathalie")

timoun.data<-saemixData(name.data="timounclean_taille2.csv",header=TRUE, sep=",", name.group=c("NUM"),name.predictors=c("age"), name.response=c("weight"),name.covariates=c("length","rang"), units=list(x="sem",y="kg",covariates=c("cm","-")))

dat<-read.table("timounclean_taille2.csv",sep=",",header=T)
dat<-dat[!is.na(dat$length),]

timoun.data.cov<-saemixData(name.data=dat,header=TRUE, sep=",",name.group=c("NUM"),name.predictors=c("age"), name.response=c("weight"),name.covariates=c("length","rang"), units=list(x="sem",y="kg",covariates=c("cm","-")))

plot(timoun.data,type="p")

log.mod<-function(psi,id,xidep) {
	age<-xidep[,1]
	int<-psi[id,1]
	sl1<-psi[id,2]
	sl2<-psi[id,3]
	return(int+sl1*age+sl2*log(age))
}
model1.log<-saemixModel(model=log.mod,description="Linear+log-linear function of age", psi0=matrix(c(3000,30,1000,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("int","sl.lin","sl.log"))), transform.par=c(0,0,0))

opt1<-list(seed=546532432,directory="fits/model1.log")

fit1.log<-saemix(model1.log,timoun.data,opt1)

# 2 lignes de cov, 2 cov (ok)
model5.log<-saemixModel(model=log.mod,description="Linear+log-linear function of age", psi0=matrix(c(3000,30,1000,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("int","sl.lin","sl.log"))), covariate.model=matrix(c(1,1,0,0,0,1),ncol=3,byrow=TRUE), transform.par=c(0,0,0))

opt1<-list(fix.seed=FALSE,directory="fits/model5.log")

fit5.log<-saemix(model5.log,timoun.data.cov,opt1)

# 3 lignes de cov, 2 cov => prend les 2 premières (ok)
model5.log.1<-saemixModel(model=log.mod,description="Linear+log-linear function of age", psi0=matrix(c(3000,30,1000,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("int","sl.lin","sl.log"))), covariate.model=matrix(c(1,1,1,1,1,0,0,0,0),ncol=3,byrow=TRUE), transform.par=c(0,0,0))

opt1<-list(fix.seed=FALSE,save=FALSE)

fit5.log.1<-saemix(model5.log.1,timoun.data.cov,opt1)

# 1 ligne de cov, 2 cov => prend les 2 premières (ok)
model5.log.2<-saemixModel(model=log.mod,description="Linear+log-linear function of age", psi0=matrix(c(3000,30,1000,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("int","sl.lin","sl.log"))), covariate.model=matrix(c(1,0,1),ncol=3,byrow=TRUE), transform.par=c(0,0,0))

opt1<-list(fix.seed=FALSE,save=FALSE)

fit5.log.2<-saemix(model5.log.2,timoun.data.cov,opt1)

# Testing it resaves if the directory exists
model1.log<-saemixModel(model=log.mod,description="Linear+log-linear function of age", psi0=matrix(c(3000,30,1000,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("int","sl.lin","sl.log"))), transform.par=c(0,0,0))

opt1<-list(seed=546532432,directory="fits/model1.log")

fit1.log<-saemix(model1.log,timoun.data,opt1)
system("ls -ltr fits/model1.log/")
system("date")

####################################################################################
# Parkinson, pb de subset

rpark<-read.table("/home/eco/cours/tps/saemix/parkinson/parkinson_subset.csv", sep=",", header=T,na=".")
rpark$SEX<-factor(rpark$SEX,labels=c("F","M"))

rsmx<-saemixData(name.data=rpark,name.group=c("SID"),name.predictors=c("time"), name.response=c("UPDR"),name.covariates=c("AGE","WT","SEX"), units=list(x="yr",y="score", covariates=c("yr","kg","-")))	

v1<-subset(rsmx,SID==12013)

v1<-subset(rsmx,(SID==13712 |SID==12013))

# Read data
park<-read.csv("/home/eco/cours/tps/saemix/parkinson/parkinson_full.csv",header=T,na=".")
park[,4]<-park[,4]/365
colnames(park)[4]<-"time"
park<-park[order(park[,2]),]

# Subset with DP=0 & LD=0
park.plac<-subset(park,LD==0 & DP==0)
ndat<-tapply(park.plac[,1],park.plac[,1],length)

park.plac<-subset(park,LD==0 & DP==0 & TC==0 & BC==0 & PG==0)
ndat<-tapply(park.plac[,2],park.plac[,2],length)
table(ndat)

####################################################################################
# Elodie, trying different distributions for a parameter without IIV

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
namdir.data<-file.path("/home/eco/work/monolix/rversion/simultests/elodie","R2A")
prefix.data<-"mc-sim-1-"
irep<-1
dat<-read.table(file.path(namdir.data,paste(prefix.data,irep,".dat",sep="")),header=T,skip=1)
dat<-cbind(dat,treat=ifelse(dat$ID<50,0,1))
xdat<-saemixData(name.data=dat,header=T,name.group="ID",name.predictors="DOSE",name.response="DV",name.covariates=c("treat"),units=list(x="mg",y="-",covariates="-"))			

pvrai<-c(30,500,2,5,0.490,0.245,0.490,0.090,4)
pfaux<-c(60,1000,1,10,0.100,0.01,0.1,0.1,1)
nsim<-100
ci.param<-pvrai
error.init<-c(sqrt(ci.param[9]),0)

elodie.mdist<-saemixModel(model=modelhill,description="Emax model with sigmoidicity (Hill model)",psi0=matrix(c(ci.param[c(4,1:2,3)],0,0,0,0), ncol=4,byrow=T, dimnames=list(NULL,c("E0","Emax","EC50","gamma"))), transform.par=c(1,1,1,0), covariate.model=matrix(c(0,0,0,0),ncol=4,byrow=T), covariance.model=matrix(c(1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0),ncol=4,byrow=T),omega.init=diag(c(ci.param[c(8,5,7)],1)),error.init=error.init)

elodie.ldist<-saemixModel(model=modelhill,description="Emax model with sigmoidicity (Hill model)",psi0=matrix(c(ci.param[c(4,1:2,3)],0,0,0,0), ncol=4,byrow=T, dimnames=list(NULL,c("E0","Emax","EC50","gamma"))), transform.par=c(1,1,1,1), covariate.model=matrix(c(0,0,0,0),ncol=4,byrow=T), covariance.model=matrix(c(1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0),ncol=4,byrow=T),omega.init=diag(c(ci.param[c(8,5,7)],1)),error.init=error.init)

opt1<-list(seed=76765767,save=FALSE,save.graphs=FALSE)

fit1<-saemix(elodie.mdist,xdat,opt1)
fit2<-saemix(elodie.ldist,xdat,opt1)

####################################################################################
##### bug Julie

# Dummy dataset

tim<-seq(1,24,2)
nsuj<-30
par<-matrix(rnorm(nsuj*3),ncol=3)
par[,1]<-2*exp(par[,1]*.3)
par[,2]<-10*exp(par[,2]*.3)
par[,3]<-1*exp(par[,3]*.3)

xpred<-NULL
for(isuj in 1:nsuj) {
	ka<-par[isuj,1]
	Vc<-par[isuj,2]
	k<-par[isuj,3]/Vc
	yp<-100/Vc*ka/(ka-k)*(exp(-k*tim)/(1-exp(-k*12))-exp(-ka*tim)/(1-exp(-ka*12)))
	yp<-yp*(1+rnorm(length(yp))*.2)
	xtab<-data.frame(id=rep(isuj,length(yp)),time=xtim,conc=yp)
	xpred<-rbind(xpred,xtab)
}
plot(xpred$time,xpred$conc)

xpred<-cbind(xpred,dose=rep(100,dim(xpred)[1]))
saemix.data.e<-saemixData(name.data=xpred,name.group="id",name.predictors=c("dose","time"), name.response="conc",name.X="time",units=list(x="hr",y="mg/L"))

model1cpt<-function(psi,id,xidep) {
	dose<-xidep[,1]
	ka<-psi[id,1]  
	Vc<-psi[id,2]
	Cl<-psi[id,3]
	k<-Cl/Vc
	tim<-xidep[,2]
	ypred<-dose/Vc*ka/(ka-k)*(exp(-k*tim)/(1-exp(-k*12))-exp(-ka*tim)/(1-exp(-ka*12)))
	return(ypred)
}

cov.matrix<-diag(1,3)

saemix.model.e<-saemixModel(model=model1cpt, description="One-compartment model oral 1",psi0=matrix(c(1,100,1),ncol=3, byrow=TRUE,dimnames=list(NULL,c("ka","Vc","Cl"))), transform.par=c(1,1,1), fixed.estim=c(1,1,1), covariance.model=cov.matrix, omega.init=cov.matrix, error.model="constant")

sopt<-list(seed=632545,nb.chains=15, nbiter.saemix = c(600, 200))
NVP.fit<-saemix(saemix.model.e,saemix.data.e,sopt)

sopt<-list(seed=632545,nb.chains=3)
NVP.fit<-saemix(saemix.model.e,saemix.data.e,sopt)

# Pb: still estimates the 3 phi !!! => debug TODO
cov.matrix<-diag(0,3)
cov.matrix[3,3]<-1

saemix.model.e<-saemixModel(model=model1cpt, description="One-compartment model oral 1",psi0=matrix(c(1,100,1),ncol=3, byrow=TRUE,dimnames=list(NULL,c("ka","Vc","Cl"))), transform.par=c(1,1,1), fixed.estim=c(0,1,1), covariance.model=cov.matrix, omega.init=cov.matrix, error.model="constant")

################################
# Test pas à pas
inp.data<-saemix.data.e
inp.model<-saemix.model.e
inp.options<-list(seed=632545,nb.chains=3)

source("debug_begin.pr")


# Running only first step
saemix.data2<-saemix.data.e
saemix.model2<-saemix.model.e
saemix.options2<-sopt

step1<-TRUE
step2<-FALSE
step3<-FALSE
step4.SA<-FALSE
step4.max<-FALSE
step4.RS<-FALSE
step5<-FALSE
niter<-1

# Running all the algorithm
step1<-step2<-step3<-step4.SA<-step4.max<-step4.RS<-step5<-TRUE
niter<-200

source("code_begin.R")

for(kiter in 1:niter) {
	source("code_iterations.R")	
}

####################################################################################
# Logit transformation


################################
# Test pas à pas
inp.data<-xdat
inp.model<-ymod
inp.options<-sopt

source("debug_begin.pr")

####################################################################################
# Orange tree

xdat<-saemixData(name.data=Orange,name.group="Tree",name.predictors="age", name.response="circumference",units=list(x="Days",y="mm"))
ymod<-saemixModel(model=logistic.model,description="Logistic growth", psi0=matrix(init.psi,ncol=3,byrow=TRUE,dimnames=list(NULL,c("phi1","phi2","phi3"))),fixed.estim=c(0,1,1),covariance.model=matrix(c(1,0,0,0,0,0,0,0,0),ncol=3,byrow=3),omega.init=init.omega,error.init=c(30,0))

logistic.model<-function(psi,id,xidep) { 
	age<-xidep[,1]  
	phi1<-psi[id,1]
	phi2<-psi[id,2]
	phi3<-psi[id,3]
	resp<-phi1/(1+exp(-(age-phi2)/phi3))
	return(resp) 
}

sopt<-list(seed=943514,save=FALSE,save.graphs=FALSE)

################################
# Test pas à pas
inp.data<-xdat
inp.model<-ymod
inp.options<-sopt

source("debug_begin.pr")

####################################################################################
