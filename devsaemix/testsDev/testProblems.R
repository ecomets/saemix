
###################################################################################### Example with sigmoid data and some fit issues (to be investigated further)

# Monolix: 
# fails when a simple Emax model is used (doesn't start)
# works with a Hill model, estimating gamma~5, omega(Emax) estimated very large
# LL by Monolix: 3585.21 (IS) [3577.49 (IS) for the model with S0] vs here 3406.60 => why 200 points difference (maybe due to obs~2000 ?)

setwd("/home/eco/work/monolix/rversion/current/libcov")
dat.oger<-saemixData(name.data="../../data/oger/DB.csv",name.group="OBS", name.predictors="WEEK", name.response="SIZE",name.covariates=c("TKI","ARA2"), units=list(x="wk",y="-",covariates=c("-","-")),sep=",")

modelemax<-function(psi,id,xidep) {
	# input:
	#   psi : matrix of parameters (3 columns, E0, Emax, EC50)
	#   id : vector of indices 
	#   xidep : dependent variables (same nb of rows as length of id)
	# returns:
	#   a vector of predictions of length equal to length of id
	week<-xidep[,1]
	e0<-psi[id,1]
	emax<-psi[id,2]
	e50<-psi[id,3]
	f<-e0+emax*week/(e50+week)
	return(f)
}

modelesigm.nobase<-function(psi,id,xidep) {
	# input:
	#   psi : matrix of parameters (3 columns, E0, Emax, EC50)
	#   id : vector of indices 
	#   xidep : dependent variables (same nb of rows as length of id)
	# returns:
	#   a vector of predictions of length equal to length of id
	week<-xidep[,1]
	emax<-psi[id,1]
	e50<-psi[id,2]
	gamma<-psi[id,3]
	f<-emax*(week**gamma)/((e50**gamma)+(week**gamma))
	return(f)
}

omod1<-saemixModel(model=modelemax,description="Emax model", psi0=matrix(c(10,2000,5,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined")

# Pb when estimating the individual parameters
saemix.options1<-list(directory="betatest/ogeremx.comb",save=save.results,map=F,lls.is=F,fim=F)
fit.oger1<-saemix(omod1,dat.oger,saemix.options1)

omod1<-saemixModel(model=modelemax,description="Emax model", psi0=matrix(c(10,2000,7,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(0,0,0), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined")

saemix.options1<-list(directory="betatest/ogeremx.comb",save=save.results,map=F,lls.is=F,fim=F)
fit.oger1<-saemix(omod1,dat.oger,saemix.options1)

plot(fit.oger1,plot.type="individual")

omod2<-saemixModel(model=modelemax,description="Emax model", psi0=matrix(c(10,2000,5,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("E0","Emax","EC50"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE))

saemix.options2<-list(directory="betatest/ogeremx.const",save=save.results)
fit.oger2<-saemix(omod2,dat.oger,saemix.options2)
plot(fit.oger2,plot.type="individual")

omod3<-saemixModel(model=modelesigm.nobase,description="Sigmoid model without baseline", psi0=matrix(c(2000,2,2,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("Emax","EC50","gamma"))), transform.par=c(1,1,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined")

saemix.options3<-list(directory="betatest/ogersigm.comb")
fit.oger3<-saemix(omod3,dat.oger,saemix.options3)

plot(fit.oger3,plot.type="individual",smooth=T)

omod4<-saemixModel(model=modelesigm.nobase,description="Sigmoid model without baseline", psi0=matrix(c(2000,2,4,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("Emax","EC50","gamma"))), transform.par=c(0,0,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(1,1,0,1,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined")

saemix.options4<-list(directory="betatest/ogercorr.comb")
fit.oger4<-saemix(omod4,dat.oger,saemix.options4)

plot(fit.oger4,plot.type="individual",smooth=T)

# No IIV on Emax - best run so far and SE obtained
omod5<-saemixModel(model=modelesigm.nobase,description="Sigmoid model without baseline", psi0=matrix(c(2000,2,4,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("Emax","EC50","gamma"))), transform.par=c(0,0,1), covariate.model=matrix(c(0,0,0,0,0,0),ncol=3,byrow=TRUE), fixed.estim=c(1,1,1), covariance.model=matrix(c(0,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="combined")

saemix.options5<-list(directory="betatest/ogerIIV1.comb")
fit.oger5<-saemix(omod5,dat.oger,saemix.options5)

plot(fit.oger5,plot.type="individual",smooth=T)

#####################################################################################
#### Nathalie

setwd("/home/eco/work/monolix/rversion/examples/nathalie")

timoun.data<-saemixData(name.data="timounclean_taille2.csv",header=TRUE, sep=",", name.group=c("NUM"),name.predictors=c("age"), name.response=c("weight"),name.covariates=c("length","rang"), units=list(x="sem",y="kg",covariates=c("cm","-")))

dat<-read.table("timounclean_taille2.csv",sep=",",header=T)
dat<-dat[!is.na(dat$length),]

log.mod<-function(psi,id,xidep) {
	age<-xidep[,1]
	int<-psi[id,1]
	sl1<-psi[id,2]
	sl2<-psi[id,3]
	return(int+sl1*age+sl2*log(age))
}

model5.log<-saemixModel(model=log.mod,description="Linear+log-linear function of age", psi0=matrix(c(3000,30,1000,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL,c("int","sl.lin","sl.log"))), covariate.model=matrix(c(1,1,1,1,1,0,0,0,0),ncol=3,byrow=TRUE), transform.par=c(0,0,0))

#####################################################################################
#### 

cd4<-read.table("/home/eco/cours/tps/saemix/forTP/cd4_MACS_TP.tab",header=T)

# Création de l'objet données
smx.dat<-saemixData(cd4,name.group="ID",name.predictors="tim.sero", name.response="CD4",name.covariates=c("age","cesd","nb.smoke","nb.partners"), units=list(x="yr",y="counts",covariates=c("yr","-","-","-")))

expmod<-function(psi,id,xidep) { 
	time<-xidep[,1]  
	E0<-psi[id,1]
	k<-psi[id,2]
	alp<-psi[id,3]
	ypred<-E0*(alp+(1-alp)*exp(-k*time))
	ypred[time<0]<-E0[time<0]
	return(ypred)
}

smx.model2<-saemixModel(model=expmod,
												description="Constant+exponential model",
												psi0=matrix(c(1000,0.8,0.2,0,0,0),ncol=3,byrow=TRUE, 
																		dimnames=list(NULL, c("E0","k","alpha"))),transform.par=c(1,1,1))

smx.options<-list(seed=89823,directory="betatest/diggleCD4/exp1")

fit.exp<-saemix(smx.model2,smx.dat,smx.options)

plot(fit.exp,plot.type="vpc")

# psi0 with only one line

smx.model3<-saemixModel(model=expmod,description="Constant+exponential model", psi0=matrix(c(1000,0.8,0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("E0","k","alpha"))),transform.par=c(1,1,1))
smx.options<-list(seed=89823,save=save.results,save.graphs=save.results)
fit.exp<-saemix(smx.model3,smx.dat,smx.options)

smx.model3<-saemixModel(model=expmod,description="Constant+exponential model", psi0=matrix(c(1000,0.8,0.2),ncol=3,byrow=TRUE,dimnames=list(NULL, c("E0","k","alpha"))),transform.par=c(1,1,1),covariate.model=matrix(c(1,rep(0,11)),ncol=3,byrow=TRUE))
smx.options<-list(seed=89823,save=save.results,save.graphs=save.results)
fit.exp<-saemix(smx.model3,smx.dat,smx.options)

smx.model3<-saemixModel(model=expmod,description="Constant+exponential model", psi0=matrix(c(1000,0.8,0.2,rep(0,3)),ncol=3,byrow=TRUE,dimnames=list(NULL, c("E0","k","alpha"))),transform.par=c(1,1,1),covariate.model=matrix(c(1,rep(0,11)),ncol=3,byrow=TRUE))
smx.options<-list(seed=89823,save=save.results,save.graphs=save.results)
fit.exp<-saemix(smx.model3,smx.dat,smx.options)
