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
source("../R/func_aux.R")
source("../R/func_plots.R")
source("../R/func_main.R")

############################################
# Algorithm initialisation
############################################

#### Theophylline example, with covariates
theo.saemix<-read.table("../data/theo.saemix.tab",header=T,na=".")
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"),name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

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
saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption",psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),covariate.model=matrix(c(0,0,1,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),covariance.model=matrix(c(1,0,0,0,1,1,0,1,1),ncol=3,byrow=TRUE),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="combined")

saemix.options<-list(seed=39546,save=FALSE,save.graphs=FALSE)

############################################
#  Main Algorithm
############################################
saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
saemix.options<-saemixObject["options"]
saemix.model<-saemixObject["model"]
saemix.data<-saemixObject["data"]
saemix.data@data<-saemix.data@data[saemix.data@data[,"mdv"]==0,]
# Initialising random generator
set.seed(saemix.options$seed)

# Initialisation - creating several lists with necessary information extracted (Uargs, Dargs, opt,varList, suffStat)
xinit<-initialiseMainAlgo(saemix.data,saemix.model,saemix.options)
Dargs<-xinit$Dargs
Uargs<-xinit$Uargs
varList<-xinit$varList
phiM<-xinit$phiM
mean.phi<-xinit$mean.phi
DYF<-xinit$DYF
opt<-xinit$opt
betas<-betas.ini<-xinit$betas
fixed.psi<-xinit$fixedpsi.ini
var.eta<-varList$diag.omega
theta0<-c(fixed.psi,var.eta[Uargs$i1.omega2],varList$pres[Uargs$ind.res])
