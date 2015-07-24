####################################################################################
####				Tests for SaemixData				####
####################################################################################
verbose<-FALSE
setwd("/home/eco/work/monolix/rversion/current/libcov/covsaemix/testsDev/")

setGeneric(name="read.saemixData",
  def=function(object){standardGeneric("read.saemixData")}
)

setGeneric(name="showall",
  def=function(object){standardGeneric("showall")}
)

# SaemixData
source("../R/SaemixData.R")

cat("Beginning tests on SaemixData class...\n")

########################### Object creation
cat("            object creation:\n")

############ FAIL
##### Missing default values - FAIL with Error
# Creating an empty object: not possible, should print a warning
y<-capture.output(x<-saemixData())
if(x=="Creation of saemixData failed") cat("    creation of empty object: test OK\n") else cat("Test failed\n")

##### Wrong column names - FAIL
# Wrong Id name : "Creation of saemixData failed"
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="wrongid"))
if(x=="Creation of saemixData failed") cat("    wrong ID: test OK\n") else cat("Test failed\n")

# Wrong predictor name : "Creation of saemixData failed"
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.predictors="wrongpred"))
if(x=="Creation of saemixData failed") cat("    wrong predictor: test OK\n") else cat("Test failed\n")

# Wrong response name : "Creation of saemixData failed"
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.response="wrongresp"))
if(x=="Creation of saemixData failed") cat("    wrong response: test OK\n") else cat("Test failed\n")

# Wrong separator : "Creation of saemixData failed"
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T,sep=",", name.group=1,name.predictors=2,name.response=3))
if(x=="Creation of saemixData failed") cat("    wrong separator: test OK\n") else cat("Test failed\n")

############ SUCCESS
# No Id name but automatic recognition
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T,verbose=F))
if(class(x)=="SaemixData") cat("    automatic recognition: test OK\n") else cat("Test failed\n")

# Column names given as numbers
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group=1,name.predictors=2,name.response=3, name.covariates=4,units=list(x="mg",y="-",covariates="-"),verbose=F))
if(class(x)=="SaemixData") cat("    column name as numbers: test OK\n") else cat("Test failed\n")

# Two predictors
y<-capture.output(x<-saemixData(name.data="../data/theo.saemix.tab",header=T, name.predictors=c("Dose","Time"),verbose=F))
if(class(x)=="SaemixData" & length(x@name.predictors)==2) cat("    two predictors: test OK\n") else cat("Test failed\n")

# Two predictors, one does not exist - dropped with a warning
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.predictors=c("dose","wrongpred"),verbose=F))
if(class(x)=="SaemixData" & length(x@name.predictors)==1) cat("    two predictors, one removed: test OK\n") else cat("Test failed\n")

# Wrong number of units for covariates
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.covariates=4,units=list(x="mg",y="-",covariates=c("-","-")),verbose=F))
if (verbose) print(x)
if(class(x)=="SaemixData") cat("    too many units: test OK\n") else cat("Test failed\n")

# No units for covariates
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.covariates=4,units=list(x="mg",y="-"),verbose=F))
if (verbose) print(x)

# No units for x
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.covariates=4,units=list(y="-"),verbose=F))
if (verbose) print(x)
if(class(x)=="SaemixData") cat("    no units for x: test OK\n") else cat("Test failed\n")

# Bimodal response
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.response="gender",verbose=F))
if(class(x)=="SaemixData") cat("    warning for binary response: test OK\n") else cat("Test failed\n")

############ SUCCESS
##### Giving all object specifications
y<-capture.output(xok<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-"),verbose=F))
if (verbose) print(xok)
if(class(xok)=="SaemixData") cat("    creation of full object: test OK\n") else cat("Test failed\n")

# Creating an object, giving only the name - file
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",verbose=F))
if (verbose) print(x)
if(class(x)=="SaemixData") cat("    creation from file name: test OK\n") else cat("Test failed\n")

# Creating an object, giving only the name - data.frame
tab1<-read.table("../data/PD1.saemix.tab",header=T)
y<-capture.output(x<-saemixData(name.data="tab1",verbose=F))
if (verbose) print(x)
if(class(x)=="SaemixData") cat("    creation from dataframe name: test OK\n") else cat("Test failed\n")

# Creating an object, giving a data.frame 
y<-capture.output(x<-saemixData(name.data=tab1,verbose=F))
if(class(x)=="SaemixData") cat("    creation from dataframe: test OK\n") else cat("Test failed\n")

############ SUCCESS
##### Censored data
tab1<-read.table("../data/PD1.saemix.tab",header=T)
tab1<-cbind(tab1,censored=0)
tab1$censored[tab1$response<5]<-1
tab1$response[tab1$response<5]<-5
y<-capture.output(x<-saemixData(name.data=tab1,verbose=F,name.cens="censored"))
if(class(x)=="SaemixData" & x@name.cens=="cens") cat("    reading censored data: test OK\n") else cat("Test failed\n")

##### Several responses
tab1<-read.table("../data/PD1.saemix.tab",header=T)
tab2<-read.table("../data/PD2.saemix.tab",header=T)
tab3<-rbind(cbind(tab1,type.rep=1),cbind(tab2,type.rep=2))
tab3<-tab3[order(tab3$subject,tab3$dose),]
y<-capture.output(x<-saemixData(name.data=tab3,verbose=F,name.ytype="type.rep"))
if(class(x)=="SaemixData" & x@name.ytype=="ytype") cat("    reading multiple responses: test OK\n") else cat("Test failed\n")

##### Missing data
tab1<-read.table("../data/PD1.saemix.tab",header=T)
tab1<-cbind(tab1,censored=0)
tab1$censored[tab1$response<5]<-1
tab1$response[tab1$response<5]<-5
tab2<-read.table("../data/PD2.saemix.tab",header=T)
tab3<-rbind(cbind(tab1,type.rep=1,missing=0),cbind(tab2,censored=0,type.rep=2,missing=0))
tab3<-tab3[order(tab3$subject,tab3$dose),]
tab3$missing[c(1,4,7)]<-1
y<-capture.output(x<-saemixData(name.data=tab3,verbose=F,name.ytype="type.rep",name.cens="censored",name.mdv="missing"))
if(class(x)=="SaemixData" & x@name.mdv=="mdv") cat("    reading missing data: test OK\n") else cat("Test failed\n")

##### Several occasions
tab1<-read.table("../data/PD1.saemix.tab",header=T)
tab1<-cbind(tab1,occas=1,censored=0)
tab1$censored[tab1$response<5]<-1
tab1$response[tab1$response<5]<-5
tab1$occas[tab1$dose==90]<-2
tab2<-read.table("../data/PD2.saemix.tab",header=T)
tab3<-rbind(cbind(tab1,type.rep=1,missing=0),cbind(tab2,occas=1,censored=0,type.rep=2,missing=0))
tab3<-tab3[order(tab3$subject,tab3$dose),]
tab3$missing[c(1,4,7)]<-1
y<-capture.output(x<-saemixData(name.data=tab3,verbose=F,name.ytype="type.rep",name.cens="censored",name.mdv="missing",name.occ="occas"))
if(class(x)=="SaemixData" & x@name.occ=="occ") cat("    reading occasions: test OK\n") else cat("Test failed\n")

##### Genetics data
dirdat<-"/home/eco/work/monolix/rversion/data/simulgenetics/"
pkpddat<-read.table(file.path(dirdat,"pkpd_withcov_full.csv"),header=T,sep='\t')

y<-capture.output(x<-saemixData(name.data=pkpddat,verbose=F,name.ytype="ytype",name.genetic.covariates = paste("gen",1:9,sep="")))

y<-capture.output(x<-saemixData(name.data=pkpddat,verbose=F,name.ytype="ytype",name.covariates = c('wt',"age","trt"),name.genetic.covariates = paste("gen",1:9,sep="")))
cat("    reading genetic covariates: ")
if(class(x)=="SaemixData" & length(x@ind.gen)==12 & sum(x@ind.gen)==9) cat("test OK\n") else cat("test failed\n")

# x<-saemixData(name.data=pkpddat,verbose=F,name.ytype="ytype",name.covariates = c('wt',"age","trt"),name.genetic.covariates = paste("gen",1:9,sep=""))

##### Covariate transformation - TODO finish, test, use transformed covariates instead of ocov (trans.cov,...)


############ SUCCESS
##### Problems with dataset

# X with missing values - remove corresponding lines
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[15,2]<-dat[41,2]<-dat[55,2]<-NA
y<-capture.output(x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-"),verbose=F))
if (verbose) print(x)
if(class(x)=="SaemixData" & dim(x@data)[1]==297) cat("    removing missing predictor values: test OK\n") else cat("Test failed\n")

# Subjects with all missing values - remove entire subject
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[7:12,3]<-NA
y<-capture.output(x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-")))
if (verbose) print(x)
if(class(x)=="SaemixData" & dim(x@data)[1]==294) cat("    removing subject with only missing values: test OK\n") else cat("Test failed\n")

# Y with missing values - set corresponding mdv items to 1
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[15,3]<-dat[41,3]<-dat[55,3]<-dat[53,3]<-NA
y<-capture.output(x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-")))
if (verbose) print(x)
if(class(x)=="SaemixData" & sum(x@data$mdv)==4) cat("    setting missing values to 1: test OK\n") else cat("Test failed\n")

# covariate with missing values
dat<-read.table("../data/PD1.saemix.tab",header=T)
dat[15,4]<-dat[41,4]<-dat[55,4]<-NA
y<-capture.output(x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-")))
if(verbose) print(summary(x))
if(class(x)=="SaemixData" & sum(is.na(x@data$gender))==3) cat("    missing covariates (no warning): test OK\n") else cat("Test failed\n")

# Covariate with character values
dat<-read.table("../data/PD1.saemix.tab",header=T)
vec<-ifelse(dat$gender==1,"W","M")
dat$gender<-vec
y<-capture.output(x<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-")))
if (verbose) print(x)
if(class(x)=="SaemixData" & typeof(x@data$gender)=="double") cat("    character type covariates recoded: test OK\n") else cat("Test failed\n")

# Covariate with factor values
dat<-read.table("../data/PD1.saemix.tab",header=T)
vec<-ifelse(dat$gender==1,"W","M")
dat$gender<-as.factor(vec)
y<-capture.output(xchar<-saemixData(name.data=dat,name.group="subject",name.predictors="dose", name.response="response",name.covariates="gender", units=list(x="mg",y="-",covariates="-")))
if (verbose) print(xchar)
if(verbose) print(table(xchar@data$gender))
#cat("Warning, a covariate model when the covariate has character/factor values will probably fail\n")
if(class(x)=="SaemixData" & typeof(x@data$gender)=="double") cat("    factor type covariates recoded: test OK\n") else cat("Test failed\n")

# 2 covariates
y<-capture.output(x<-saemixData(name.data="../data/theo.saemix.tab",name.group="Id",name.predictors=c("Time","Dose"), name.response="Concentration",name.covariates=c("Weight","Sex"), units=list(x="mg",y="-",covariates=c("kg","-"))))
if (verbose) print(x)
if(class(x)=="SaemixData" & typeof(x@data$Sex)=="double") cat("    character type covariates recoded: test OK\n") else cat("Test failed\n")

# Error on a covariate name
y<-capture.output(x<-saemixData(name.data="../data/theo.saemix.tab",name.group="Id",name.predictors=c("Time","Dose"), name.response="Concentration",name.covariates=c("Height","Sex"), units=list(x="mg",y="-",covariates=c("kg","-"))))
if (verbose) print(x)
if(class(x)=="SaemixData" & identical(x@name.covariates,c("Sex"))) cat("    error in one covariate name: test OK\n") else cat("Test failed\n")

cok<-readline(prompt="Finished testing creation of object, press any key\n")

########################### Methods

############ Covariate transformation
##### Continuous covariates
cat("            method transform:\n")
dat<-read.table("../data/cow.saemix.tab",header=T)
vec<-ifelse(dat$twin==1,"Single","Twin")
dat$twin<-vec
vec<-as.factor(dat$birthrank)
dat$birthrank<-as.factor(vec)
y<-capture.output(x<-saemixData(name.data=dat,name.group="cow",name.predictors=c("time"), name.response="weight",name.covariates=c("birthyear","twin","birthrank"), units=list(x="d",y="kg",covariates=c("yr","-","-"))))

x2<-transform.cov(x,birthyear,verbose=T)
x2<-transform.cov(x,birthyear,centering="mean",verbose=T)
x2<-transform.cov(x,birthyear,centering=1998,verbose=T)

x2<-transform.cov(x,birthyear,centering=1900,transformation=function(x) log(x),verbose=T)

if (verbose) print(xchar)
if(verbose) print(table(xchar@data$gender))
#cat("Warning, a covariate model when the covariate has character/factor values will probably fail\n")
if(class(x)=="SaemixData" & typeof(x@data$gender)=="double") cat("    factor type covariates recoded: test OK\n") else cat("Test failed\n")

##### Categorical covariates

if(verbose) print(table(x@data$Sex))

############ SUCCESS
##### Method show
cat("            method show:\n")

y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-")))
show(x)

showall(x)

##### Method subset
subset(x,gender==1 & dose==0)

print(summary(subset(x,subject<3)))

##### Method plot

cat("            method plot:\n")
y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-")))

plot(x)

plot(x,main="Raw data",type="p",col="Blue")

plot(x,main="Raw data",type="b")

plot(x,main="Raw data",type="b",col="Red",ylog=TRUE)

# Individual plots
plot(x,type="b",individual=TRUE,col="Blue",nmax=30)

plot(x,type="b",individual=TRUE,pch=3,lty=2,ilist=1:4)

plot(x,type="b",individual=TRUE,col="Blue",ilist=1:24,limit=FALSE)

plot(x,type="l",individual=TRUE,col="Red",sample=T)

cat("End of tests on SaemixData class.\n")

####################################################################################
####				Tests for SaemixRepData				####
####################################################################################

cat("Beginning tests on SaemixRepData class...\n")

y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-")))

xrep<-new(Class="SaemixRepData",data=x)
show(xrep)

cat("End of tests on SaemixRepData class.\n")

####################################################################################
####				Tests for SaemixSimData				####
####################################################################################

cat("Beginning tests on SaemixSimData class...\n")

y<-capture.output(x<-saemixData(name.data="../data/PD1.saemix.tab",header=T, name.group="subject",name.predictors="dose",name.response="response", name.covariates="gender",units=list(x="mg",y="-",covariates="-")))

y<-capture.output(xsim<-new(Class="SaemixSimData",data=x))
show(xsim)

cat("End of tests on SaemixSimData class.\n")

####################################################################################
####					Mirror plots				####
####################################################################################

# ECO TODO
cat("Beginning second series of tests on SaemixSimData class, filling the data...\n")
cat("ECO TODO \n")

y<-capture.output(x<-saemixData(name.data="../data/theo.saemix.tab",header=T, name.predictors=c("Dose","Time"),name.X="Time"))

xsim<-new(Class="SaemixSimData",data=x)
show(xsim)

#cat("End of mirror plots.\n")

####################################################################################
