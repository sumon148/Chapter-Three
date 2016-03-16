# ------------------------------------------------------------------------------
# Script 1: Chapter Three
# Subject: Simulation study using ELL.2L, ELL.3L, EBP, & MQ Estimators 
# using Bangladesh data: Census 2001 and HIES 2000
# Description: Homoskedastic (HM) level-one and level-two errors
# ------------------------------------------------------------------------------
rm(list=ls(all=TRUE))
# ------------------------------------------------------------------------------
# Loading Packages
# ------------------------------------------------------------------------------
library(samplingbook)
library(nlme)
library(mvtnorm)
library(foreach)
library(doMC)
library(MASS)
library(np)
registerDoMC(20)
#-------------------------------------------------------------------------------
# Census and Survey Data
#-------------------------------------------------------------------------------
load("Thesis R Script/Census.C.Rdata") # Census Data
load("Thesis R Script/HIES2000.Rdata") # Survey Data
load("Thesis R Script/sample.design.Rdata") # Sampling Design Structure
load("Thesis R Script/stuc.stratum.cluster.U.M.Rdata") # Stratum Structure with Cluster Size
#-------------------------------------------------------------------------------
# Order by Sub-district (UPZ), Cluster (EA), Households (HH)
sample.design<-sample.design[order(sample.design$Stratum.ID,sample.design$Area.ID),]
stuc.stratum.cluster.U.M<-stuc.stratum.cluster.U.M[order(stuc.stratum.cluster.U.M$ID.Stratum,stuc.stratum.cluster.U.M$ID.Area.RMO),]
Census.C<-Census.C[order(Census.C$ID.UPZ,Census.C$EA,Census.C$ID.HH),]
HIES2000<-HIES2000[order(HIES2000$ID.UPZ,HIES2000$psu,HIES2000$hhold),]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Functions of ELL Methodology
#-------------------------------------------------------------------------------
ELL.PB.HM.2L.FGT.Estimator<- function(beta,var.beta,var.com.1,var.com.2,ID.D,ID.C,X.U,t,HH.Size,No.Boot){
  
  # This function is based on ParallelPackage 
  
  # Validation
  #beta<-beta.gls.2;var.beta<-var.beta.gls.2;var.com.1<-est.sigma2.2$sigma2.1;var.com.2<-est.sigma2.2$sigma2.2;
  #ID.D<-Census.data.fitted$ID.UPZ;ID.C<-Census.data.fitted$psu;
  #X.U<-x.matrix.U; t<-Census.data.fitted$lpovln;HH.Size<-Census.data.fitted$HH.size.U
  
  
  
  # This function is for estimating FGT estimates under ELL HM
  # Considering HH Size, No HH size put HH.Size=NULL
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector : It should be in original form
  # ID.D, ID.C,X.U,t,HH.Size all correspond to HH-level information
  # All variables are ordered according to ID.D
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT0=rbind(obj1$FGT0,obj2$FGT0),
      FGT1=rbind(obj1$FGT1,obj2$FGT1),
      FGT2=rbind(obj1$FGT2,obj2$FGT2)
    )
  }
  
  N<-length(ID.D)
  N.c<-as.vector(table(ID.C))
  N.d<-as.vector(table(ID.D))
  C<-length(unique(ID.C))
  
  
  r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
    
    beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
    eta.l<-rnorm(C,0,sqrt(var.com.2))
    eps.l<-rnorm(N,0,sqrt(var.com.1))
    
    if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(eta.l,N.c)+eps.l
    if (! is.null(X.U)) z.l<-as.matrix(cbind(rep(1,N),X.U))%*%t(beta.l)+rep(eta.l,N.c)+eps.l # z.l is in logarithm scale
    
    y.l<-exp(z.l) # y.l is in original scale
    
    if (is.null(HH.Size)) {
      
      index.0<-ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,mean)
      FGT1<-tapply(index.1,ID.D,mean)
      FGT2<-tapply(index.2,ID.D,mean)
      
    }
    
    if (! is.null(HH.Size)) {
      
      index.0<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT1<-tapply(index.1,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT2<-tapply(index.2,ID.D,sum)/tapply(HH.Size,ID.D,sum)
    }
    
    list(FGT0=FGT0,FGT1=FGT1,FGT2=FGT2)
  }
  
  
  F0.F11<-colMeans(r.FGT$FGT0)
  F0.F11.MSE<-apply(r.FGT$FGT0,2,sd)
  F0.F11.Q.02.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F0.F11.Q.97.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F1.F11<-colMeans(r.FGT$FGT1)
  F1.F11.MSE<-apply(r.FGT$FGT1,2,sd)
  F1.F11.Q.02.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F1.F11.Q.97.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F2.F11<-colMeans(r.FGT$FGT2)
  F2.F11.MSE<-apply(r.FGT$FGT2,2,sd)
  F2.F11.Q.02.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F2.F11.Q.97.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  list(F0.F11=F0.F11,F0.F11.MSE=F0.F11.MSE,F0.F11.Q.02.5=F0.F11.Q.02.5,F0.F11.Q.97.5=F0.F11.Q.97.5,
       F1.F11=F1.F11,F1.F11.MSE=F1.F11.MSE,F1.F11.Q.02.5=F1.F11.Q.02.5,F1.F11.Q.97.5=F1.F11.Q.97.5,
       F2.F11=F2.F11,F2.F11.MSE=F2.F11.MSE,F2.F11.Q.02.5=F2.F11.Q.02.5,F2.F11.Q.97.5=F2.F11.Q.97.5)
  
}
ELL.PB.HM.3L.FGT.Estimator<- function(beta,var.beta,var.com.1,var.com.2,var.com.3,ID.D,ID.C,X.U,t,HH.Size,No.Boot){
  
  # This function is based on ParallelPackage 
  # This function is for estimating FGT estimates
  # Basic ELL Parametric method considering Homoskedastic Variance components
  # t is a value or a vector : It should be in original form
  # y.l is in original scale
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT0=rbind(obj1$FGT0,obj2$FGT0),
      FGT1=rbind(obj1$FGT1,obj2$FGT1),
      FGT2=rbind(obj1$FGT2,obj2$FGT2)
    )
  }
  
  
  N<-length(ID.D)
  N.c<-as.vector(table(ID.C))
  C<-length(unique(ID.C))
  N.d<-as.vector(table(ID.D))
  D<-length(unique(ID.D))
  
  r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
    
    beta.l<-mvtnorm::rmvnorm(1,beta,var.beta)
    
    u.l<-rnorm(D,0,sqrt(var.com.3))
    eta.l<-rnorm(C,0,sqrt(var.com.2))
    eps.l<-rnorm(N,0,sqrt(var.com.1))
    
    if (is.null(X.U)) z.l<-cbind(rep(1,N))%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l
    if (! is.null(X.U)) z.l<-as.matrix(cbind(rep(1,N),X.U))%*%t(beta.l)+rep(u.l,N.d)+rep(eta.l,N.c)+eps.l # z.l is in logarithm scale
    
    y.l<-exp(z.l) # y.l is in original scale
    
    if (is.null(HH.Size)) {
      
      index.0<-ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,mean)
      FGT1<-tapply(index.1,ID.D,mean)
      FGT2<-tapply(index.2,ID.D,mean)
      
    }
    
    if (! is.null(HH.Size)) {
      
      index.0<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^0
      index.1<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^1
      index.2<-HH.Size*ifelse(y.l<t,1,0)*((t-y.l)/t)^2
      
      FGT0<-tapply(index.0,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT1<-tapply(index.1,ID.D,sum)/tapply(HH.Size,ID.D,sum)
      FGT2<-tapply(index.2,ID.D,sum)/tapply(HH.Size,ID.D,sum)
    }
    
    list(FGT0=FGT0,FGT1=FGT1,FGT2=FGT2)
  }
  
  F0.F11<-colMeans(r.FGT$FGT0)
  F0.F11.MSE<-apply(r.FGT$FGT0,2,sd)
  F0.F11.Q.02.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F0.F11.Q.97.5<-apply(r.FGT$FGT0,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F1.F11<-colMeans(r.FGT$FGT1)
  F1.F11.MSE<-apply(r.FGT$FGT1,2,sd)
  F1.F11.Q.02.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F1.F11.Q.97.5<-apply(r.FGT$FGT1,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  F2.F11<-colMeans(r.FGT$FGT2)
  F2.F11.MSE<-apply(r.FGT$FGT2,2,sd)
  F2.F11.Q.02.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.025,na.rm=TRUE))
  F2.F11.Q.97.5<-apply(r.FGT$FGT2,2,function(x) quantile(x,0.975,na.rm=TRUE))
  
  list(F0.F11=F0.F11,F0.F11.MSE=F0.F11.MSE,F0.F11.Q.02.5=F0.F11.Q.02.5,F0.F11.Q.97.5=F0.F11.Q.97.5,
       F1.F11=F1.F11,F1.F11.MSE=F1.F11.MSE,F1.F11.Q.02.5=F1.F11.Q.02.5,F1.F11.Q.97.5=F1.F11.Q.97.5,
       F2.F11=F2.F11,F2.F11.MSE=F2.F11.MSE,F2.F11.Q.02.5=F2.F11.Q.02.5,F2.F11.Q.97.5=F2.F11.Q.97.5)
  
}
#-------------------------------------------------------------------------------
# Generation of Population Data Set
#-------------------------------------------------------------------------------
X<-cbind(Census.C$electric,Census.C$ilattr_1,Census.C$ilattr_3,Census.C$iwater1,Census.C$ibuild_3,
         Census.C$ibuild_4,Census.C$owner,Census.C$workst2p,Census.C$workst3p,Census.C$iincom_3,
         Census.C$num_hh,Census.C$num_hh2,Census.C$hhhprmed,Census.C$literatep,Census.C$child5p,
         Census.C$rural,Census.C$mhhsize,Census.C$depratio,Census.C$paginc,Census.C$idiv_1,
         Census.C$idiv_2,Census.C$idiv_4,Census.C$idiv_5,Census.C$femalep.rural,Census.C$ibuild_3.rural,
         Census.C$owner.ibuild_4,Census.C$workst3p.rural,Census.C$num_hh.rural,Census.C$num_hh2.rural,
         Census.C$hhhprmed.rural) ; dim(X)
N<-dim(Census.C)[1]
Z<-cbind(rep(1,N),X)
lme.3.2<-lme(fixed=ln_p_expnd~1+electric+ilattr_1+ilattr_3+iwater1+ibuild_3+ibuild_4+owner+
               workst2p+workst3p+iincom_3+num_hh+num_hh2+hhhprmed+literatep+child5p+rural+mhhsize+
               depratio+paginc+idiv_1+idiv_2+idiv_4+idiv_5+femalep.rural+ibuild_3.rural+owner.ibuild_4+
               workst3p.rural+num_hh.rural+num_hh2.rural+hhhprmed.rural, random=~1|ID.UPZ/psu,data=HIES2000)
summary.lme.3.2<-summary(lme.3.2)
BETA.3<-summary.lme.3.2$tTable[,1]
Var.BETA.3<-vcov(lme.3.2)
sigma2.3<-sqrt(round(as.numeric(c(VarCorr(lme.3.2)[5,1],VarCorr(lme.3.2)[4,1],VarCorr(lme.3.2)[2,1])),5))
#-------------------------------------------------------------------------------
# Census and Survey structure 
#-------------------------------------------------------------------------------
ID.Area.U<-unique(Census.C$ID.UPZ) 
No.Area.U<-length(unique(Census.C$ID.UPZ))
N.D<-as.vector(table(Census.C$ID.UPZ))
ID.EA.U<-unique(Census.C$EA) 
No.Cluster.U<-length(unique(Census.C$EA))
N.C<-as.vector(table(Census.C$EA))

Area.RMO.s<-sample.design$Area.ID
No.Clust.U<-sample.design$No.Clust
No.Clust.s<-sample.design$No.Clust.s
n.c<-sample.design$HH.s
n.EA<-sample.design$n.EA
#-------------------------------------------------------------------------------
# Simulation Work : ELL.2L and ELL.3L using Parallel Process : FGT estimates & their MSE
# NoSim X NoBoot: 500 X 500 - FGT estimates (With MSE estimates)
# NoSim X NoBoot: 100 X 100 - MSE of FGT estimates (With MSE estimates)
#-------------------------------------------------------------------------------
data.U<-Census.C
data.U$t<-data.U$upovln
ptm <- proc.time()
NoSim<-5
NoBoot<-5

# Matrix for storing result
FGT.0.True<-matrix(0,NoSim,No.Area.U)  
FGT.1.True<-matrix(0,NoSim,No.Area.U)
FGT.2.True<-matrix(0,NoSim,No.Area.U)
FGT.0.True.s<-matrix(0,NoSim,No.Area.U)  
FGT.1.True.s<-matrix(0,NoSim,No.Area.U)
FGT.2.True.s<-matrix(0,NoSim,No.Area.U)
FGT0.ELL2<-matrix(0,NoSim,No.Area.U)  
FGT1.ELL2<-matrix(0,NoSim,No.Area.U)
FGT2.ELL2<-matrix(0,NoSim,No.Area.U)
SE.FGT0.ELL2<-matrix(0,NoSim,No.Area.U)
SE.FGT1.ELL2<-matrix(0,NoSim,No.Area.U)
SE.FGT2.ELL2<-matrix(0,NoSim,No.Area.U)
FGT0.ELL3<-matrix(0,NoSim,No.Area.U)  
FGT1.ELL3<-matrix(0,NoSim,No.Area.U)
FGT2.ELL3<-matrix(0,NoSim,No.Area.U)
SE.FGT0.ELL3<-matrix(0,NoSim,No.Area.U)
SE.FGT1.ELL3<-matrix(0,NoSim,No.Area.U)
SE.FGT2.ELL3<-matrix(0,NoSim,No.Area.U)
FGT0.CR.I.ELL2<-matrix(0,NoSim,No.Area.U)
FGT1.CR.I.ELL2<-matrix(0,NoSim,No.Area.U)
FGT2.CR.I.ELL2<-matrix(0,NoSim,No.Area.U)
FGT0.CR.I.ELL3<-matrix(0,NoSim,No.Area.U)
FGT1.CR.I.ELL3<-matrix(0,NoSim,No.Area.U)
FGT2.CR.I.ELL3<-matrix(0,NoSim,No.Area.U)
#-------------------------------------------------------------------------------#
# Simulation Start from Here 
#-------------------------------------------------------------------------------#
set.seed(2015)
for (s in 1:NoSim){
  cat(date(),"Iteration number",s,"starting","\n",fill=T)  
  e1.3<-rnorm(N,0,sigma2.3[1])
  e2<-rnorm(No.Cluster.U,0,sigma2.3[2])
  e2.3<-rep(e2,N.C)
  e3<-rnorm(No.Area.U,0,sigma2.3[3])
  e3.3<-rep(e3,N.D)
  BETA.3.s<-mvtnorm::rmvnorm(1,BETA.3,Var.BETA.3)
  data.U$lexp.3<-Z%*%t(as.matrix(BETA.3.s))+e3.3+e2.3+e1.3
  data.U$y<-exp(data.U$lexp.3)
  
  data.U$index.0<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^0
  data.U$index.1<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^1
  data.U$index.2<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^2
  
  FGT.0.True[s,]<-tapply(data.U$index.0,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.1.True[s,]<-tapply(data.U$index.1,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.2.True[s,]<-tapply(data.U$index.2,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  
  # Selection of clusters from each small area RMO 
  Cluster.s.Selected<-NULL
  for (i in 1:length(Area.RMO.s)){
    Cluster.s.Selected<-c(Cluster.s.Selected,pps.sampling(stuc.stratum.cluster.U.M$N.C[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], sample.design$No.Clust.s[sample.design$Area.ID==Area.RMO.s[i]], id = stuc.stratum.cluster.U.M$ID.EA.U[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], method = 'tille', return.PI = FALSE)$sample)
  }
  HH.Struc.s<-data.frame(Area.RMO.s=rep(Area.RMO.s,No.Clust.s),Cluster.s=Cluster.s.Selected,n.EA=rep(n.EA,No.Clust.s))
  # Selection of HHs from each selected clusters 
  HH.s.selected<-NULL
  for (h in 1:length(HH.Struc.s$Cluster.s)){
    HH.s.selected<-c(HH.s.selected,sample(data.U$ID.HH[data.U$EA==HH.Struc.s$Cluster.s[h]],HH.Struc.s$n.EA[HH.Struc.s$Cluster.s==HH.Struc.s$Cluster.s[h]]))
  }
  # selection of sample data set from the Census data set 
  data.s<-subset(data.U,data.U$ID.HH%in%HH.s.selected) 
  data.r<-subset(data.U,!data.U$ID.HH%in%HH.s.selected) 
  # Design Based estimate of poverty Indicator
  data.s$index.0<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^0
  data.s$index.1<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^1
  data.s$index.2<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^2
  # Calculation sample (Design) Poverty Estimate 
  FGT.0.True.s[s,]<-tapply(data.s$index.0,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  FGT.1.True.s[s,]<-tapply(data.s$index.1,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  FGT.2.True.s[s,]<-tapply(data.s$index.2,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  # Fitting 3-L Mixed Model and Estimation of Variance components 
  fit.lme.3<-lme(fixed=lexp.3~1+electric+ilattr_1+ilattr_3+iwater1+ibuild_3+ibuild_4+owner+
                   workst2p+workst3p+iincom_3+num_hh+num_hh2+hhhprmed+literatep+child5p+rural+mhhsize+
                   depratio+paginc+idiv_1+idiv_2+idiv_4+idiv_5+femalep.rural+ibuild_3.rural+owner.ibuild_4+
                   workst3p.rural+num_hh.rural+num_hh2.rural+hhhprmed.rural, random=~1|ID.UPZ/EA,data=data.s)
  lemda.3.lme<-nlme::VarCorr(fit.lme.3)
  summary.lme.3<-summary(fit.lme.3)
  beta.3<-summary.lme.3$tTable[,1]
  vcov.beta.3<-vcov(fit.lme.3)
  est.lemda.3<-sqrt(as.numeric(c(lemda.3.lme[5,1],lemda.3.lme[4,1],lemda.3.lme[2,1])))
  # Fitting 2-L Mixed Model and Estimation of Variance components 
  fit.lme.2<-lme(fixed=lexp.3~1+electric+ilattr_1+ilattr_3+iwater1+ibuild_3+ibuild_4+owner+
                   workst2p+workst3p+iincom_3+num_hh+num_hh2+hhhprmed+literatep+child5p+rural+mhhsize+
                   depratio+paginc+idiv_1+idiv_2+idiv_4+idiv_5+femalep.rural+ibuild_3.rural+owner.ibuild_4+
                   workst3p.rural+num_hh.rural+num_hh2.rural+hhhprmed.rural, random=~1|EA,data=data.s)
  lemda.2.lme<-nlme::VarCorr(fit.lme.2)
  summary.lme.2<-summary(fit.lme.2)
  beta.2<-summary.lme.2$tTable[,1]
  vcov.beta.2<-vcov(fit.lme.2)
  est.lemda.2<-sqrt(as.numeric(c(lemda.2.lme[2,1],lemda.2.lme[1,1])))
  Data.Boot<-data.frame(ID.UPZ=data.U$ID.UPZ,lpovln=data.U$lpovln,num_hh=data.U$num_hh)
  ELL2<-ELL.PB.HM.2L.FGT.Estimator(beta=beta.2,var.beta=vcov.beta.2,var.com.1=est.lemda.2[1]^2,var.com.2=est.lemda.2[2]^2,ID.D=data.U$ID.UPZ,ID.C=data.U$EA,X.U=as.matrix(X),t=data.U$t,HH.Size=data.U$num_hh,No.Boot=NoBoot)
  ELL3<-ELL.PB.HM.3L.FGT.Estimator(beta=beta.3,var.beta=vcov.beta.3,var.com.1=est.lemda.3[1]^2,var.com.2=est.lemda.3[2]^2,var.com.3=est.lemda.3[3]^2,ID.D=data.U$ID.UPZ,ID.C=data.U$EA,X.U=as.matrix(X),t=data.U$t,HH.Size=data.U$num_hh,No.Boot=NoBoot)
  FGT0.ELL2[s,]<-ELL2$F0.F11 
  FGT1.ELL2[s,]<-ELL2$F1.F11
  FGT2.ELL2[s,]<-ELL2$F2.F11
  SE.FGT0.ELL2[s,]<-ELL2$F0.F11.MSE
  SE.FGT1.ELL2[s,]<-ELL2$F1.F11.MSE
  SE.FGT2.ELL2[s,]<-ELL2$F2.F11.MSE
  FGT0.CR.I.ELL2[s,]<-(FGT.0.True[s,]>=ELL2$F0.F11.Q.02.5 & FGT.0.True[s,]<=ELL2$F0.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
  FGT1.CR.I.ELL2[s,]<-(FGT.1.True[s,]>=ELL2$F1.F11.Q.02.5 & FGT.1.True[s,]<=ELL2$F1.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
  FGT2.CR.I.ELL2[s,]<-(FGT.2.True[s,]>=ELL2$F2.F11.Q.02.5 & FGT.2.True[s,]<=ELL2$F2.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
  FGT0.ELL3[s,]<-ELL3$F0.F11 
  FGT1.ELL3[s,]<-ELL3$F1.F11 
  FGT2.ELL3[s,]<-ELL3$F2.F11 
  SE.FGT0.ELL3[s,]<-ELL3$F0.F11.MSE 
  SE.FGT1.ELL3[s,]<-ELL3$F1.F11.MSE 
  SE.FGT2.ELL3[s,]<-ELL3$F2.F11.MSE 
  FGT0.CR.I.ELL3[s,]<- (FGT.0.True[s,]>=ELL3$F0.F11.Q.02.5 & FGT.0.True[s,]<=ELL3$F0.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
  FGT1.CR.I.ELL3[s,]<- (FGT.1.True[s,]>=ELL3$F1.F11.Q.02.5 & FGT.1.True[s,]<=ELL3$F1.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
  FGT2.CR.I.ELL3[s,]<- (FGT.2.True[s,]>=ELL3$F2.F11.Q.02.5 & FGT.2.True[s,]<=ELL3$F2.F11.Q.97.5)*1 # Converting the logical matrix into numerical matrix
}
run.time.min <- (proc.time()-ptm)[1]/60
Simulation.ELL2.ELL3<-list(NoSim=NoSim,NoBoot=NoBoot,ID.Area.U=ID.Area.U,FGT.0.True=FGT.0.True,FGT.1.True=FGT.1.True,FGT.2.True=FGT.2.True,
                           FGT0.ELL2=FGT0.ELL2,FGT1.ELL2=FGT1.ELL2,FGT2.ELL2=FGT2.ELL2,SE.FGT0.ELL2=SE.FGT0.ELL2,SE.FGT1.ELL2=SE.FGT1.ELL2,SE.FGT2.ELL2=SE.FGT2.ELL2,
                           FGT0.ELL3=FGT0.ELL3,FGT1.ELL3=FGT1.ELL3,FGT2.ELL3=FGT2.ELL3,SE.FGT0.ELL3=SE.FGT0.ELL3,SE.FGT1.ELL3=SE.FGT1.ELL3,SE.FGT2.ELL3=SE.FGT2.ELL3,
                           FGT0.CR.I.ELL2=FGT0.CR.I.ELL2,FGT1.CR.I.ELL2=FGT1.CR.I.ELL2,FGT2.CR.I.ELL2=FGT2.CR.I.ELL2,
                           FGT0.CR.I.ELL3=FGT0.CR.I.ELL3,FGT1.CR.I.ELL3=FGT1.CR.I.ELL3,FGT2.CR.I.ELL3=FGT2.CR.I.ELL3)
#save(Simulation.ELL2.ELL3,file="Simulation.ELL2.ELL3.UPOVLN.Rdata")
#-------------------------------------------------------------------------------#
# Simulation End Here 
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------
# Functions of EBP Methodology
#-------------------------------------------------------------------------------
# Sample design matrix for ith area 
#-------------------------------------------------------------------------------
Zi  <- function(i){as.matrix(cbind(rep(1,ni.s[i]),data.s$electric[data.s$ID.UPZ==ID.Area.U[i]],data.s$ilattr_1[data.s$ID.UPZ==ID.Area.U[i]],data.s$ilattr_3[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$iwater1[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_3[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_4[data.s$ID.UPZ==ID.Area.U[i]],data.s$owner[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$workst2p[data.s$ID.UPZ==ID.Area.U[i]],data.s$workst3p[data.s$ID.UPZ==ID.Area.U[i]],data.s$iincom_3[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$num_hh2[data.s$ID.UPZ==ID.Area.U[i]],data.s$hhhprmed[data.s$ID.UPZ==ID.Area.U[i]],data.s$literatep[data.s$ID.UPZ==ID.Area.U[i]],data.s$child5p[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$mhhsize[data.s$ID.UPZ==ID.Area.U[i]],data.s$depratio[data.s$ID.UPZ==ID.Area.U[i]],data.s$paginc[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$idiv_1[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_2[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_4[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_5[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$femalep.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_3.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$owner.ibuild_4[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$workst3p.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh2.rural[data.s$ID.UPZ==ID.Area.U[i]],
                  data.s$hhhprmed.rural[data.s$ID.UPZ==ID.Area.U[i]]))}
#-------------------------------------------------------------------------------
# Non-Sample design matrix for ith area
#-------------------------------------------------------------------------------
Zr.i  <- function(i){as.matrix(cbind(rep(1,ni.r[i]),data.r$electric[data.r$ID.UPZ==ID.Area.U[i]],data.r$ilattr_1[data.r$ID.UPZ==ID.Area.U[i]],data.r$ilattr_3[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$iwater1[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_3[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_4[data.r$ID.UPZ==ID.Area.U[i]],data.r$owner[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$workst2p[data.r$ID.UPZ==ID.Area.U[i]],data.r$workst3p[data.r$ID.UPZ==ID.Area.U[i]],data.r$iincom_3[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$num_hh2[data.r$ID.UPZ==ID.Area.U[i]],data.r$hhhprmed[data.r$ID.UPZ==ID.Area.U[i]],data.r$literatep[data.r$ID.UPZ==ID.Area.U[i]],data.r$child5p[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$mhhsize[data.r$ID.UPZ==ID.Area.U[i]],data.r$depratio[data.r$ID.UPZ==ID.Area.U[i]],data.r$paginc[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$idiv_1[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_2[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_4[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_5[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$femalep.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_3.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$owner.ibuild_4[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$workst3p.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh2.rural[data.r$ID.UPZ==ID.Area.U[i]],
                  data.r$hhhprmed.rural[data.r$ID.UPZ==ID.Area.U[i]]))}
#-------------------------------------------------------------------------------
tapply.order<-function(value,area,FUN,target.order){
  # value: The values
  # area: target area
  # target.order: The order u wish to see the outcome ....
  
  raw.output<-tapply(value,area,FUN)
  data<-data.frame(key=names(raw.output), value=raw.output)
  ordered.value<-data[match(target.order, data$key),]$value
  return(ordered.value)
}
#-------------------------------------------------------------------------------
# Functions of FGT Estimates via EBP Methodology
#-------------------------------------------------------------------------------
EBP.Parallel<-function(data.s,data.r,beta.2,vcov.beta.2,estsigma2u,estsigma2e,No.Boot){
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT0=rbind(obj1$FGT0,obj2$FGT0),
      FGT1=rbind(obj1$FGT1,obj2$FGT1),
      FGT2=rbind(obj1$FGT2,obj2$FGT2)
    )
  }
  
  
  ID.Area.U<-unique(data.s$ID.UPZ)
  No.Area.U<-length(ID.Area.U)
  ni.s<-tapply.order(data.s$logy,data.s$ID.UPZ,length,unique(data.s$ID.UPZ))
  ni.r<-tapply.order(data.r$logy,data.r$ID.UPZ,length,unique(data.r$ID.UPZ))
  
  Gamma.i<-estsigma2u/(estsigma2u+estsigma2e/ni.s)
  
  # Function for selecting the vector of sampled response variable for ith sampled HHS ===================
  
  yi  <- function(i) {matrix(c(data.s$logy[data.s$ID.UPZ==ID.Area.U[i]]),ni.s[i],1)}  
  
  # Function for selecting the vector of sampled response variable for ith non-sampled HHS ===================
  yr.i  <- function(i) {matrix(c(data.r$logy[data.r$ID.UPZ==ID.Area.U[i]]),ni.r[i],1)}
  
  # Inverse of sample variance-covariance matrix of V.ds ==============================================
  V.ds.inv <- function(i) {solve(estsigma2e*(diag(rep(1,ni.s[i])))+estsigma2u*(matrix(1,ni.s[i],ni.s[i])))}
  
  r.FGT <- foreach(icount(No.Boot), .combine=f2) %dopar% {
    
    beta.b<-beta.2
    mu.dr.i<- function(i) {Zi.r.f[[i]]%*%as.matrix(beta.b)+(estsigma2u)*matrix(rep(1,ni.r[i]),ncol=1)%*%t(matrix(rep(1,ni.s[i]),ncol=1))%*%V.ds.inv(i)%*%(yi(i)-Zi.f[[i]]%*%as.matrix(beta.b))}
    v.d<- function(i) {rnorm(1,0,sqrt((estsigma2u)*(1-Gamma.i[i])))}
    e.dr<-function(i) {rnorm((ni.r[i]),0,sqrt(estsigma2u))}
    y.dr.i<- function(i) {mu.dr.i(i)+v.d(i)*matrix(1,ni.r[i],ncol=1)+matrix(e.dr(i),ncol=1)}
    
    # Combine the sampled and non-sampled part of the i=th sampled area 
    
    y.hat.s<-NULL
    y.hat.r<-NULL
    
    for (d in 1:No.Area.U){
      y.hat.s=c(y.hat.s,as.vector(exp(yi(d))))
      y.hat.r=c(y.hat.r,as.vector(exp(y.dr.i(d))))
    }
    
    data.s$y.hat<-y.hat.s
    data.r$y.hat<-y.hat.r
    
    data.u<-rbind(data.s,data.r)
    data.u<-data.u[order(data.u$ID.HH),]
    
    data.u$index.0<-data.u$num_hh*ifelse(data.u$y.hat<data.u$t,1,0)*((data.u$t-data.u$y.hat)/data.u$t)^0
    data.u$index.1<-data.u$num_hh*ifelse(data.u$y.hat<data.u$t,1,0)*((data.u$t-data.u$y.hat)/data.u$t)^1
    data.u$index.2<-data.u$num_hh*ifelse(data.u$y.hat<data.u$t,1,0)*((data.u$t-data.u$y.hat)/data.u$t)^2
    
    # Poverty Estimate ===========================================================================
    FGT0<-tapply(data.u$index.0,data.u$ID.UPZ,sum)/tapply(data.u$num_hh,data.u$ID.UPZ,sum)
    FGT1<-tapply(data.u$index.1,data.u$ID.UPZ,sum)/tapply(data.u$num_hh,data.u$ID.UPZ,sum)
    FGT2<-tapply(data.u$index.2,data.u$ID.UPZ,sum)/tapply(data.u$num_hh,data.u$ID.UPZ,sum)
    
    list(FGT0=FGT0,FGT1=FGT1,FGT2=FGT2)
  }
  
  F0.F11<-colMeans(r.FGT$FGT0)
  F1.F11<-colMeans(r.FGT$FGT1)
  F2.F11<-colMeans(r.FGT$FGT2)
  
  list(F0.F11=F0.F11,F1.F11=F1.F11,F2.F11=F2.F11)
  
}
#-------------------------------------------------------------------------------
# Functions of MSE of FGT Estimates via EBP Methodology
#-------------------------------------------------------------------------------
MSE.EBP.Parallel<-function(Data.Boot,Z,beta.2,vcov.beta.2,estsigma2u,estsigma2e,HH.sampled,NoBoot.1,NoBoot.2){
  
  ID.Area.U<-unique(Data.Boot$ID.UPZ)
  No.Area.U<-length(ID.Area.U)
  N<-length(Data.Boot$ID.UPZ)
  N.D<-tapply.order(Data.Boot$ID.HH,Data.Boot$ID.UPZ,length,unique(Data.Boot$ID.UPZ))
  
  # Storing the results for each bootstrap population 
  FGT.0.True.BS<-matrix(0,NoBoot.1,No.Area.U) # True FGT for bootstrap population
  FGT.1.True.BS<-matrix(0,NoBoot.1,No.Area.U)
  FGT.2.True.BS<-matrix(0,NoBoot.1,No.Area.U)
  
  FGT.0.EBP.BS<-matrix(0,NoBoot.1,No.Area.U) # Estimated FGT for bootstrap population
  FGT.1.EBP.BS<-matrix(0,NoBoot.1,No.Area.U)
  FGT.2.EBP.BS<-matrix(0,NoBoot.1,No.Area.U) 
  
  # Second Bootstrap start from here =========================================================
  
  for(b in 1:NoBoot.1){
    cat(date(),"Iteration number",b,"starting","\n",fill=T) 
    e1.2<-rnorm(N,0,sqrt(estsigma2e))
    e2<-rnorm(No.Area.U,0,sqrt(estsigma2u))
    e2.2<-rep(e2,N.D)
    
    Data.Boot$logy<-Z%*%as.matrix(beta.2)+e2.2+e1.2
    Data.Boot$y<-exp(Data.Boot$logy)
    
    Data.Boot$index.0<-Data.Boot$num_hh*ifelse(Data.Boot$y<Data.Boot$t,1,0)*((Data.Boot$t-Data.Boot$y)/Data.Boot$t)^0
    Data.Boot$index.1<-Data.Boot$num_hh*ifelse(Data.Boot$y<Data.Boot$t,1,0)*((Data.Boot$t-Data.Boot$y)/Data.Boot$t)^1
    Data.Boot$index.2<-Data.Boot$num_hh*ifelse(Data.Boot$y<Data.Boot$t,1,0)*((Data.Boot$t-Data.Boot$y)/Data.Boot$t)^2
    
    # Calculation of True Poverty Estimate 
    
    FGT.0.True.BS[b,]<-tapply(Data.Boot$index.0,Data.Boot$ID.UPZ,sum)/tapply(Data.Boot$num_hh,Data.Boot$ID.UPZ,sum)
    FGT.1.True.BS[b,]<-tapply(Data.Boot$index.1,Data.Boot$ID.UPZ,sum)/tapply(Data.Boot$num_hh,Data.Boot$ID.UPZ,sum)
    FGT.2.True.BS[b,]<-tapply(Data.Boot$index.2,Data.Boot$ID.UPZ,sum)/tapply(Data.Boot$num_hh,Data.Boot$ID.UPZ,sum)
    
    # Reapply the EBP method to estimate the poverty estimate from separating =========================
    
    # selection of sample data set from the Census data set 
    
    data.s.BS<-subset(Data.Boot,Data.Boot$ID.HH%in%HH.sampled) 
    data.r.BS<-subset(Data.Boot,!Data.Boot$ID.HH%in%HH.sampled) 
    
    # Fitting 2-L Mixed Model and Estimation of Variance components 
    
    fit.lme.2.BS<-lme(fixed=logy~1+electric+ilattr_1+ilattr_3+iwater1+ibuild_3+ibuild_4+owner+
                        workst2p+workst3p+iincom_3+num_hh+num_hh2+hhhprmed+literatep+child5p+rural+mhhsize+
                        depratio+paginc+idiv_1+idiv_2+idiv_4+idiv_5+femalep.rural+ibuild_3.rural+owner.ibuild_4+
                        workst3p.rural+num_hh.rural+num_hh2.rural+hhhprmed.rural, random=~1|ID.UPZ,data=data.s.BS)
    lemda.2.lme.BS<-nlme::VarCorr(fit.lme.2.BS)
    summary.lme.2.BS<-summary(fit.lme.2.BS)
    beta.2.BS<-summary.lme.2.BS$tTable[,1]
    vcov.beta.2.BS<-vcov(fit.lme.2.BS)
    
    estsigma2u.b<-as.numeric(lemda.2.lme.BS[2,1])
    estsigma2e.b<-as.numeric(lemda.2.lme.BS[1,1])
    
    # ================================================================================================
    # Estimation of Poverty Estimate Using EBP Methodology 
    # ================================================================================================
    
    FGT.EBP.BS<-EBP.Parallel(data.s.BS,data.r.BS,beta.2.BS,vcov.beta.2.BS,estsigma2u.b,estsigma2e.b,NoBoot.2)
    
    FGT.0.EBP.BS[b,]<-FGT.EBP.BS$F0.F11
    FGT.1.EBP.BS[b,]<-FGT.EBP.BS$F1.F11
    FGT.2.EBP.BS[b,]<-FGT.EBP.BS$F2.F11
  }
  
  EST.MSE.FGT.0.EBP<-colMeans((FGT.0.True.BS-FGT.0.EBP.BS)^2)
  EST.MSE.FGT.1.EBP<-colMeans((FGT.1.True.BS-FGT.1.EBP.BS)^2)
  EST.MSE.FGT.2.EBP<-colMeans((FGT.2.True.BS-FGT.2.EBP.BS)^2)
  
  
  
  list(F0.F11.MSE=EST.MSE.FGT.0.EBP,F1.F11.MSE=EST.MSE.FGT.1.EBP,F2.F11.MSE=EST.MSE.FGT.2.EBP)
  
}
#-------------------------------------------------------------------------------
# Simulation Work : Parallel Process : Estimate of FGT Indicators
# NoSim X NoBoot.2: 500 X 500 - FGT Estimates
#-------------------------------------------------------------------------------
ptm <- proc.time()
NoSim<-5
NoBoot<-5
#-------------------------------------------------------------------------------#
# Matrix for storing result
FGT.0.True<-matrix(0,NoSim,No.Area.U)  
FGT.1.True<-matrix(0,NoSim,No.Area.U)
FGT.2.True<-matrix(0,NoSim,No.Area.U)
#FGT.0.True.s<-matrix(0,NoSim,No.Area.U)  
#FGT.1.True.s<-matrix(0,NoSim,No.Area.U)
#FGT.2.True.s<-matrix(0,NoSim,No.Area.U)
FGT0.EBP<-matrix(0,NoSim,No.Area.U)  
FGT1.EBP<-matrix(0,NoSim,No.Area.U)
FGT2.EBP<-matrix(0,NoSim,No.Area.U)
#SE.FGT0.EBP<-matrix(0,NoSim,No.Area.U)
#SE.FGT1.EBP<-matrix(0,NoSim,No.Area.U)
#SE.FGT2.EBP<-matrix(0,NoSim,No.Area.U)
#FGT0.CR.I.EBP<-matrix(0,NoSim,No.Area.U)
#FGT1.CR.I.EBP<-matrix(0,NoSim,No.Area.U)
#FGT2.CR.I.EBP<-matrix(0,NoSim,No.Area.U)
#-------------------------------------------------------------------------------#
# Simulation start from Here 
set.seed(2015)
data.U<-Census.C
data.U$t<-data.U$upovln

for (s in 1:NoSim){
  
  cat(date(),"Iteration number",s,"starting","\n",fill=T)  
  
  e1.3<-rnorm(N,0,sigma2.3[1])
  
  e2<-rnorm(No.Cluster.U,0,sigma2.3[2])
  e2.3<-rep(e2,N.C)
  
  e3<-rnorm(No.Area.U,0,sigma2.3[3])
  e3.3<-rep(e3,N.D)
  
  # data.U$lexp.3<-Z%*%as.matrix(BETA.3)+e3.3+e2.3+e1.3
  BETA.3.s<-mvtnorm::rmvnorm(1,BETA.3,Var.BETA.3)
  # lexp.3<-Z%*%t(as.matrix(BETA.3.s))+e3.3+e2.3+e1.3
  
  data.U$logy<-Z%*%t(BETA.3.s)+e3.3+e2.3+e1.3
  
  data.U$y<-exp(data.U$logy)
  
  data.U$index.0<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^0
  data.U$index.1<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^1
  data.U$index.2<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^2
  
  
  # Calculation of True Poverty Estimate 
  
  FGT.0.True[s,]<-tapply(data.U$index.0,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.1.True[s,]<-tapply(data.U$index.1,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.2.True[s,]<-tapply(data.U$index.2,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  
  
  # Selection of clusters from each small area RMO 
  
  Cluster.s.Selected<-NULL
  for (i in 1:length(Area.RMO.s)){
    Cluster.s.Selected<-c(Cluster.s.Selected,pps.sampling(stuc.stratum.cluster.U.M$N.C[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], sample.design$No.Clust.s[sample.design$Area.ID==Area.RMO.s[i]], id = stuc.stratum.cluster.U.M$ID.EA.U[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], method = 'tille', return.PI = FALSE)$sample)
  }
  
  HH.Struc.s<-data.frame(Area.RMO.s=rep(Area.RMO.s,No.Clust.s),Cluster.s=Cluster.s.Selected,n.EA=rep(n.EA,No.Clust.s))
  
  # Selection of HHs from each selected clusters 
  
  HH.s.selected<-NULL
  for (h in 1:length(HH.Struc.s$Cluster.s)){
    HH.s.selected<-c(HH.s.selected,sample(Census.C$ID.HH[Census.C$EA==HH.Struc.s$Cluster.s[h]],HH.Struc.s$n.EA[HH.Struc.s$Cluster.s==HH.Struc.s$Cluster.s[h]]))
  }
  
  # selection of sample data set from the Census data set 
  
  data.s<-subset(data.U,data.U$ID.HH%in%HH.s.selected) 
  data.r<-subset(data.U,!data.U$ID.HH%in%HH.s.selected) 
  ni.s<-as.vector(table(data.s$ID.UPZ)) ## Sampled HH size for the sampled areas
  n.s<-sum(ni.s)
  ni.r<-as.vector(table(data.r$ID.UPZ)) ## Sampled HH size for the sampled areas
  n.r<-sum(ni.r)
  
  # Calculation of sample Poverty Estimate 
  
  #  data.s$index.0<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^0
  #  data.s$index.1<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^1
  #  data.s$index.2<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^2
  
  #  FGT.0.True.s[s,]<-tapply(data.s$index.0,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  #  FGT.1.True.s[s,]<-tapply(data.s$index.1,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  #  FGT.2.True.s[s,]<-tapply(data.s$index.2,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  
  # Sample design matrix for ith area ==================================================== #
  Zi  <- function(i){
    as.matrix(cbind(rep(1,ni.s[i]),data.s$electric[data.s$ID.UPZ==ID.Area.U[i]],data.s$ilattr_1[data.s$ID.UPZ==ID.Area.U[i]],data.s$ilattr_3[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$iwater1[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_3[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_4[data.s$ID.UPZ==ID.Area.U[i]],data.s$owner[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$workst2p[data.s$ID.UPZ==ID.Area.U[i]],data.s$workst3p[data.s$ID.UPZ==ID.Area.U[i]],data.s$iincom_3[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$num_hh2[data.s$ID.UPZ==ID.Area.U[i]],data.s$hhhprmed[data.s$ID.UPZ==ID.Area.U[i]],data.s$literatep[data.s$ID.UPZ==ID.Area.U[i]],data.s$child5p[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$mhhsize[data.s$ID.UPZ==ID.Area.U[i]],data.s$depratio[data.s$ID.UPZ==ID.Area.U[i]],data.s$paginc[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$idiv_1[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_2[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_4[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_5[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$femalep.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_3.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$owner.ibuild_4[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$workst3p.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh2.rural[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$hhhprmed.rural[data.s$ID.UPZ==ID.Area.U[i]]))}
  
  # Non-Sample design matrix for ith area ==================================================== #
  Zr.i  <- function(i){
    as.matrix(cbind(rep(1,ni.r[i]),data.r$electric[data.r$ID.UPZ==ID.Area.U[i]],data.r$ilattr_1[data.r$ID.UPZ==ID.Area.U[i]],data.r$ilattr_3[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$iwater1[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_3[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_4[data.r$ID.UPZ==ID.Area.U[i]],data.r$owner[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$workst2p[data.r$ID.UPZ==ID.Area.U[i]],data.r$workst3p[data.r$ID.UPZ==ID.Area.U[i]],data.r$iincom_3[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$num_hh2[data.r$ID.UPZ==ID.Area.U[i]],data.r$hhhprmed[data.r$ID.UPZ==ID.Area.U[i]],data.r$literatep[data.r$ID.UPZ==ID.Area.U[i]],data.r$child5p[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$mhhsize[data.r$ID.UPZ==ID.Area.U[i]],data.r$depratio[data.r$ID.UPZ==ID.Area.U[i]],data.r$paginc[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$idiv_1[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_2[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_4[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_5[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$femalep.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_3.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$owner.ibuild_4[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$workst3p.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh2.rural[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$hhhprmed.rural[data.r$ID.UPZ==ID.Area.U[i]]))}
  
  Zi.f <- vector("list", length = No.Area.U)
  for (i in 1:No.Area.U) Zi.f[[i]] <- Zi(i)  
  
  Zi.r.f <- vector("list", length = No.Area.U)
  for (i in 1:No.Area.U) Zi.r.f[[i]] <- Zr.i(i) 
  
  
  # Fitting 2-L Mixed Model and Estimation of Variance components 
  
  fit.lme.2<-lme(fixed=logy~1+electric+ilattr_1+ilattr_3+iwater1+ibuild_3+ibuild_4+owner+
                   workst2p+workst3p+iincom_3+num_hh+num_hh2+hhhprmed+literatep+child5p+rural+mhhsize+
                   depratio+paginc+idiv_1+idiv_2+idiv_4+idiv_5+femalep.rural+ibuild_3.rural+owner.ibuild_4+
                   workst3p.rural+num_hh.rural+num_hh2.rural+hhhprmed.rural, random=~1|ID.UPZ,data=data.s)
  lemda.2.lme<-nlme::VarCorr(fit.lme.2)
  summary.lme.2<-summary(fit.lme.2)
  beta.2<-summary.lme.2$tTable[,1]
  vcov.beta.2<-vcov(fit.lme.2)
  
  estsigma2u<-as.numeric(lemda.2.lme[2,1])
  estsigma2e<-as.numeric(lemda.2.lme[1,1])
  
  # ================================================================================================#
  # Estimation of Poverty Estimate Using EBP Methodology 
  # ================================================================================================#
  
  FGT.EBP<-EBP.Parallel(data.s,data.r,beta.2,vcov.beta.2,estsigma2u,estsigma2e,NoBoot)
  
  #  Data.Boot<-data.U
  #  SE.FGT.EBP<-MSE.EBP.Parallel(Data.Boot,Z,beta.2,vcov.beta.2,estsigma2u,estsigma2e,HH.s.selected,NoBoot.1,NoBoot.2)
  
  FGT0.EBP[s,]<-FGT.EBP$F0.F11
  FGT1.EBP[s,]<-FGT.EBP$F1.F11
  FGT2.EBP[s,]<-FGT.EBP$F2.F11
  #  SE.FGT0.EBP[s,]<-SE.FGT.EBP$F0.F11.MSE
  #  SE.FGT1.EBP[s,]<-SE.FGT.EBP$F1.F11.MSE
  #  SE.FGT2.EBP[s,]<-SE.FGT.EBP$F2.F11.MSE
  
  #  FGT0.CR.I.EBP[s,]<-(FGT.0.True[s,]>=(FGT.EBP$F0.F11-2*SE.FGT.EBP$F0.F11.MSE) & FGT.0.True[s,]<=FGT.EBP$F0.F11+2*SE.FGT.EBP$F0.F11.MSE)*1 # Converting the logical matrix into numerical matrix
  #  FGT1.CR.I.EBP[s,]<-(FGT.1.True[s,]>=(FGT.EBP$F1.F11-2*SE.FGT.EBP$F1.F11.MSE) & FGT.1.True[s,]<=(FGT.EBP$F1.F11+2*SE.FGT.EBP$F1.F11.MSE))*1 # Converting the logical matrix into numerical matrix
  #  FGT2.CR.I.EBP[s,]<-(FGT.2.True[s,]>=(FGT.EBP$F2.F11-2*SE.FGT.EBP$F2.F11.MSE) & FGT.2.True[s,]<=(FGT.EBP$F2.F11+2*SE.FGT.EBP$F2.F11.MSE))*1 # Converting the logical matrix into numerical matrix
  
}
run.time.min <- (proc.time()-ptm)[1]/60
EBP.FGT.Estimate.UPOVLN<-list(NoSim=NoSim,ID.Area.U=ID.Area.U,
                              FGT.0.True=FGT.0.True,FGT.1.True=FGT.1.True,FGT.2.True=FGT.2.True,
                              FGT0.EBP=FGT0.EBP,FGT1.EBP=FGT1.EBP,FGT2.EBP=FGT2.EBP)
#save(EBP.FGT.Estimate.UPOVLN,file="ASC-IMS 2014/EBP.FGT.Estimate.UPOVLN.Rdata")
#-------------------------------------------------------------------------------
# Simulation Work : Parallel Process : Estimated MSE of FGT Estimates
# NoSim X NoBoot.1 X NoBoot.2: 100 X 100 X 100 - MSE FGT Estimates
#-------------------------------------------------------------------------------
ptm <- proc.time()

NoSim<-2
NoBoot.1<-5
NoBoot.2<-5
# Matrix for storing result
FGT.0.True<-matrix(0,NoSim,No.Area.U)  
FGT.1.True<-matrix(0,NoSim,No.Area.U)
FGT.2.True<-matrix(0,NoSim,No.Area.U)
#FGT.0.True.s<-matrix(0,NoSim,No.Area.U)  
#FGT.1.True.s<-matrix(0,NoSim,No.Area.U)
#FGT.2.True.s<-matrix(0,NoSim,No.Area.U)
#FGT0.EBP<-matrix(0,NoSim,No.Area.U)  
#FGT1.EBP<-matrix(0,NoSim,No.Area.U)
#FGT2.EBP<-matrix(0,NoSim,No.Area.U)
SE.FGT0.EBP<-matrix(0,NoSim,No.Area.U)
SE.FGT1.EBP<-matrix(0,NoSim,No.Area.U)
SE.FGT2.EBP<-matrix(0,NoSim,No.Area.U)
#FGT0.CR.I.EBP<-matrix(0,NoSim,No.Area.U)
#FGT1.CR.I.EBP<-matrix(0,NoSim,No.Area.U)
#FGT2.CR.I.EBP<-matrix(0,NoSim,No.Area.U)
#-------------------------------------------------------------------------------#
# Simulation start from Here 

set.seed(2015)
data.U<-Census.C
data.U$t<-data.U$upovln
for (s in 1:NoSim){
  cat(date(),"Iteration number",s,"starting","\n",fill=T)  
  e1.3<-rnorm(N,0,sigma2.3[1])
  e2<-rnorm(No.Cluster.U,0,sigma2.3[2])
  e2.3<-rep(e2,N.C)
  e3<-rnorm(No.Area.U,0,sigma2.3[3])
  e3.3<-rep(e3,N.D)
  
  # data.U$lexp.3<-Z%*%as.matrix(BETA.3)+e3.3+e2.3+e1.3
  BETA.3.s<-mvtnorm::rmvnorm(1,BETA.3,Var.BETA.3)
  # lexp.3<-Z%*%t(as.matrix(BETA.3.s))+e3.3+e2.3+e1.3
  
  data.U$logy<-Z%*%t(BETA.3.s)+e3.3+e2.3+e1.3
  data.U$y<-exp(data.U$logy)
  data.U$index.0<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^0
  data.U$index.1<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^1
  data.U$index.2<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^2
  # Calculation of True Poverty Estimate 
  FGT.0.True[s,]<-tapply(data.U$index.0,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.1.True[s,]<-tapply(data.U$index.1,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.2.True[s,]<-tapply(data.U$index.2,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  # Selection of clusters from each small area RMO 
  Cluster.s.Selected<-NULL
  for (i in 1:length(Area.RMO.s)){
    Cluster.s.Selected<-c(Cluster.s.Selected,pps.sampling(stuc.stratum.cluster.U.M$N.C[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], sample.design$No.Clust.s[sample.design$Area.ID==Area.RMO.s[i]], id = stuc.stratum.cluster.U.M$ID.EA.U[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], method = 'tille', return.PI = FALSE)$sample)
  }
  HH.Struc.s<-data.frame(Area.RMO.s=rep(Area.RMO.s,No.Clust.s),Cluster.s=Cluster.s.Selected,n.EA=rep(n.EA,No.Clust.s))
  # Selection of HHs from each selected clusters 
  HH.s.selected<-NULL
  for (h in 1:length(HH.Struc.s$Cluster.s)){
    HH.s.selected<-c(HH.s.selected,sample(Census.C$ID.HH[Census.C$EA==HH.Struc.s$Cluster.s[h]],HH.Struc.s$n.EA[HH.Struc.s$Cluster.s==HH.Struc.s$Cluster.s[h]]))
  }
  # selection of sample data set from the Census data set 
  data.s<-subset(data.U,data.U$ID.HH%in%HH.s.selected) 
  data.r<-subset(data.U,!data.U$ID.HH%in%HH.s.selected) 
  ni.s<-as.vector(table(data.s$ID.UPZ)) ## Sampled HH size for the sampled areas
  n.s<-sum(ni.s)
  ni.r<-as.vector(table(data.r$ID.UPZ)) ## Sampled HH size for the sampled areas
  n.r<-sum(ni.r)
  # Calculation of sample Poverty Estimate 
  #  data.s$index.0<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^0
  #  data.s$index.1<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^1
  #  data.s$index.2<-data.s$num_hh*ifelse(data.s$y<data.s$t,1,0)*((data.s$t-data.s$y)/data.s$t)^2
  #  FGT.0.True.s[s,]<-tapply(data.s$index.0,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  #  FGT.1.True.s[s,]<-tapply(data.s$index.1,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  #  FGT.2.True.s[s,]<-tapply(data.s$index.2,data.s$ID.UPZ,sum)/tapply(data.s$num_hh,data.s$ID.UPZ,sum)
  # Sample design matrix for ith area ==================================================== #
  Zi  <- function(i){
    as.matrix(cbind(rep(1,ni.s[i]),data.s$electric[data.s$ID.UPZ==ID.Area.U[i]],data.s$ilattr_1[data.s$ID.UPZ==ID.Area.U[i]],data.s$ilattr_3[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$iwater1[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_3[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_4[data.s$ID.UPZ==ID.Area.U[i]],data.s$owner[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$workst2p[data.s$ID.UPZ==ID.Area.U[i]],data.s$workst3p[data.s$ID.UPZ==ID.Area.U[i]],data.s$iincom_3[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$num_hh2[data.s$ID.UPZ==ID.Area.U[i]],data.s$hhhprmed[data.s$ID.UPZ==ID.Area.U[i]],data.s$literatep[data.s$ID.UPZ==ID.Area.U[i]],data.s$child5p[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$mhhsize[data.s$ID.UPZ==ID.Area.U[i]],data.s$depratio[data.s$ID.UPZ==ID.Area.U[i]],data.s$paginc[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$idiv_1[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_2[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_4[data.s$ID.UPZ==ID.Area.U[i]],data.s$idiv_5[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$femalep.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$ibuild_3.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$owner.ibuild_4[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$workst3p.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh.rural[data.s$ID.UPZ==ID.Area.U[i]],data.s$num_hh2.rural[data.s$ID.UPZ==ID.Area.U[i]],
                    data.s$hhhprmed.rural[data.s$ID.UPZ==ID.Area.U[i]]))}
  # Non-Sample design matrix for ith area ==================================================== #
  Zr.i  <- function(i){
    as.matrix(cbind(rep(1,ni.r[i]),data.r$electric[data.r$ID.UPZ==ID.Area.U[i]],data.r$ilattr_1[data.r$ID.UPZ==ID.Area.U[i]],data.r$ilattr_3[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$iwater1[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_3[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_4[data.r$ID.UPZ==ID.Area.U[i]],data.r$owner[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$workst2p[data.r$ID.UPZ==ID.Area.U[i]],data.r$workst3p[data.r$ID.UPZ==ID.Area.U[i]],data.r$iincom_3[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$num_hh2[data.r$ID.UPZ==ID.Area.U[i]],data.r$hhhprmed[data.r$ID.UPZ==ID.Area.U[i]],data.r$literatep[data.r$ID.UPZ==ID.Area.U[i]],data.r$child5p[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$mhhsize[data.r$ID.UPZ==ID.Area.U[i]],data.r$depratio[data.r$ID.UPZ==ID.Area.U[i]],data.r$paginc[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$idiv_1[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_2[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_4[data.r$ID.UPZ==ID.Area.U[i]],data.r$idiv_5[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$femalep.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$ibuild_3.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$owner.ibuild_4[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$workst3p.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh.rural[data.r$ID.UPZ==ID.Area.U[i]],data.r$num_hh2.rural[data.r$ID.UPZ==ID.Area.U[i]],
                    data.r$hhhprmed.rural[data.r$ID.UPZ==ID.Area.U[i]]))}
  Zi.f <- vector("list", length = No.Area.U)
  for (i in 1:No.Area.U) Zi.f[[i]] <- Zi(i)  
  Zi.r.f <- vector("list", length = No.Area.U)
  for (i in 1:No.Area.U) Zi.r.f[[i]] <- Zr.i(i) 
  #-------------------------------------------------------------------------------#
  # Fitting 2-L Mixed Model and Estimation of Variance components 
  fit.lme.2<-lme(fixed=logy~1+electric+ilattr_1+ilattr_3+iwater1+ibuild_3+ibuild_4+owner+
                   workst2p+workst3p+iincom_3+num_hh+num_hh2+hhhprmed+literatep+child5p+rural+mhhsize+
                   depratio+paginc+idiv_1+idiv_2+idiv_4+idiv_5+femalep.rural+ibuild_3.rural+owner.ibuild_4+
                   workst3p.rural+num_hh.rural+num_hh2.rural+hhhprmed.rural, random=~1|ID.UPZ,data=data.s)
  lemda.2.lme<-nlme::VarCorr(fit.lme.2)
  summary.lme.2<-summary(fit.lme.2)
  beta.2<-summary.lme.2$tTable[,1]
  vcov.beta.2<-vcov(fit.lme.2)
  estsigma2u<-as.numeric(lemda.2.lme[2,1])
  estsigma2e<-as.numeric(lemda.2.lme[1,1])
  # ================================================================================================#
  # Estimation of Poverty Estimate Using EBP Methodology 
  # ================================================================================================#
  #FGT.EBP<-EBP.Parallel(data.s,data.r,beta.2,vcov.beta.2,estsigma2u,estsigma2e,NoBoot)
  Data.Boot<-data.U
  SE.FGT.EBP<-MSE.EBP.Parallel(Data.Boot,Z,beta.2,vcov.beta.2,estsigma2u,estsigma2e,HH.s.selected,NoBoot.1,NoBoot.2)
  #FGT0.EBP[s,]<-FGT.EBP$F0.F11
  #FGT1.EBP[s,]<-FGT.EBP$F1.F11
  #FGT2.EBP[s,]<-FGT.EBP$F2.F11
  SE.FGT0.EBP[s,]<-SE.FGT.EBP$F0.F11.MSE
  SE.FGT1.EBP[s,]<-SE.FGT.EBP$F1.F11.MSE
  SE.FGT2.EBP[s,]<-SE.FGT.EBP$F2.F11.MSE
  #  FGT0.CR.I.EBP[s,]<-(FGT.0.True[s,]>=(FGT.EBP$F0.F11-2*SE.FGT.EBP$F0.F11.MSE) & FGT.0.True[s,]<=FGT.EBP$F0.F11+2*SE.FGT.EBP$F0.F11.MSE)*1 # Converting the logical matrix into numerical matrix
  #  FGT1.CR.I.EBP[s,]<-(FGT.1.True[s,]>=(FGT.EBP$F1.F11-2*SE.FGT.EBP$F1.F11.MSE) & FGT.1.True[s,]<=(FGT.EBP$F1.F11+2*SE.FGT.EBP$F1.F11.MSE))*1 # Converting the logical matrix into numerical matrix
  #  FGT2.CR.I.EBP[s,]<-(FGT.2.True[s,]>=(FGT.EBP$F2.F11-2*SE.FGT.EBP$F2.F11.MSE) & FGT.2.True[s,]<=(FGT.EBP$F2.F11+2*SE.FGT.EBP$F2.F11.MSE))*1 # Converting the logical matrix into numerical matrix
}

run.time.min <- (proc.time()-ptm)[1]/60
EBP.FGT.EMSE<-list(ID.Area.U=ID.Area.U,
                   FGT.0.True=FGT.0.True,FGT.1.True=FGT.1.True,FGT.2.True=FGT.2.True,
                   SE.FGT0.EBP=sqrt(SE.FGT0.EBP),SE.FGT1.EBP=sqrt(SE.FGT1.EBP),SE.FGT2.EBP=sqrt(SE.FGT2.EBP))
#save(EBP.FGT.EMSE,file="Thesis R Script/EBP.FGT.EMSE.Rdata")
#-------------------------------------------------------------------------------
# Functions of MQ SAE Methodology
#-------------------------------------------------------------------------------
QRLM <- function (x, y, case.weights = rep(1, nrow(x)), var.weights = rep(1, nrow(x)), ..., w = rep(1, nrow(x)), init = "ls", psi = psi.huber, scale.est = c("MAD", "Huber", "proposal 2"), k2 = 1.345, method = c("M", "MM"), maxit = 20, acc = 1e-04, test.vec = "resid", q = 0.5)
{
  irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r*w,1,length(r)) %*% x)/sqrt(matrix(w,1,length(r)) %*% (x^2))))/sqrt(sum(w*r^2))
  }
  method <- match.arg(method)
  nmx <- deparse(substitute(x))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)
  if (is.null(colnames(x))) 
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x)) 
    stop("x is singular: singular fits are not implemented in rlm")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec))) 
    stop("invalid testvec")
  if (length(var.weights) != nrow(x)) 
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0)) 
    stop("Negative var.weights value")
  if (length(case.weights) != nrow(x)) 
    stop("Length of case.weights must equal number of observations")
  w <- (w * case.weights)/var.weights
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi)) 
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0)) 
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      if (init == "ls") 
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts") 
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init)) 
        coef <- init$coef
      else coef <- init
      resid <- y - x %*% coef
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...))) 
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else stop("method is unknown")
  done <- FALSE
  conv <- NULL
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM") 
    scale <- mad(resid/sqrt(var.weights), 0)
  theta <- 2 * pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec)) 
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD") 
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights))) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec)) 
        convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done) 
        break
    }
    if (!done) 
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
  }
  list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt, coefficients = qest)
}
#-------------------------------------------------------------------------------
# COMPUTE THE QUANTILE ORDER
# COMPUTING OF THE QUANTILE-ORDERS
"zerovalinter"<-function(y, x) {
  if(min(y) > 0) {
    xmin <- x[y == min(y)]
    if(length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  
  else {
    if(max(y) < 0) {
      xmin <- x[y == max(y)]
      if(length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if(length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if(length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if(length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if(length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
      xmin <- x1
      if(abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}
#-------------------------------------------------------------------------------
# Function for Finding the Quantile Orders by Linear Interpolation
# Assumes that "zerovalinter" function has been already loaded
"gridfitinter"<-function(y,expectile,Q)  {
  # computing of the expectile-order of each observation of y by interpolation
  nq<-length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile        
  vectordest <- apply(diff, 1, zerovalinter,Q)    
  
  #print(vectordest)
  #qord<-list(ord=c(vectordest))
  #qord
}
#-------------------------------------------------------------------------------
mq.coef<-function(myx,myy,myregioncode,maxiter=100){
  
  #This function estimate the m-quantile regression coefficients
  #myx<- x sample matrix of auxiliary variables
  #myy<- y vector
  #mynumauxvar<- number of auxiliary variables (include constant)
  #myregioncode<- area code for y and x units, data must be ordered by area code
  #maxiter<- OPTIONAL, number of maximum iteration for ob algorithm
  
  myar<-unique(myregioncode)
  myareas<-length(myar)
  mysamplesize<-sum(as.numeric(table(myregioncode)))
  mynumauxvar<-dim(myx)[2]
  
  ob<-QRLM(myx, myy, maxit=maxiter,q=sort(c(seq(0.001,0.999,0.005),0.5,0.994,0.01,0.02,0.96,0.98)),k=1.345) ## Calculation of quantiles for specified probability
  qo<-matrix(c(gridfitinter(myy,ob$fitted.values,ob$q.values)),nrow=mysamplesize,ncol=1) ## linera interpolation for calculating MQ coefficient for individual observed y's
  qmat<-matrix(c(qo,myregioncode),nrow=length(myregioncode),ncol=2) ## Arrange according to area
  mqo<-aggregate(qmat[,1],list(d2=qmat[,2]),mean)[,2] ## Calculation of Area-specific MQ coefficients
  saq<-matrix(c(mqo,myar),nrow=myareas,ncol=2) ### Arrange according to area
  saq<-rbind(saq,c(0.5,9999)) # Include one row
  ob1<-QRLM(myx, myy,maxit = maxiter,q=c(mqo[1:myareas]),k=1.345) ## Calculation of MQ for the obsereved y's
  ob2<-QRLM(myx, myy,maxit = maxiter,q=0.5,k=1.345) ## Calculation of MQ for the non-obsereved y's
  mycoef<-matrix(c(t(ob1$coefficients)),nrow=myareas,ncol=mynumauxvar) #Calculation of MQ regression coef by area # need to be ordered by area
  mycoef<-t(mycoef)
  mycoef2<-t(matrix(c(t(ob2$coefficients)),nrow=1,ncol=mynumauxvar))
  list(q.mean=mycoef,q.unit=qmat,q.50=mycoef2) ## q.mean=Regression coefficient, q.unit= Individual MQ coefficient
}
#-------------------------------------------------------------------------------
# Functions of FGT Estimates via MQ Methodology
#-------------------------------------------------------------------------------
MQ.Parallel<-function(data.s,data.r,L){
  
  f2 <- function(obj1,obj2) {
    # f2 is for single quantile
    z <- list(
      FGT.0.iMQ=rbind(obj1$FGT.0.iMQ,obj2$FGT.0.iMQ),
      FGT.1.iMQ=rbind(obj1$FGT.1.iMQ,obj2$FGT.1.iMQ),
      FGT.2.iMQ=rbind(obj1$FGT.2.iMQ,obj2$FGT.2.iMQ)
    )
  }
  
  row.names(data.s)<-c(1:dim(data.s)[1])  
  row.names(data.r)<-c(1:dim(data.r)[1])  
  
  # Variables in data.s: Dependetnt variable, explanatory variables, household size, poverty line for all HHs,
  # Identification variables from top administrative units to HHs
  # Variables in data.r: All variables except Dependetnt variable
  
  Area.ID.s<-data.s$ID.UPZ
  ID.Area.U<-unique(data.s$ID.UPZ)
  No.Area.U<-length(ID.Area.U)
  
  ni.s<-as.vector(table(data.s$ID.UPZ)) ## Sampled HH size for the sampled areas
  n.s<-sum(ni.s)
  
  ni.r<-as.vector(table(data.r$ID.UPZ)) ## Sampled HH size for the sampled areas
  n.r<-sum(ni.r)
  
  y.s<-data.s$logy
  x.s<-cbind(rep(1,n.s),data.s$electric,data.s$ilattr_1,data.s$ilattr_3,data.s$iwater1,data.s$ibuild_3,
             data.s$ibuild_4,data.s$owner,data.s$workst2p,data.s$workst3p,data.s$iincom_3,
             data.s$num_hh,data.s$num_hh2,data.s$hhhprmed,data.s$literatep,data.s$child5p,
             data.s$rural,data.s$mhhsize,data.s$depratio,data.s$paginc,data.s$idiv_1,
             data.s$idiv_2,data.s$idiv_4,data.s$idiv_5,data.s$femalep.rural,data.s$ibuild_3.rural,
             data.s$owner.ibuild_4,data.s$workst3p.rural,data.s$num_hh.rural,data.s$num_hh2.rural,
             data.s$hhhprmed.rural)
  
  # Function for selecting the vector of sampled response variable for ith sampled HHS 
  
  yi  <- function(i){matrix(c(data.s$logy[data.s$ID.UPZ==ID.Area.U[i]]),ni.s[i],1)}  
  
  # Sample design matrix for All area 
  
  Zs<-as.matrix(cbind(rep(1,n.s),data.s$electric,data.s$ilattr_1,data.s$ilattr_3,
                      data.s$iwater1,data.s$ibuild_3,data.s$ibuild_4,data.s$owner,
                      data.s$workst2p,data.s$workst3p,data.s$iincom_3,data.s$num_hh,
                      data.s$num_hh2,data.s$hhhprmed,data.s$literatep,data.s$child5p,
                      data.s$rural,data.s$mhhsize,data.s$depratio,data.s$paginc,
                      data.s$idiv_1,data.s$idiv_2,data.s$idiv_4,data.s$idiv_5,
                      data.s$femalep.rural,data.s$ibuild_3.rural,data.s$owner.ibuild_4,
                      data.s$workst3p.rural,data.s$num_hh.rural,data.s$num_hh2.rural,
                      data.s$hhhprmed.rural))
  
  # Sample design matrix for ith area 
  
  row.names(data.s)<-c(1:dim(data.s)[1])  
  
  Zi.s.f <- vector("list", length = No.Area.U)
  for (i in 1:No.Area.U) Zi.s.f[[i]] <- Zs[as.numeric(row.names(data.s[data.s$ID.UPZ==ID.Area.U[i],])),]  
  
  # Non-Sample design matrix for ALL area 
  
  Zr  <- as.matrix(cbind(rep(1,n.r),data.r$electric,data.r$ilattr_1,data.r$ilattr_3,
                         data.r$iwater1,data.r$ibuild_3,data.r$ibuild_4,data.r$owner,
                         data.r$workst2p,data.r$workst3p,data.r$iincom_3,data.r$num_hh,
                         data.r$num_hh2,data.r$hhhprmed,data.r$literatep,data.r$child5p,
                         data.r$rural,data.r$mhhsize,data.r$depratio,data.r$paginc,
                         data.r$idiv_1,data.r$idiv_2,data.r$idiv_4,data.r$idiv_5,
                         data.r$femalep.rural,data.r$ibuild_3.rural,data.r$owner.ibuild_4,
                         data.r$workst3p.rural,data.r$num_hh.rural,data.r$num_hh2.rural,
                         data.r$hhhprmed.rural))
  # Non-Sample design matrix for ith area 
  
  row.names(data.r)<-c(1:dim(data.r)[1])  
  
  
  Zi.r.f <- vector("list", length = No.Area.U)
  for (i in 1:No.Area.U) Zi.r.f[[i]] <- Zr[as.numeric(row.names(data.r[data.r$ID.UPZ==ID.Area.U[i],])),]  
  
  # Fit the  M-quantile model to get MQ regression coefficient for all areas
  
  tmp=mq.coef(x.s,y.s,Area.ID.s)
  beta.mq<-tmp$q.mean ## MQ Regression coeffficents for covered areas
  beta.mq.50<-tmp$q.50 ## MQ Regression coeffficents for uncovered areas
  
  # MQ regression coefficient for ith area
  p<-dim(x.s)[2]
  beta.i<-function(i) {matrix(beta.mq[,i],p,1)}
  
  # Predicted sample y for ith area
  ys.mq<-function(i) {Zi.s.f[[i]]%*%beta.i(i)}
  
  # Calculation of MQ sample residuals 
  res.mq<-NULL
  
  for(i in 1:No.Area.U){
    res.d.mq<-yi(i)-ys.mq(i)
    res.mq<-c(res.mq,res.d.mq)
  }
  
  
  r.FGT <- foreach(icount(L), .combine=f2) %dopar% {
    
    # Generating the y's for non-sample population units of i-th area 
    
    mu.dr.i<- function(i) {Zi.r.f[[i]]%*%beta.i(i)}
    e.dr<-function(i) sample(res.mq,ni.r[i],replace=TRUE)
    y.dr.i<- function(i) {mu.dr.i(i)+e.dr(i)}
    
    # Combine the sampled and non-sampled part of the i=th area 
    
    y.hat.s<-NULL
    y.hat.r<-NULL
    
    for (d in 1:No.Area.U){
      y.hat.s=c(y.hat.s,as.vector(exp(yi(d))))
      y.hat.r=c(y.hat.r,as.vector(exp(y.dr.i(d))))
    }
    
    data.s$y.hat<-y.hat.s
    data.r$y.hat<-y.hat.r
    
    data.u<-rbind(data.s,data.r)
    data.u<-data.u[order(data.u$ID.UPZ),]
    
    
    data.u$index.0<-data.u$num_hh*ifelse(data.u$y.hat<data.u$t,1,0)*((data.u$t-data.u$y.hat)/data.u$t)^0
    data.u$index.1<-data.u$num_hh*ifelse(data.u$y.hat<data.u$t,1,0)*((data.u$t-data.u$y.hat)/data.u$t)^1
    data.u$index.2<-data.u$num_hh*ifelse(data.u$y.hat<data.u$t,1,0)*((data.u$t-data.u$y.hat)/data.u$t)^2
    
    # Poverty Estimate ===========================================================================
    FGT.0.iMQ<-tapply(data.u$index.0,data.u$ID.UPZ,sum)/tapply(data.u$num_hh,data.u$ID.UPZ,sum)
    FGT.1.iMQ<-tapply(data.u$index.1,data.u$ID.UPZ,sum)/tapply(data.u$num_hh,data.u$ID.UPZ,sum)
    FGT.2.iMQ<-tapply(data.u$index.2,data.u$ID.UPZ,sum)/tapply(data.u$num_hh,data.u$ID.UPZ,sum)
    
    list(FGT.0.iMQ=FGT.0.iMQ,FGT.1.iMQ=FGT.1.iMQ,FGT.2.iMQ=FGT.2.iMQ)
  }
  
  
  f.MQ.0<-colMeans(r.FGT$FGT.0.iMQ)
  f.MQ.1<-colMeans(r.FGT$FGT.1.iMQ)
  f.MQ.2<-colMeans(r.FGT$FGT.2.iMQ)
  
  
  # Restriction if happens
  
  f.MQ.0[which(f.MQ.0>1)]<-1 
  f.MQ.1[which(f.MQ.1>1)]<-1
  f.MQ.2[which(f.MQ.2>1)]<-1
  
  res<-list(L=L,Area.Code=ID.Area.U,FGT0.MQ=f.MQ.0,FGT1.MQ=f.MQ.1,FGT2.MQ=f.MQ.2,res.mq=res.mq,mycoef=beta.mq,beta.mq.50=beta.mq.50)
  res
  
}
#-------------------------------------------------------------------------------
# Functions of MSE of FGT Estimates via MQ Methodology
#-------------------------------------------------------------------------------
MSE.MQ.Parallel<-function(data.U,Zi.f,beta.mq,res.s.centered,stuc.stratum.cluster.U.M,sample.design,myB,myR,L){
  
  
  ID.Area.U<-unique(data.U$ID.UPZ) 
  No.Area.U<-length(unique(data.U$ID.UPZ))
  N.D<-as.vector(table(data.U$ID.UPZ))
  
  
  myq.true.boot.0<-matrix(0,myB,No.Area.U)
  myq.true.boot.1<-matrix(0,myB,No.Area.U)
  myq.true.boot.2<-matrix(0,myB,No.Area.U)
  
  myarea.q.boot.0<-matrix(0,myB,No.Area.U)
  myarea.q.boot.1<-matrix(0,myB,No.Area.U)
  myarea.q.boot.2<-matrix(0,myB,No.Area.U)
  
  CI.boot.0<-matrix(0,No.Area.U,2)
  CI.boot.1<-matrix(0,No.Area.U,2)
  CI.boot.2<-matrix(0,No.Area.U,2)
  
  #  myarea.q.boot.0.v2<-matrix(0,myB,No.Area.U)
  #  myarea.q.boot.1.v2<-matrix(0,myB,No.Area.U)
  #  myarea.q.boot.2.v2<-matrix(0,myB,No.Area.U)
  
  myarea.q.boot.r.0<-array(0,dim=c(myR,No.Area.U,myB))
  myarea.q.boot.r.1<-array(0,dim=c(myR,No.Area.U,myB))
  myarea.q.boot.r.2<-array(0,dim=c(myR,No.Area.U,myB))
  
  # Bootstrap start from here 
  
  for (b in 1:myB){
    
    cat(date(),"Bootstrap pop",b,"\n")
    
    
    logy.boot<-NULL
    
    for (i in 1:No.Area.U){
      logy.boot.i<-Zi.f[[i]]%*%beta.mq[,i]+sample(res.s.centered,N.D[i],replace=TRUE)
      logy.boot<-c(logy.boot,logy.boot.i)
    }
    
    
    logy.boot[logy.boot<0]<-0
    data.U$logy=logy.boot
    
    data.U$index.0.b<-data.U$num_hh*ifelse(exp(data.U$logy)<data.U$t,1,0)*((data.U$t-exp(data.U$logy))/data.U$t)^0
    data.U$index.1.b<-data.U$num_hh*ifelse(exp(data.U$logy)<data.U$t,1,0)*((data.U$t-exp(data.U$logy))/data.U$t)^1
    data.U$index.2.b<-data.U$num_hh*ifelse(exp(data.U$logy)<data.U$t,1,0)*((data.U$t-exp(data.U$logy))/data.U$t)^2
    
    # True Poverty rate ===========================================================================
    
    myq.true.boot.0[b,]<-as.vector(tapply(data.U$index.0.b,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum))
    myq.true.boot.1[b,]<-as.vector(tapply(data.U$index.1.b,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum))
    myq.true.boot.2[b,]<-as.vector(tapply(data.U$index.2.b,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum))
    
    ### Starts the bootstrap drawing sample R times ===========================================
    
    #    myarea.q.boot.r.0<-matrix(0,myR,No.Area.U)
    #    myarea.q.boot.r.1<-matrix(0,myR,No.Area.U)
    #    myarea.q.boot.r.2<-matrix(0,myR,No.Area.U)
    
    
    for(rr in 1:myR){
      
      # cat(date(),"Bootstrap pop",b,"   Bootstrap sample",rr,"\n")
      
      # Selection of clusters from each small area RMO
      
      Cluster.s.Selected<-NULL
      for (i in 1:length(Area.RMO.s)){
        Cluster.s.Selected<-c(Cluster.s.Selected,pps.sampling(stuc.stratum.cluster.U.M$N.C[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], sample.design$No.Clust.s[sample.design$Area.ID==Area.RMO.s[i]], id = stuc.stratum.cluster.U.M$ID.EA.U[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], method = 'tille', return.PI = FALSE)$sample)
      }
      
      HH.Struc.s<-data.frame(Area.RMO.s=rep(Area.RMO.s,No.Clust.s),Cluster.s=Cluster.s.Selected,n.EA=rep(n.EA,No.Clust.s))
      
      # Selection of HHs from each selected clusters 
      
      HH.s.selected<-NULL
      for (h in 1:length(HH.Struc.s$Cluster.s)){
        HH.s.selected<-c(HH.s.selected,sample(data.U$ID.HH[data.U$EA==HH.Struc.s$Cluster.s[h]],HH.Struc.s$n.EA[HH.Struc.s$Cluster.s==HH.Struc.s$Cluster.s[h]]))
      }
      
      # selection of sample data set from the Census data set 
      
      myR.data.s<-subset(data.U,data.U$ID.HH%in%HH.s.selected) 
      myR.data.r<-subset(data.U,!data.U$ID.HH%in%HH.s.selected) 
      
      FGT.MQ.rr<-MQ.Parallel(myR.data.s,myR.data.r,L)
      
      myarea.q.boot.r.0[rr,,b]<-FGT.MQ.rr$FGT0.MQ
      myarea.q.boot.r.1[rr,,b]<-FGT.MQ.rr$FGT1.MQ
      myarea.q.boot.r.2[rr,,b]<-FGT.MQ.rr$FGT2.MQ
      
      #      myarea.q.boot.r.0[rr,]<-FGT.MQ.rr$FGT0.MQ
      #      myarea.q.boot.r.1[rr,]<-FGT.MQ.rr$FGT1.MQ
      #      myarea.q.boot.r.2[rr,]<-FGT.MQ.rr$FGT2.MQ
      
      
    }
    
    myarea.q.boot.0[b,]<-colMeans(myarea.q.boot.r.0[,,b])
    myarea.q.boot.1[b,]<-colMeans(myarea.q.boot.r.1[,,b])
    myarea.q.boot.2[b,]<-colMeans(myarea.q.boot.r.2[,,b])
    
    #    myarea.q.boot.0.v2[b,]<-colMeans((myarea.q.boot.r.0[,,b]-matrix(myarea.q.boot.0[b,],myR,No.Area.U,byrow=TRUE))^2)
    #    myarea.q.boot.1.v2[b,]<-colMeans((myarea.q.boot.r.1[,,b]-matrix(myarea.q.boot.1[b,],myR,No.Area.U,byrow=TRUE))^2)
    #    myarea.q.boot.2.v2[b,]<-colMeans((myarea.q.boot.r.2[,,b]-matrix(myarea.q.boot.2[b,],myR,No.Area.U,byrow=TRUE))^2)
    
    #    myarea.q.boot.0[b,]<-colMeans(myarea.q.boot.r.0)
    #    myarea.q.boot.1[b,]<-colMeans(myarea.q.boot.r.1)
    #    myarea.q.boot.2[b,]<-colMeans(myarea.q.boot.r.2)
    
    #    myarea.q.boot.0.v2[b,]<-colMeans((myarea.q.boot.r.0-matrix(myarea.q.boot.0[b,],myR,No.Area.U,byrow=TRUE))^2)
    #    myarea.q.boot.1.v2[b,]<-colMeans((myarea.q.boot.r.1-matrix(myarea.q.boot.1[b,],myR,No.Area.U,byrow=TRUE))^2)
    #    myarea.q.boot.2.v2[b,]<-colMeans((myarea.q.boot.r.2-matrix(myarea.q.boot.2[b,],myR,No.Area.U,byrow=TRUE))^2)
    
  }
  
  # Bootstrap stop here
  
  # Averaging the B bootstrap results: True & Estimates =====================================
  
  myq.true.boot.0.m<-colMeans(myq.true.boot.0)
  myq.true.boot.1.m<-colMeans(myq.true.boot.1)
  myq.true.boot.2.m<-colMeans(myq.true.boot.2)
  
  myarea.q.boot.0.m<-colMeans(myarea.q.boot.0)
  myarea.q.boot.1.m<-colMeans(myarea.q.boot.1)
  myarea.q.boot.2.m<-colMeans(myarea.q.boot.2)
  
  bias.boot.0<-myarea.q.boot.0.m-myq.true.boot.0.m
  bias.boot.1<-myarea.q.boot.1.m-myq.true.boot.1.m
  bias.boot.2<-myarea.q.boot.2.m-myq.true.boot.2.m
  
  # Bias and variance calculation =======================================
  
  aux.0<-matrix(0,myB,No.Area.U)
  aux.1<-matrix(0,myB,No.Area.U)
  aux.2<-matrix(0,myB,No.Area.U)
  
  for (b in 1:myB){
    aux.0[b,]<-colMeans((myarea.q.boot.r.0[,,b]-matrix(myarea.q.boot.0[b,],myR,No.Area.U,byrow=TRUE))^2)
    aux.1[b,]<-colMeans((myarea.q.boot.r.1[,,b]-matrix(myarea.q.boot.1[b,],myR,No.Area.U,byrow=TRUE))^2)
    aux.2[b,]<-colMeans((myarea.q.boot.r.2[,,b]-matrix(myarea.q.boot.2[b,],myR,No.Area.U,byrow=TRUE))^2)
  }
  
  
  var.boot.0<-colMeans(aux.0)
  var.boot.1<-colMeans(aux.1)
  var.boot.2<-colMeans(aux.2)
  
  #  var.boot.0<-colMeans(myarea.q.boot.0.v2)
  #  var.boot.1<-colMeans(myarea.q.boot.1.v2)
  #  var.boot.2<-colMeans(myarea.q.boot.2.v2)
  
  
  # MSe and CI calculation =======================================
  
  RMSE.boot.0<-sqrt(bias.boot.0^2+var.boot.0)
  RMSE.boot.1<-sqrt(bias.boot.1^2+var.boot.1)
  RMSE.boot.2<-sqrt(bias.boot.2^2+var.boot.2)
  
  
  for (d in 1:No.Area.U){
    CI.boot.0[d,]<-quantile(c(myarea.q.boot.r.0[,d,]),prob=c(0.025,0.975))
    CI.boot.1[d,]<-quantile(c(myarea.q.boot.r.1[,d,]),prob=c(0.025,0.975))
    CI.boot.2[d,]<-quantile(c(myarea.q.boot.r.2[,d,]),prob=c(0.025,0.975))
  }
  
  
  list(myB=myB,myR=myR,L=L,ID.Area.U=ID.Area.U,SE.FGT0.MQ=RMSE.boot.0,SE.FGT1.MQ=RMSE.boot.1,SE.FGT2.MQ=RMSE.boot.2,
       CI.FGT0.MQ=CI.boot.0,CI.FGT1.MQ=CI.boot.1,CI.FGT2.MQ=CI.boot.2)
  
}
#-------------------------------------------------------------------------------
# Simulation Work : Parallel Work - FGT estimates with MSE Calculation
# NoSim X L: 500 X 500 - FGT Estimates
# NoSim X myB X myR X L: 100 X 100 X 5 X 100 - MSE of FGT Estimates
#-------------------------------------------------------------------------------
load("Thesis R Script/Census.C.Rdata") # Census Data
data.U<-Census.C
data.U<-data.U[order(data.U$ID.UPZ,data.U$EA,data.U$ID.HH),]
P<-length(BETA.3)
n<-sum(n.c)

set.seed(2015)
ptm <- proc.time()
NoSim<-2
L<-2
myB=2
myR=2
# Matrix for storing result
FGT.0.True<-matrix(0,NoSim,No.Area.U)
FGT.1.True<-matrix(0,NoSim,No.Area.U)
FGT.2.True<-matrix(0,NoSim,No.Area.U)

FGT0.MQ<-matrix(0,NoSim,No.Area.U)
FGT1.MQ<-matrix(0,NoSim,No.Area.U)
FGT2.MQ<-matrix(0,NoSim,No.Area.U)

SE.FGT0.MQ<-matrix(0,NoSim,No.Area.U)
SE.FGT1.MQ<-matrix(0,NoSim,No.Area.U)
SE.FGT2.MQ<-matrix(0,NoSim,No.Area.U)

CI.FGT0.MQ<-array(0,dim=c(No.Area.U,2,NoSim))
CI.FGT1.MQ<-array(0,dim=c(No.Area.U,2,NoSim))
CI.FGT2.MQ<-array(0,dim=c(No.Area.U,2,NoSim))

data.U$t<-data.U$upovln

Z.U<-as.matrix(cbind(rep(1,N),data.U$electric,data.U$ilattr_1,data.U$ilattr_3,
                     data.U$iwater1,data.U$ibuild_3,data.U$ibuild_4,data.U$owner,
                     data.U$workst2p,data.U$workst3p,data.U$iincom_3,data.U$num_hh,
                     data.U$num_hh2,data.U$hhhprmed,data.U$literatep,data.U$child5p,
                     data.U$rural,data.U$mhhsize,data.U$depratio,data.U$paginc,
                     data.U$idiv_1,data.U$idiv_2,data.U$idiv_4,data.U$idiv_5,
                     data.U$femalep.rural,data.U$ibuild_3.rural,data.U$owner.ibuild_4,
                     data.U$workst3p.rural,data.U$num_hh.rural,data.U$num_hh2.rural,
                     data.U$hhhprmed.rural))
# Sample design matrix for ith area
row.names(data.U)<-c(1:dim(data.U)[1])

Zi.f <- vector("list", length = No.Area.U)
for (i in 1:No.Area.U) Zi.f[[i]] <-Z.U[as.numeric(row.names(data.U[data.U$ID.UPZ==ID.Area.U[i],])),]
# Simulation start from here ----------------------#
for (s in 1:NoSim){
  cat(date(),"Iteration number",s,"starting","\n",fill=T)
  e1.3<-rnorm(N,0,sigma2.3[1])
  e2<-rnorm(No.Cluster.U,0,sigma2.3[2])
  e2.3<-rep(e2,N.C)
  e3<-rnorm(No.Area.U,0,sigma2.3[3])
  e3.3<-rep(e3,N.D)
  BETA.3.s<-mvtnorm::rmvnorm(1,BETA.3,Var.BETA.3)
  data.U$logy<-Z%*%t(BETA.3.s)+e3.3+e2.3+e1.3
  data.U$y<-exp(data.U$logy)
  data.U$index.0<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^0
  data.U$index.1<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^1
  data.U$index.2<-data.U$num_hh*ifelse(data.U$y<data.U$t,1,0)*((data.U$t-data.U$y)/data.U$t)^2
  # Calculation of True Poverty Estimate
  FGT.0.True[s,]<-tapply(data.U$index.0,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.1.True[s,]<-tapply(data.U$index.1,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  FGT.2.True[s,]<-tapply(data.U$index.2,data.U$ID.UPZ,sum)/tapply(data.U$num_hh,data.U$ID.UPZ,sum)
  # Selection of clusters from each small area RMO
  Cluster.s.Selected<-NULL
  for (i in 1:length(Area.RMO.s)){
    Cluster.s.Selected<-c(Cluster.s.Selected,pps.sampling(stuc.stratum.cluster.U.M$N.C[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], sample.design$No.Clust.s[sample.design$Area.ID==Area.RMO.s[i]], id = stuc.stratum.cluster.U.M$ID.EA.U[stuc.stratum.cluster.U.M$ID.Area.RMO==Area.RMO.s[i]], method = 'tille', return.PI = FALSE)$sample)
  }
  HH.Struc.s<-data.frame(Area.RMO.s=rep(Area.RMO.s,No.Clust.s),Cluster.s=Cluster.s.Selected,n.EA=rep(n.EA,No.Clust.s))
  # Selection of HHs from each selected clusters
  HH.s.selected<-NULL
  for (h in 1:length(HH.Struc.s$Cluster.s)){
    HH.s.selected<-c(HH.s.selected,sample(data.U$ID.HH[data.U$EA==HH.Struc.s$Cluster.s[h]],HH.Struc.s$n.EA[HH.Struc.s$Cluster.s==HH.Struc.s$Cluster.s[h]]))
  }
  # selection of sample data set from the Census data set
  data.s<-subset(data.U,data.U$ID.HH%in%HH.s.selected)
  data.r<-subset(data.U,!data.U$ID.HH%in%HH.s.selected)
  FGT.MQ<-MQ.Parallel(data.s,data.r,L)
  beta.mq<-FGT.MQ$mycoef
  res.s.centered<-sort(FGT.MQ$res.mq-mean(FGT.MQ$res.mq)) # Centered Unconditional
  data.U.MSE<-data.U
  SE.FGT.MQ<-MSE.MQ.Parallel(data.U=data.U.MSE,Zi.f,beta.mq,res.s.centered,stuc.stratum.cluster.U.M,sample.design,myB,myR,L)
  FGT0.MQ[s,]<-FGT.MQ$FGT0.MQ
  FGT1.MQ[s,]<-FGT.MQ$FGT1.MQ
  FGT2.MQ[s,]<-FGT.MQ$FGT2.MQ
  SE.FGT0.MQ[s,]<-SE.FGT.MQ$SE.FGT0.MQ
  SE.FGT1.MQ[s,]<-SE.FGT.MQ$SE.FGT1.MQ
  SE.FGT2.MQ[s,]<-SE.FGT.MQ$SE.FGT2.MQ
  CI.FGT0.MQ[,,s]<-SE.FGT.MQ$CI.FGT0.MQ
  CI.FGT1.MQ[,,s]<-SE.FGT.MQ$CI.FGT1.MQ
  CI.FGT2.MQ[,,s]<-SE.FGT.MQ$CI.FGT2.MQ
}
run.time.min <- (proc.time()-ptm)[1]/60
MSE.MQ.FGT.Estimate<-list(NoSim=NoSim,L=L,B=myB,R=myR,ID.Area.U=ID.Area.U,
                              FGT.0.True=FGT.0.True,FGT.1.True=FGT.1.True,FGT.2.True=FGT.2.True,
                              FGT0.MQ=FGT0.MQ,FGT1.MQ=FGT1.MQ,FGT2.MQ=FGT2.MQ,
                              SE.FGT0.MQ=SE.FGT0.MQ,SE.FGT1.MQ=SE.FGT1.MQ,SE.FGT2.MQ=SE.FGT2.MQ,
                              CI.FGT0.MQ=CI.FGT0.MQ,CI.FGT1.MQ=CI.FGT1.MQ,CI.FGT2.MQ=CI.FGT2.MQ)
#save(MSE.MQ.FGT.Estimate,file="Thesis R Script/MSE.MQ.FGT.Estimate.UPOVLN.Rdata")


