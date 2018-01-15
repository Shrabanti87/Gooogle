# Gooogle (Group Regularization for Zero Inflated Count Regression Models)

## Introduction
Zero-inflated count data are omnipresent in many fields including health care research and actuarial science. Zero-inflated Poisson (ZIP) and Zero-inflated Negative Binomial (ZINB) regression are commonly used to model these outcomes. However, when the features to be associated possess an inherent grouping structure, traditional variable selection approaches are known to produce nonsensical results. In order to be able to perform group variable selection in ZIP/ZINB models, we extend various commonly used group regularizations such as group LASSO and group bridge to ZIP/ZINB models. These mixture models typically include a logistic component to model the presence of excess zeros and a Poisson/negative Binomial component to model the count data. 

The details of the statistical model are as follows:

<img src="misc/github.pdf" width="600" align="center">

With the above formulation, we are able to achieve bi-level variable selection both zero and count models. The tuning parameter of the final model can be chosen according to the minimum AIC/BIC criteria.  

You can install our Gooogle package from Github
```r
install.packages("devtools")
devtools::install_github("himelmallick/Gooogle")
library(Gooogle)
```

## Basic Usage

```r
gooogle(data=data,yvar=yvar,xvars=xvars,zvars=xvars,group=rep(1,14),dist="poisson",penalty="gBridge")
```

- **data**: the dataset (in data frame) to be used for the analysis. 
- **yvar**:  the outcome variable name. 
- **xvars**: the vector of variable names to be included in count model.
- **zvars**: the vector of variable names to be included in zero model.
- **group**: the vector of integers indicating the grouping structure among predictors. 
- **dist**: the distribution of count model (Poisson or Negative Binomial).   
- **penalty**: the penalty to be applied for regularization. For group selection, it is one of grLasso, grMCP or grSCAD while for bi-level selection it is gBridge.  

For greatest efficiency and least ambiguity, it is best if group is a vector of consecutive integers. If any coefficients are to be included in the model without being penalized, their grouping index should be zero. 

The gooogle function will return a list containing the following objects:
- **coefficients**: a list containing the estimates for the count and zero inflation model.  
- **aic**: AIC for the fitted model.  
- **bic**: BIC for the fitted model.  
- **loglik**: Log-likelihood for the fitted model.

## Examples

#### Simulated data 1 
We construct a simulation function similar to Huang et al. (2009) [ https://doi.org/10.1093/biomet/asp020] which can be used to simulate a dataset according to the zero-inflated negative binomial model. We use 40 covariates in 5 groups with 8 covariates in each group. For this example, we assume the covariates in the count part (X) and in the zero inflation part (Z) to be the same (i.e set Z=X).

```r
library(mpath)

data.sim.func<-function(n,p,ngrp,beta,gamma,rho,family,seedval)
{
  set.seed(seedval)
  
  R<-matrix(rnorm(n*p),n,p)
  V<-matrix(0,ngrp,ngrp)
  for(i in 1:ngrp)
  {
    for(j in 1:ngrp)
    {
      V[i,j]=rho^(abs(i-j))
    }
  }
  Z<-mvrnorm(n,mu=rep(0,ngrp),Sigma=V)
  
  X<-matrix(0,n,p)
  size=rep(p/ngrp,ngrp)
  for(g in 1:ngrp)
  {
    for (j in 1:size[g])
    {
      X[,(g-1)*size[g]+j]<-(Z[,g]+R[,(g-1)*size[g]+j])/sqrt(2)
    }
  }
  X<-scale(X)
  colnames(X)<-paste("X",c(1:ncol(X)),sep="")
  xvars=colnames(X)
  zvars=xvars
  
  ## capture zero inflation ##
  y<-rzi(n,x=X,z=X,a=beta,b=gamma,family=family)
  
  data<-cbind.data.frame(y,X)
  return(list(data=data,yvar="y",xvars=xvars,zvars=zvars))
}

data.sim<-data.sim.func(n=500,p=40,ngrp=5,rho=0.4,family="poisson",beta=c(5,-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24)),gamma=c(0.5,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24)),seedval=1)
```

The simulation function returns a dataset with X (the predictor matrix corresponding to count model), Z (the predictor matrix corresponding to excess zero model) and the outcome with zero abundance (Y) simulated according to the above parameter settings. 

We can do both group level and bi-level selection in zero inflated count data using gooogle function. If we want to do a group level selection we can use any of the penalties among grLasso, grMCP or grSCAD. Below is an example of using gooogle in the simulated dataset using gBridge.

```r
data=data.sim$data
yvar<-data.sim$yvar
xvars<-data.sim$xvars
zvars<-data.sim$zvars
group=rep(1:5,each=8)

fit.gooogle <- gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty="gBridge")
fit.gooogle
```

<!---
Similarly we can do a bi-level selection on the simulated data using gBridge penalty in the gooogle function.

```r
fit.gooogle <- gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty="gBridge")
fit.gooogle
```
-->

#### Simulated data 2  
We consider another simulation study following Park and Yoon (2011) containing a mixture of continuous and categorical predictors which are used to simulate a dataset according to the zero-inflated Negative-Binomial model. We use 30 covariates of 10 continuous divided into 6 groups and 20 categorical divided into 5 equal groups. We assume the covariates in the count part (X) and in the zero part (Z) to be the same (i.e set Z=X).

```r
library(mpath)
library(dummies)

data.sim.func<-function(n,size,beta,gamma,rho,family,seedval) 
{
  set.seed(seedval)
  
  p1=6;p2=5;ngrp1=6;ngrp2=5;
  R<-matrix(rnorm(n*p1),n,p1)
  w<- rnorm(n,0,1)
  x1<-matrix(0,n,p1)
  
  for(g in 1:ngrp1)
  {
    x1[,g]<- (R[,g]+w)/(sqrt(2))
  }
  xc<-cbind(x1[,1],x1[,2],x1[,3],(x1[,3])^2,(x1[,3])^3,x1[,4],x1[,5],x1[,6],(x1[,6])^2,(x1[,6])^3)
  
  V<-matrix(0,ngrp2,ngrp2)
  
  for(i in 1:ngrp2)
  {
    for(j in 1:ngrp2)
    {
      V[i,j]=rho^(abs(i-j))
    }
  }
  
  x2<-matrix(mvrnorm(n*p2,mu=rep(0,ngrp2),Sigma=V),n,p2)
  
  for(i in 1:n)
  {
    for (j in 1:p2)
    {
      if (x2[i,j] < qnorm(0.2,0,1)){
        x2[i,j]=0
      } else if (qnorm(0.2,0,1) < x2[i,j] && x2[i,j] < qnorm(0.4,0,1)){
        x2[i,j]=1
      } else if (qnorm(0.4,0,1) < x2[i,j] && x2[i,j] < qnorm(0.6,0,1)){
        x2[i,j]=2
      } else if (qnorm(0.6,0,1) < x2[i,j] && x2[i,j] < qnorm(0.8,0,1)){
        x2[i,j]=3
      } else {
        x2[i,j]=4
      }
    }
  }
  
  xd<-NULL
  
  for(j in 1:p2)
  {
    xd<-cbind(xd,dummy(x2[,j]))
  }
  xd<-xd[,-seq(p2,p2*ngrp2,ngrp2)]
  
  X<-cbind(scale(xc),xd)
  colnames(X)<-paste("X",c(1:ncol(X)),sep="")
  xvars=colnames(X)
  zvars=xvars
  
  y<-rzi(n,x=X,z=X,a=beta,b=gamma,family=family)
  
  data<-cbind.data.frame(y,X)
  return(list(data=data,yvar="y",xvars=xvars,zvars=zvars,zeroinfl=zeroinfl))
}

betag1<-c(0)
betag2<-c(0)
betag3<-c(-0.1,0.2,0.1)
betag4<-c(0)
betag5<-c(0)
betag6<-c(2/3,-1,1/3)
betag7<-c(-2,-1,1,2)
betag8<-c(0,0,0,0)
betag9<-c(0,0,0,0)
betag10<-rep(0,4)
betag11<-c(0,0,0,0)


beta<-c(5,betag1,betag2,betag3,betag4,betag5,betag6,betag7,betag8,betag9,betag10,betag11)

gamma<-c(-0.15,beta[-1])
data.sim<-data.sim.func(n=400,beta=beta,gamma=gamma,rho=0.4,family="negbin",seedval=1)
    
```
The simulation function returns a dataset with X (the predictor matrix corresponding to count model), Z (the predictor matrix corresponding to zero inflation model) and the outcome with zero abundance (Y) simulated according to the above parameter settings. 

We can do both group level and bi-level selection in zero inflated count data using gooogle function. If we want to do a group level selection we can use any of the penalties among grLasso, grMCP or grSCAD. Below is an example of using gooogle in the simulated dataset using gBridge.

```r
data<-data.sim$data
yvar<-data.sim$yvar
xvars<-data.sim$xvars
zvars<-data.sim$zvars
size=c(1,1,3,1,1,3,4,4,4,4,4)
ngrp<-length(size)

# group
group<-NULL
for (k in 1:ngrp)
{
  group<-c(group,rep(k,size[k]))
}

fit.gooogle <- gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty="gBridge")
fit.gooogle
```

<!---
Similarly we can do a bi-level selection on the simulated data using gBridge penalty in the gooogle function.

```r
fit.gooogle <- gooogle(data=data.sim,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty="gBridge")
fit.gooogle
```
-->

#### Real data  
Let's try one example on the real data for which we are using docvisit dataset from library zic. Similar to previous studies (Jochmann, 2013), we express each continuous predictor as a group of three cubic spline variables, resulting in 24 candidate predictors with 5 triplets and and 9 singleton groups.

####################
```r

library(zic)
library(splines)

data("docvisits")

n<-nrow(docvisits)
age<-bs(docvisits$age,3)[1:n,]
hlth<-bs(docvisits$health,3)[1:n,]
hdeg<-bs(docvisits$hdegree,3)[1:n,]
schl<-bs(docvisits$schooling,3)[1:n,]
hhin<-bs(docvisits$hhincome,3)[1:n,]

attach(docvisits)
doc.spline<-cbind.data.frame(docvisits$docvisits,age,hlth,hdeg,schl,hhin,handicap,married,children,self,civil,bluec,employed,public,addon)

names(doc.spline)[1:16]<-c("docvisits",paste("age",1:3,sep=""),paste("health",1:3,sep=""),paste("hdegree",1:3,sep=""),paste("schooling",1:3,sep=""),paste("hhincome",1:3,sep=""))
data<-doc.spline
```
#####################################################################

Considering the grouping structure among the variables for age, health, hdegree, schooling and hhincome, we can use our algorithm to perform group level or bi level variable selection. Below is an example of implementation of gooogle function using gbridge penalty.

```r
group=c(rep(1:5,each=3),(6:14))

yvar<-names(data)[1]
xvars<-names(data)[-1]
zvars<-xvars

fit.gooogle <- gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty="gBridge")
fit.gooogle
```

The code below uses 100 iterations of 5 fold CV on the docvisit data to calculate MAE and MASE and demonstrates the performance of Gooogle methods as compared to EM LASSO.

```r
library(cvTools)
library(forecast)
library(mpath)

ITER<- 100
K<-5
set.seed(123)

penalties<-c("grLasso", "grMCP", "grSCAD", "gBridge")

folds <- cvFolds(nrow(data), K = K, R = ITER)
cv<-cbind.data.frame(folds$which,folds$subsets)
names(cv)<-c("fold",paste("iter",1:ITER,sep=""))

################## compute MAE and MASE for different penalties and likelihood #############
measures.list<-list()
coeff.list<-list()

for(penalty in penalties)
{
  measures<-NULL
  coeff.count.vec<-NULL
  coeff.zero.vec<-NULL

  for (k in 1:ITER)
  {
    print(c(penalty,k))

    y.test.comp.vec<-NULL
    y.pred.comp.vec<-NULL
    y.train.comp.vec<-NULL

    time.taken<-0
    for (i in 1:K)
    {
      train.idx<-folds$subsets[folds$which!=i,k]
      train<-data[train.idx,]
      test<-data[-train.idx,]

      ptm<-proc.time()
      fit<-gooogle(data=train,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty=penalty)

      time.taken<-time.taken+round((proc.time()-ptm)[3],3)

      z.test<-as.matrix(cbind(1,test[,zvars]))
      phi.hat<-1/(1+exp(-z.test%*%fit$coefficients$zero))
      x.test<-as.matrix(cbind(1,test[,xvars]))
      lam.hat<-exp(x.test%*%fit$coefficients$count)
      y.pred<-(1-phi.hat)*lam.hat
      y.test<-test[,yvar]
      y.train<-train[,yvar]

      coeff.count<-fit$coefficients$count
      coeff.zero<-fit$coefficients$zero

      y.test.comp.vec<-c(y.test.comp.vec,y.test)
      y.pred.comp.vec<-c(y.pred.comp.vec,y.pred)
      y.train.comp.vec<-c(y.train.comp.vec,y.train)

      coeff.count.vec<-rbind(coeff.count.vec,c(penalty=penalty,iteration=k,fold=i,coeff.count))
      coeff.zero.vec<-rbind(coeff.zero.vec,c(penalty=penalty,iteration=k,fold=i,coeff.zero))
    }

    forecast <- structure(list(mean=y.pred.comp.vec, fitted=y.test.comp.vec, x=y.train.comp.vec), class='forecast')

    measures<-rbind(measures,c(iter=k,accuracy(forecast,y.test)[2,c(3,6)],time.taken=time.taken))
  }
  measures.list[[paste(penalty,sep="")]]<-measures
  coeff.list[[paste(penalty,sep="")]]<-list(count=coeff.count.vec,zero=coeff.zero.vec)
}

# --------------- Compute same measures for EM ---------------- #

fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+")),sep=""))

coeff.count.vec_EM<-NULL
coeff.zero.vec_EM<-NULL
measures_EM<-NULL

for (k in 1:ITER)
{
  print(c("EM_Lasso",k))

  y.test.comp.vec<-NULL
  y.pred.comp.vec<-NULL
  y.train.comp.vec<-NULL

  time.taken<-0
  for (i in 1:K)
  {
    train.idx<-folds$subsets[folds$which!=i,k]
    train<-data[train.idx,]
    test<-data[-train.idx,]

    ptm<-proc.time()
    fit.em<-zipath(formula=fit.formula,data=train,family="negbin")
    time.taken<-time.taken+round((proc.time()-ptm)[3],3)

    bic.idx<-which.min(fit.em$bic)
    gammahat<-fit.em$coefficients$zero[,bic.idx]
    betahat<-fit.em$coefficients$count[,bic.idx]
    z.test<-as.matrix(cbind(1,test[,zvars]))
    phi.hat<-1/(1+exp(-z.test%*%gammahat))
    x.test<-as.matrix(cbind(1,test[,xvars]))
    lam.hat<-exp(x.test%*%betahat)
    y.pred<-(1-phi.hat)*lam.hat
    y.test<-test[,yvar]
    y.train<-train[,yvar]
    coeff.count_EM<-betahat
    coeff.zero_EM<-gammahat

    y.test.comp.vec<-c(y.test.comp.vec,y.test)
    y.pred.comp.vec<-c(y.pred.comp.vec,y.pred)
    y.train.comp.vec<-c(y.train.comp.vec,y.train)
    coeff.count.vec_EM<-rbind(coeff.count.vec_EM,c(penalty=penalty,iteration=k,fold=i,coeff.count_EM))
    coeff.zero.vec_EM<-rbind(coeff.zero.vec_EM,c(penalty=penalty,iteration=k,fold=i,coeff.zero_EM))
  }

  forecast <- structure(list(mean=y.pred.comp.vec, fitted=y.test.comp.vec, x=y.train.comp.vec), class='forecast')

  measures_EM<-rbind(measures_EM,c(iter=k,accuracy(forecast,y.test)[2,c(3,6)],time=time.taken))
}
coeff.list[["EM"]]<-list(count=coeff.count.vec_EM,zero=coeff.zero.vec_EM)

measures.list[["EM"]]<-measures_EM

########### doc.cv_dist is the final output containing the measures for all methods including EM and coefficients are saved for each fold within each iteration #########################

doc.cv<-list(measures=measures.list,coeff=coeff.list,cv=cv)
```

## Citation

Huang, J., S. Ma, H. Xie, and C. Zhang (2009). A group bridge approach for variable selection. Biometrika 96(2), 339–355.

Park, C. and Y. Yoon (2011). Bridge regression: adaptivity and group selection. Journal of
Statistical Planning and Inference 141 (11), 3506{3519.

Jochmann, M. (2013). What belongs where? variable selection for zero-in
ated count models with an application to the demand for health care. Computational Statistics 28, 1947{1964.


## Contact
Feel free to contact us at <schatterjee@niu.edu> and/or <gg0658@wayne.edu>
