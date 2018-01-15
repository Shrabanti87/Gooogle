# Gooogle (Group Regularization for Zero Inflated Count Regression Models)

## Introduction
Zero-inflated count data are omnipresent in many fields including health care research and actuarial science. Zero-inflated Poisson (ZIP) and Zero-inflated Negative Binomial (ZINB) regression are commonly used to model these outcomes. However, when the features to be associated possess an inherent grouping structure, traditional variable selection approaches are known to produce nonsensical results. In order to be able to perform group variable selection in ZIP/ZINB models, we extend various commonly used group regularizations such as group LASSO and group bridge to ZIP/ZINB models. These mixture models typically include a logistic component to model the presence of excess zeros and a Poisson/negative Binomial component to model the count data. 

The details of the statistical model are as follows:

<img src="misc/github.png" width="600" align="center">

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

## Citation

Huang, J., S. Ma, H. Xie, and C. Zhang (2009). A group bridge approach for variable selection. Biometrika 96(2), 339–355.

Park, C. and Y. Yoon (2011). Bridge regression: adaptivity and group selection. Journal of
Statistical Planning and Inference 141 (11), 3506{3519.

Jochmann, M. (2013). What belongs where? variable selection for zero-in
ated count models with an application to the demand for health care. Computational Statistics 28, 1947{1964.


## Contact
Feel free to contact us at <schatterjee@niu.edu> and/or <gg0658@wayne.edu>
