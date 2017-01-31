# Gooogle (Group Regularization for Zero Inflated Count Regression Models)

## Introduction
Zero inflated count data are common in many fields including health care research, insurance industries etc. which are modeled by ZIP/ZINB regression models. When the associated features/covariates possess an inherent grouping structure (statistically/mechanistically correlated) the traditional variable selection approaches are known to perform poorly. Thus they have to be extended to the grouped predictors selection in order to obtain sparse group solution along with identifying the important variables at both the group and individual levels. For group variable selection in ZIP/ZINB models we consider different commonly used group regularizations such as group LASSO, group SCAD and other group level regularization methods and group bridge, sparse group LASSO and GEL which perform bi-evel selection. The ZIP/ZINB model includes a logistic component to model the presence of excess zeros and a Poisson/negative Binomial component to model the count data. For both the models, log link function is used to  regress the mean of the count model on a set of predictors while the mixture proportion parameter is regresed on another (or same as the count part) set of covariates via the logistic regression.  

The details of the statistical model are as follows:
<img src="misc/model.png" width="600" align="center">

Using LSA approximation on the ZIP/ZINB likelihood we obtain regularized estimates of the regression parameters for zero and count models. The tuning parameter for the final model corresponds to the minimum AIC/BIC values.  

You can install our Gooogle package from Github
```r
install.packages("devtools")
devtools::install_github("/Gooogle")
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
- **penalty**: the penalty to be applied for regularization. For group selection, it is one of grLasso, grMCP or grSCAD while for bi-level selection it is gBridge, gel, cMCP or SGL.  

For greatest efficiency and least ambiguity, it is best if group is a vector of consecutive integers. If any coefficients are to be included in the model without being penalized, their grouping index should be zero. 

The gooogle function will return a list containing the following objects:
- **coefficients**: a list containing the estimates for the count and zero model.  
- **aic**: AIC for the fitted model.  
- **bic**: BIC for the fitted model.  
- **loglik**: Log-likelihood for the fitted model.

## Examples

#### Simulated data  
We construct a simulation function similar to Huang et al. (2009) [ https://doi.org/10.1093/biomet/asp020] which can be used to simulate a dataset according to the zero-inflated Poisson model. We use 40 covariates in 5 groups with 8 covariates in each group. For this example, we assume the covariates in the count part (X) and in the zero part (Z) to be the same (i.e set Z=X).

```r
    data.sim<-data.func.sim1(n=500,p=40,ngrp=5,rho=0.4,family="poisson",
    beta=c(1, -1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24)),
    gamma=c(-1,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24)))
```

The simulation function returns a dataset with X (the predictor matrix corresponding to count model), Z (the predictor matrix corresponding to zero model) and the outcome with zero abundance (Y) simulated according to the above parameter settings. 

We can do both group level and bi-level selection in zero inflated count data using gooogle function. If we want to do a group level selection we can use any of the penalties among grLasso, grMCP or grSCAD. Below is an example of using gooogle in the simulated dataset using grLasso.

```r
yvar<-data.sim$yvar
xvars<-data.sim$xvars
zvars<-data.sim$zvars
group<-gamma<-c(-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
fit.gooogle <- gooogle(data=data.sim,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="poisson",penalty="grLasso")
fit.gooogle
```
Similarly we can do a bi-level selection on the simulated data using gBridge penalty in the gooogle function.

```r
fit.gooogle <- gooogle(data=data.sim,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="poisson",penalty="gBridge")
fit.gooogle
```

#### Real data  
Let's try one example on the real data. I am using an auto insurance claim dataset used by Qian et al. (2016) (http://www.tandfonline.com/doi/full/10.1080/10618600.2015.1005213). The scaled dataset along with the grouping indices is available in their HDtweedie package in R. The dataset contains the aggregate claim loss of an auto insurance policy as the response and 21 predictors of which 11 are continuous and 10 are categorical. 

Type ```auto```.
```r
# claim loss	CAR_TYPE_2	CAR_TYPE_3	CAR_TYPE_4	CAR_TYPE_5	CAR_TYPE_6	MAX_EDUC_2	MAX_EDUC_3	MAX_EDUC_4	MAX_EDUC_5	KIDSDRIV	KIDSDRIV2	KIDSDRIV3	TRAVTIME
0	0	1	0	0	0	0	1	0	0	0	-0.166666667	0	-0.748425538
0	0	0	0	0	1	1	0	0	0	0	-0.166666667	0	-0.494449479
19	0	0	0	1	0	0	0	0	1	0	-0.166666667	0	0.140490669
0	0	0	1	0	0	1	0	0	0	0	-0.166666667	0	0.775430816
0	0	0	0	1	0	0	0	0	0	0	-0.166666667	0	-0.049991376
2	0	0	0	0	1	1	0	0	0	1	0.333333333	0.2	0.648442787
0	0	0	0	1	0	1	0	0	0	0	-0.166666667	0	0.013502639
0	0	0	1	0	0	0	0	1	0	0	-0.166666667	0	-1.192883642
0	0	0	0	1	0	0	0	0	0	0	-0.166666667	0	0.902418846
45	0	1	0	0	0	0	1	0	0	0	-0.166666667	0	0.711936802
0	0	0	0	0	0	0	0	0	0	0	-0.166666667	0	1.156394905
0	0	0	1	0	0	0	1	0	0	0	-0.166666667	0	-1.129389627
0	0	0	0	1	0	0	1	0	0	0	-0.166666667	0	-0.494449479
0	0	1	0	0	0	1	0	0	0	0	-0.166666667	0	-0.24047342
0	0	1	0	0	0	0	0	0	1	0	-0.166666667	0	1.981817097
0	0	1	0	0	0	0	0	0	1	0	-0.166666667	0	-1.002401597
0	0	0	0	1	0	0	1	0	0	0	-0.166666667	0	0.838924831
0	0	0	0	1	0	0	0	0	1	0	-0.166666667	0	-1.510353716
0	0	0	0	1	0	0	1	0	0	0	-0.166666667	0	0.902418846
0	1	0	0	0	0	1	0	0	0	0	-0.166666667	0	-0.748425538
0	0	0	1	0	0	1	0	0	0	0	-0.166666667	0	0.394466728
0	0	1	0	0	0	0	0	1	0	0	-0.166666667	0	-1.510353716
0	0	1	0	0	0	0	1	0	0	0	-0.166666667	0	-0.621437509
0	0	1	0	0	0	0	0	0	0	0	-0.166666667	0	-0.303967435
8	0	0	0	1	0	0	1	0	0	0	-0.166666667	0	0.203984683
9	0	0	1	0	0	0	0	1	0	2	1.833333333	3.4	-0.303967435
4	0	0	0	1	0	0	0	0	1	0	-0.166666667	0	-0.36746145
0	0	0	0	1	0	0	0	0	1	0	-0.166666667	0	-1.446859701
0	0	1	0	0	0	1	0	0	0	1	0.333333333	0.2	-0.557943494
8	0	0	0	1	0	1	0	0	0	0	-0.166666667	0	-0.430955464
...
...
```

Note that for each of the 11 continuous predictors the polynomials (up to order 3)  have been created to treat them as single group. For the categorical variables (more than 2 levels) dummy variables are created according to the number of levels of each.     


We can run the gooogle function to this real data to fit a ZIP model using group bridge regularization. 

```r
fit.poisson<-gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,samegrp.overlap=T,crit="BIC"
dist="poisson",penalty="gBridge")
fit.poisson
```
The final model is selected based on the minimum BIC value. From the output we see that out of 56 variables only 8 variables distributed across 4 groups (JOBCLASS, MVRPTS, REVOLKED and  AREA) for count model and 5 variables distributed across 3 groups (JOBCLASS, MVRPTS and AREA) for the zero model have been selected in the final model. The output also gives the minimum BIC value and the corresponding AIC and the log likelihood for the final fitted model.

## Citation


## Contact
Feel free to contact us at <schatterjee@niu.edu> and/or <gg0658@wayne.edu>
