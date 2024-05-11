# BANSMA: Bayesian Inference Framework for Multi-layered Mediation Analysis
In this project, we introduce a Bayesian inference framework for multi-layered mediation analysis, capable of handling continuous, binary, and polychotomous ordinal outcomes using probit models. From the perspective of mediation analysis, we decompose the joint effect into effects attributable to individual mediators within the framework of interventional mediation analysis. Our method, BANSMA, is carefully developed based on the [BANS](https://github.com/MinJinHa/BANS) framework.

## Setup
### Installing dependency packages
Since the BANSMA program is built upon the BANS framework, it is necessary to install packages related to BANS.
~~~
install.packages("rms")
install.packages("RcppArmadillo")
install.packages("RcppEigen")
install.packages(
    './BANS_package/BANS_1.0.tar.gz',
    repos = NULL,
    type = 'source'
)
~~~

## Functions in BANSMA

### `fit.Outcome`
The `fit.Outcome` function, situated under the known ordered layer, returns the fitting output that models the response variable located in the last layer using the remaining variables.

### `fit.Mediator`
Within the known ordered layer, the `fit.Mediator` function outputs the fitting that models the multi-layered relationships among the variables, excluding the response variable located in the last layer.

### `getIndEffect`
The `getIndEffect` function takes the outputs of `fit.Outcome` and `fit.Mediator` as inputs and returns an object containing the indirect effects. Through this object, access to the values of IE, IED, and IEC is facilitated.

An example usage is as follows:
~~~
bansma.out <- getIndEffect(palist=palist,chlist=chlist,fit.med=medfit,fit.out=outfit)
    
IE <- as.matrix(bansma.out$IE)
IED <- as.matrix(bansma.out$IED)
IEC <- as.matrix(bansma.out$IEC)
~~~
