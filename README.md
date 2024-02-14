# Multi-Step_Encompassing
This repository contains the codes for using the predictive accuracy comparison tests developed in Pitarakis, J. (2023). "Direct Multi-Step Forecast based Comparison
of Nested Models via an Encompassing Test", arXiv. DOI: https://doi.org/10.48550/arXiv.2312.16099.

# Background
Pitarakis propose a novel approach to testing whether multi-step forecasts from one model encompass the forecasts from a rival larger model when the two models have a nested structure.


Informally: consider two competing (nested) models given by Model 1 (smaller one) and Model 2 (larger one) and suppose that you wish to test the null that 
$E(MSE1(MSE1-MSE2))=0$. 
Once one obtains the out-of-sample forecast errors associated with each model (see below) one feeds them in the test statistic dbar which is shown to be asymptotically normal distributed. There is also an adjusted version of each of these test statistics discussed below.

The program builds forecasts in a recursive manner. The codes available are adapted for h-steps ahead prediction: 

$y_{t+h} = β_0 + βx_t + u_{t+h}$, (1)

and suppose that there are n observations. 

# The test statistics dBar
Suppose you have obtained the sequence of forecast errors associated with the two models. These are then used as inputs to obtain the numerator of the test statistics. Note that this statistic requires 1 ad-hoc inputs from the user ($mu_0$). So the code should expect these as inputs as well in addition to $\hat{e}_1$ and $\hat{e}_2$. 

It is usually advisable to use the residuals from the larger model to form the $\eta$ (i.e. use $\hat{e}_2$ instead of $\hat{e}_1$). 

In default setting, 3 conditions are considered, where the test statistic would be calculated under conditional homskedasticity and heteroskedasticity.

For the version of the test statistic that corrects for heteroskedasticity one needs to use a robust Newey-West type estimator as in Deng and Perron (2008) for instance.

Alternatively, one could consider adjust the long-run variance via Andrews (1991) approach. 

In addition, we consider both the situation where the error terms are demeaned or not.

To sum-up: given $\hat{e}_{1}$, $\hat{e}_{2}$, $mu_{0}$ the program should output the quantity $d_{bar}$ for in total 6 situation (demeaned or not, and 3 different types of long-run variance)

**Note that the corresponding p-values will be returned along with the test statistics, with the null hypothesis being the larger model does not add information in prediction.**

# Code Example

**Please first download and install the package *EncompassTest* from Github (via *EncompassTest_0.2.tar.gz*). The package is also available on CRAN.**

Suppose we have a dataset *data* (see the *dummy_data_2.xlsx* file in data folder for example) that contains 3 variables, $y$, $x_1$ and $x_2$ of a length of 250. 

For illustration, we consider 1 step ahead prediction and set $\pi_0 = 0.25$ (user could change it manually within 0 and 1 deciding what fraction of sample should be used). In this case, this means *round(250 x 0.25)* numbers of recursive residuals would be computed. **Note here, R would round up .5 to the nearest EVEN number, while MATLAB would push it far away from 0. Take our case for example, R would round it to 62, while MATLAB would round it to 63.**

`library(EncompassTest)`<br />
`ehat1 = recursive_hstep_fast(y,x1,pi0,1)`<br />
`ehat2 = recursive_hstep_fast(y,cbind(x1,x2),pi0,1)`<br />

This will give back two series of recursive errors with a length of $(250-62-1+1)=188$. 

`SEtest =  pred_encompass_dnorm(ehat1,ehat2,mu0=0.2)`<br />


**It is noted here, mu0 correpond to the user-define split-based sample used for the computation of test statistics. Normally speaking, all of them should be something within 0 and 1 and mu0=0.5 would reduce the test to the class sample mean, which suffers from variance degeneracy problem.** 


If we set mu0=0.2, running the codes above and we should be able to obtain results as follows:

| Statistics name        | Test           | P-values  |
| ---------------------- |:--------------:| ---------:|
| $T1$                   | 0.952          | 0.171     |
| $T1_{nw}$              | 0.864          | 0.194     |
| $T1_{alrv}$            | 0.952          | 0.171     |
| $T1^{d}$               | 1.25           | 0.106     |
| $T1^d_{nw}$            | 1.29           | 0.098     |
| $T1^d_{alrv}$          | 1.36           | 0.088     |

That is, based on our test, the the larger model do contain additional information that helps in prediction.

#' $T1^d_{alrv}$     #demeaned d statistics with Andrews quadratic kernel long-run variance
#' $T1_{nw}$     #d statistics with Newey-West type long run variance

