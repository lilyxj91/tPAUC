
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tPAUC

<!-- badges: start -->
<!-- badges: end -->

The goal of tPAUC is to …

## Installation

You can install the development version of tPAUC like so:

``` r
# devtools::install_github("lilyxj91/tPAUC",force = TRUE)
```

## Example

``` r
library(tPAUC)
## basic example code
```

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tidyr)
library(fastglm)
#> Loading required package: bigmemory
dat <- simdata
n=nrow(dat)
Y0<-dat$Y
C0<-dat$delta
M0<-dat$M
VXS<-dat[c("xd","xc")]
datijp<-data_crossingdc(n,Y0,C0,M0,VXS,1)
#> Adding missing grouping variables: `k`

this_u = 0.4
ptm = Sys.time()
datijp<-datijp%>%filter(fp<=this_u)
YK = datijp$yk
XS<-cbind(int = 1,
          t5 = YK^(0.5),
          t6 = YK,
          t7 = YK^(2),
          xdi = datijp$xdk,
          xci = datijp$xck)
YS<-datijp$Ikj
rho01<-datijp$rho01
rho02<-datijp$rho02
rhoweight<-datijp$rhoweight
ordic=datijp$k-1
ordjc=datijp$j-1
m<-fastglm(XS,YS,weights = rhoweight,family=binomial(link="logit"),maxit=10000L)
#> Warning in eval(family$initialize): non-integer #successes in a binomial glm!
beta.hat = m$coefficients
L=cov_cal_long(beta.hat,M0,rho01,rho02,t(XS),ordic,ordjc)
V=solve(L$sigma1)%*%(L$sigma2)%*%solve(L$sigma1)

select=c(5:length(beta.hat))
se.store = diag(V)[select]
for(i in 1:length(select)){
  se <- sqrt(max(0.00000001,diag(V)[select[i]]))
  se.store[i] <- se
}
beta = beta.hat[select]
wald = beta/se.store
pvalue = 2*(1-pnorm(abs(wald)))
table<-round(as.matrix(cbind(
  beta,
  se.store,
  wald,
  pvalue)),3)
colnames(table)<-c("Estimate","SE","Wald","P-value")
table[,4]<-ifelse(table[,4]<0.001,"<0.001",table[,4])
table
#>     Estimate SE      Wald    P-value 
#> xdi "0.541"  "0.127" "4.259" "<0.001"
#> xci "1.122"  "0.213" "5.259" "<0.001"


print(Sys.time()-ptm)
#> Time difference of 35.73711 secs
```
