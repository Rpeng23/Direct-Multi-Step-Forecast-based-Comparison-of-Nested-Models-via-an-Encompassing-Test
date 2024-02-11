#' Direct Multi-Step Forecast based Comparison of Nested Models via an Encompassing Test
#'
#' @param  e1hat: n by 1 vector of out of sample forecast errors from model 1
#' @param  e2hat: n by 1 vector of out of sample forecast errors from model 2
#' @param  mu0:   Fraction of the sample, which should be within 0 and 1. (0.5 should be avoid)
#' @return A list of normalised d statistics and corresponding P values will be produced.
#' @references Pitarakis, J. Y. (2023). Direct Multi-Step Forecast based Comparison of Nested Models via an Encompassing Test. arXiv preprint arXiv:2312.16099.
#' @references Andrews, D. W. (1991). Heteroskedasticity and autocorrelation consistent covariance matrix estimation. Econometrica: Journal of the Econometric Society, 817-858.
#' @examples
#' e1<- rnorm(15);
#' e2<- rnorm(15);
#' temp1 <- pred_encompass_dnorm(e1,e2,mu0=0.2)
#' temp1$T1_d_alrv     #normalised d statistics for raw data using long-run variance with Andrews quadratic spectral kernel HAC
#' temp1$Pval_T1_d_alrv  #P value of it
#' @export

pred_encompass_dnorm=function(e1hat,e2hat,mu0){
  
  

  #Stopping criteria
  if(is.null(e1hat)){
    stop("e1hat must be a numeric series")
  }else if(!is.numeric(e1hat)){
    stop("e1hat must be a numeric series")
  }
  
  if(is.null(e2hat)){
    stop("e2hat must be a numeric series")
  }else if(!is.numeric(e2hat)){
    stop("e2hat must be a numeric series")
  }
  
   e1hat <- as.matrix(e1hat)
   e2hat <- as.matrix(e2hat)
  
  if(dim(e1hat)[1] != dim(e2hat)[1]){
    stop("e1hat and e2hat must have same length")
  }else if(dim(e1hat)[2]>1 | dim(e2hat)[2]>1){
    stop("e1hat and e2hat must be one dimension")
  }
  
  if(is.null(mu0)){
    stop("mu0 must be between 0 and 1")
  }else if(!is.numeric(mu0)){
    stop("mu0 must be between 0 and 1")
  }
  
   if(mu0 == 0.5 ){
     stop("mu0 cannot equal to 0.5")
   }
   
  
  #e1hat <- as.matrix(e1hat)
  n=dim(e1hat)[1]
  
  time_vec = c(1:n)
  m0 = round(n*mu0)
  w = (time_vec<=m0)
  
  e1sq = e1hat^2
  e2sq = e2hat^2
  e12 = e1hat*e2hat
  sigsq2 = mean(e2sq)
  e2sq_demean = e2sq-sigsq2
  
  fm = ((1-2*mu0)^2)/(4*mu0*(1-mu0))
  nw_lags = min(floor(1.2*n^(1/3)),n)
  
  #d= e1sq - 0.5*((1/mu0)*e12*w+(1/(1-mu0))*e12*(1-w))
  d = (e1sq-sigsq2) - 0.5*((1/mu0)*(e12-sigsq2)*w+(1/(1-mu0))*(e12-sigsq2)*(1-w))
  
  sigsq_d = t(d-mean(d))%*%(d-mean(d))/n
  sigsq_d_nw = NW_lrv(d,nw_lags,1)
  sigsq_d_alrv = andrews_lrv(d)
  
  phihatsq = t(e2sq_demean)%*%(e2sq_demean)/n
  phihatsq_nw = NW_lrv(e2sq_demean,nw_lags,1)
  phihatsq_alrv = andrews_lrv(e2sq_demean)
  
  T1 = sqrt(n)*mean(d)/sqrt(fm*phihatsq)
  T1_nw = sqrt(n)*mean(d)/sqrt(fm*phihatsq_nw)
  T1_alrv = sqrt(n)*mean(d)/sqrt(fm*phihatsq_alrv)
  
  T1_d = sqrt(n)*mean(d)/sqrt(sigsq_d)
  T1_d_nw = sqrt(n)*mean(d)/sqrt(sigsq_d_nw)
  T1_d_alrv = sqrt(n)*mean(d)/sqrt(sigsq_d_alrv)
  
  Pval_T1 = pnorm(T1, lower.tail = FALSE)
  Pval_T1_nw = pnorm(T1_nw, lower.tail = FALSE)
  Pval_T1_alrv = pnorm(T1_alrv, lower.tail = FALSE)
  
  Pval_T1_d = pnorm(T1_d, lower.tail = FALSE)
  Pval_T1_d_nw = pnorm(T1_d_nw, lower.tail = FALSE)
  Pval_T1_d_alrv = pnorm(T1_d_alrv, lower.tail = FALSE)
  
  return(list(T1=T1,T1_nw=T1_nw,T1_alrv=T1_alrv,T1_d=T1_d,T1_d_nw=T1_d_nw,T1_d_alrv=T1_d_alrv, 
              Pval_T1=Pval_T1,Pval_T1_nw=Pval_T1_nw,Pval_T1_alrv=Pval_T1_alrv,
              Pval_T1_d=Pval_T1_d,Pval_T1_d_nw=Pval_T1_d_nw,Pval_T1_d_alrv=Pval_T1_d_alrv) )
  
}