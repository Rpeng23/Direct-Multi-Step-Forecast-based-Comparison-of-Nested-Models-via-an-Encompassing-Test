WB_HCPI_Q<-read.csv("C:/Users/rp1e23/Desktop/WB_HCPI_Q.csv",header= TRUE)

q0=4
h=4
q=q0
m = max(h,q)+2
n=213

usa_cpi_q = WB_HCPI_Q$USA_HCPIQ
globadata = WB_HCPI_Q

globadata = globadata[,-1]

Y_usa = 100*(log(usa_cpi_q[m:n])-log(usa_cpi_q[(m-h):(n-h)]))

globalcpi = matrix(0,n-m+1,23)

for (i in 1:23) {
  globalcpi[,i] = 100*(log(globadata[m:n,i])-log(globadata[(m-h):(n-h),i]))
}


Y_usa <- as.matrix(Y_usa)

nx=dim(Y_usa)[1]
nc=dim(Y_usa)[2] 



XMAT_usa = matrix(0,nx,(q+1))
XMAT_global = matrix(0,nx,(q+1))

for (i in 1:23) {
  for (k in 0:q0) {
    XMAT_global[,k+1]= XMAT_global[,k+1]+400*(log(globadata[(m-k):(n-k),i])-log(globadata[(m-k-1):(n-k-1),i]))
  }
}

XMAT_global = XMAT_global/23


for (k in 0:q0){
XMAT_usa[,k+1]=400*(log(usa_cpi_q[(m-k):(n-k)])-log(usa_cpi_q[(m-k-1):(n-k-1)]))

}

ehat_usa = matrix(0, nrow = nx-h-round(0.25*nx)+1,ncol = q+1)

for (k in 0:q){
ehat_usa[,(k+1)]=recursive_hstep_fast(Y_usa,XMAT_usa[,1:(k+1)],0.25,h)
}

ehat_usa <- as.matrix(ehat_usa)

nehat=dim(ehat_usa)[1]
cehat=dim(ehat_usa)[2] 


cov_ehat_usa = t(ehat_usa)%*%ehat_usa/nehat


qhat_usa = which.min(log(diag(cov_ehat_usa))+((log(nehat)/nehat)*t(1:(q+1))))


#Encompassing using the above lags

zehat1_usa=recursive_hstep_fast(Y_usa,XMAT_usa[,1:qhat_usa],0.25,h)
zehat2_usa=recursive_hstep_fast(Y_usa,cbind(XMAT_usa[,1:qhat_usa],XMAT_global[,1:4]),0.25,h)


mu_vec = c(0.40,0.45)
nmu=length(mu_vec)

stat_usa = list()


for (m in 1:nmu){
stat_usa[[m]] = pred_encompass_dnorm(zehat1_usa,zehat2_usa,mu_vec[m])

}

pval_usa1 = stat_usa[[1]]$Pval_T1 #0.6868066
pval_usa2 = stat_usa[[2]]$Pval_T1 #0.964937

rmserat_usa = (sqrt(t(zehat2_usa)%*%zehat2_usa/nehat)/sqrt(t(zehat1_usa)%*%zehat1_usa/nehat)) #1.10368

