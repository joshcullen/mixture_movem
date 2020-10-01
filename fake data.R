rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(6)

#basic settings
nbehavior=6
ncat.data=c(8,8)
ndata.types=length(ncat.data)

#get parameters
phi=list()
for (i in 1:ndata.types){
  phi[[i]]=rdirichlet(nbehavior,alpha=rep(0.1,ncat.data[i]))
}
theta.true=theta=rep(1/nbehavior,nbehavior)
phi.true=phi

#look at these parameters
par(mfrow=c(ceiling(nbehavior/2),2),mar=rep(1,4))
for (i in 1:nbehavior) plot(phi.true[[1]][i,],type='h',main=i)
for (i in 1:nbehavior) plot(phi.true[[2]][i,],type='h',main=i)

#generate data
nobs=10000
ztmp=rmultinom(nobs,size=1,prob=theta)
z=rep(NA,nobs)
res=matrix(NA,nobs,2)
for (i in 1:nobs){
  z[i]=which(ztmp[,i]==1)
  y1=rmultinom(1,size=1,prob=phi[[1]][z[i],])
  y2=rmultinom(1,size=1,prob=phi[[2]][z[i],])
  res[i,]=c(which(y1==1),which(y2==1))
}
table(z)
colnames(res)=c('y1','y2')
z.true=z

for (i in 1:2){
  print(c(i,length(unique(res[,i]))))
}

#export fake data
setwd('U:\\GIT_models\\mixture_movem')
write.csv(res,'fake data.csv',row.names=F)
