rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(319)

#basic settings
nbehavior=5
ncat.data=c(15,15)
ndata.types=length(ncat.data)

#get parameters
phi=list()
for (i in 1:ndata.types){
  tmp=matrix(NA,nbehavior,ncat.data[i])
  seq1=1:ncat.data[i]
  mus=seq(from=2,to=ncat.data[i],length.out=nbehavior)
  diff1=mus[2]-mus[1]
  for (j in 1:nbehavior){
    aux=dnorm(seq1,mean=mus[j],sd=diff1/4)
    tmp[j,]=aux/sum(aux)
  }
  phi[[i]]=tmp
  image(phi[[i]])
}
theta.true=theta=rep(1/nbehavior,nbehavior)
phi.true=phi

#look at these parameters
par(mfrow=c(ceiling(nbehavior/2),2),mar=rep(1,4))
for (i in 1:nbehavior) plot(phi.true[[1]][i,],type='h',main=i)
round(apply(phi.true[[1]],2,max),3)

par(mfrow=c(ceiling(nbehavior/2),2),mar=rep(1,4))
for (i in 1:nbehavior) plot(phi.true[[2]][i,],type='h',main=i)
round(apply(phi.true[[2]],2,max),3)

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
