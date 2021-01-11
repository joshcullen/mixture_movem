seq1=1:ngibbs
seq1=5000:ngibbs
plot(model1$loglikel[seq1],type='l')

compare=function(true1,estim1){
  rango=range(c(true1,estim1))
  plot(true1,estim1,ylim=rango,xlim=rango)
  lines(rango,rango,col='red')
}

par(mfrow=c(1,1))
theta.estim=model1$theta[ngibbs,]
# theta.estim=store.theta[ngibbs,]
plot(theta.estim,type='h')

#get z.estim (mode)
nobs=nrow(model1$z.posterior)
z.estim=rep(NA,nobs)
for (i in 1:nobs){
  ind=which(model1$z.posterior[i,]==max(model1$z.posterior[i,]))
  z.estim[i]=ind
}

k=data.frame(z.post=z.estim,z.map=model1$z.MAP)
table(k)

#compare with true
tmp=data.frame(zestim=z.estim,ztrue=z.true)
tmp1=table(tmp); tmp1
ordem=c(5,3,1,4,2)
tmp1[ordem,]

#look at phi's
par(mfrow=c(ceiling(ndata.types/2),2),mar=rep(1,4))
for (j in 1:ndata.types){
  phi1.estim=matrix(model1$phi[[j]][ngibbs,],nmaxclust,ncat.data[j])
  # phi1.estim=matrix(store.phi[[j]][ngibbs,],nmaxclust,ncat.data[j])
  compare(phi.true[[j]],phi1.estim[ordem,])
}

#look at gamma
gamma1=model1$gamma1
plot(table(gamma1),type='h')