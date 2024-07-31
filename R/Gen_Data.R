Gen_Data=function(nu,xi,R){
  m=length(R)
  w=u=v=t=vv=RR=numeric(m)
  w=runif(m,0,1)
  for(i in 1:m){
    RR[i]=0
    for(j in (m-i+1):m) RR[i]=RR[i]+R[j]
    v[i]=w[i]^(1/(i+RR[i]))
  }
  for(i in 1:m){
    vv[i]=1
    for(j in (m-i+1):m)
      vv[i]=vv[i]*v[j]
    u[i]=1-vv[i]
    t[i]=qnaka(u[i],scale=xi,shape=nu)
  }
  return(sort(t))
}
