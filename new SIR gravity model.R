numofcities = 20
set.seed(1)
T = 100
p = 0.9
gamma=0.2
c = 5
vmax=0.85
vbase=0.15
V = function(It,Nt) ((vmax-vbase)*(It/Nt))+vbase
phi=0
numofinfected=0
S = matrix(data = 0, nrow = T, ncol = numofcities)
S[1,] = 1000
S[1,1]=999
I = matrix(data = 0, nrow = T, ncol = numofcities)
I[1,1]=1
R = matrix(data = 0, nrow = T, ncol = numofcities)
N=matrix(data=1000,nrow=1,ncol=numofcities)
N[1,1] = 500
Nmax=max(N)
xvals=runif(n=numofcities,min=0,max=100)
yvals=runif(n=numofcities,min=0,max=100)

distances = function(x1,y1,x2,y2) sqrt((x2-x1)^2+(y2-y1)^2) 
D = matrix(data = 0, nrow = numofcities, ncol = numofcities)
for(i in 1:numofcities){
  for(j in 1:numofcities){
    D[i,j]=distances(xvals[i],yvals[i],xvals[j],yvals[j])
  }
}

beta = p*c
alpha = 2
E=matrix(data = 0,nrow=1,ncol=numofcities)
for(t in 1:(T-1)){
  E=exp(-beta*((1/(D+1)^alpha)%*%((I[t,]/N[1,])*(N[1,]/Nmax))))
  S[t+1,] = (S[t,]*E*(1-V(I[t,],N[1,]))) + (phi*R[t,])
  I[t+1,] = S[t,]*(1-E)*(1-V(I[t,],N[1,])) + I[t,]*(1-gamma)
  R[t+1,] = ((1-phi)*R[t,])+(I[t,]*gamma)+V(I[t,],N[1,])*S[t,]
}



mycol=colorRamp(c(rgb(0,0,1), rgb(1,0,0)))
time=1
maxfrac=max(I[time,]/N[1,])
cols=mycol(I[time,]/N[1,]/maxfrac)/255
plot(xvals,yvals,pch=21,bg=rgb(cols),cex=2)

#matplot(I,type="l",col=1,lty=1)

