n = 144
p = 10
R = 3

# Simulate coordinates on a grid
coords = as.matrix(expand.grid(1:sqrt(n),1:sqrt(n))) + matrix(runif(2*n,-.05,.05),n,2)

# Choose centroids using K-means clustering
fit = kmeans(coords,centers=R)
centroids.true = fit$centers
clusters.true = fit$cluster

# Use the true matrix W to choose phi.beta 
# Set sigma2.beta to make unit diagonal on covariance for beta
W.true = 1/as.matrix(dist(centroids.true)); diag(W.true) = 0
L = eigen(W.true)$values
phi.beta.true = .9*(1/max(L))
sigma2.beta.true = 1/solve(diag(1,R)-phi.beta.true*W.true)[1,1]

# Specify theta based on true centroids
theta.true = .85*(-log(1/.95-1)/mean(rowSums(W.true)))

# Simulate Data
data = Simulate.Data.Probit(n=n,p=p,R=R,phi=rep(5,R),coords=coords,beta.mean=0,sigma2.beta=sigma2.beta.true,
                            phi.beta=phi.beta.true,centroids=centroids.true,var.select=FALSE,theta=theta.true,pj=.5)

# Unpack objects from data
z.all = data$z.all
X.all = data$X.all
y.all = data$y.all
beta.true = data$beta
mu.true = data$mu
gamma.true = data$gamma
tau2.true = data$tau2
g.true = data$g
phi.true = data$phi
S.list.true = data$S.list
S.inv.list.true = data$S.inv.list.true
theta.true = data$theta
D = data$D
K = data$K
CHs.true = Calc.CHs(coords,clusters.true)
CH.inds.true = Calc.CH.inds(coords,clusters.true,CHs.true)

# Split data into training and testing
holdout = .2
train.ind = sort(sample(1:n,size=round((1-holdout)*n,0),replace=FALSE))
test.ind = c(1:n)[!1:n%in%train.ind]

z.train = z.all[train.ind,]
y.train = y.all[train.ind,]
X.train = X.all[train.ind,]
coords.train = coords[train.ind,]
D.train = D[train.ind,train.ind]

z.test = z.all[-train.ind,]
y.test = y.all[-train.ind,]
X.test = X.all[-train.ind,]
coords.test = coords[-train.ind,]
D.test = D[-train.ind,-train.ind]

# Initialize with random clusters
initial.centroids = matrix(randomLHS(R,2),R,2)*(sqrt(n)-1)+1
clusters.curr = apply(distance(coords.train,initial.centroids),1,which.min)
CHs.curr = Calc.CHs(coords.train,clusters.curr)
CH.inds.curr = Calc.CH.inds(coords.train,clusters.curr,CHs.curr)

# Plot true clusters vs. Initial clusters
par(mfrow=c(2,2))

plot(NA,NA,xlim=c(1,sqrt(n)),ylim=c(1,sqrt(n)),main="True")
text(coords[,1],coords[,2],labels=1:n,cex=.75,col=c("black","red")[CH.inds.true+1])
lapply(CHs.true,function(c) polygon(c,lty=2))
text(centroids.true[,1],centroids.true[,2],1:R,col="blue",cex=1.25)

plot(NA,NA,xlim=c(1,sqrt(n)),ylim=c(1,sqrt(n)),main="Initial")
text(coords.train[,1],coords.train[,2],labels=train.ind,cex=.75,col=c("black","red")[CH.inds.curr+1])
lapply(CHs.curr,function(c) polygon(c,lty=2))

# Create storage matrices/arrays
I = 20000
y.store = matrix(NA,I,length(y.train)); y.store[1,] = z.train*2-1
beta.store = array(NA,c(I,p,R)); beta.store[1,,] = 1
gamma.store = array(NA,c(I,p,R)); gamma.store[1,,] = rep(1,p)
theta.store = rep(NA,I); theta.store[1] = 0
mu.store = matrix(NA,I,R); mu.store[1,] = 1
phi.store = matrix(NA,I,R); phi.store[1,] = rep(.1,R)
phi.beta.store = rep(NA,I); phi.beta.store[1] = 0
sigma2.beta.store = rep(NA,I); sigma2.beta.store[1] = 1
y.preds = matrix(NA,I,length(y.test)); y.preds[1,] = 0

y.curr = matrix(y.store[1,],ncol=1)
mu.curr = mu.store[1,]
beta.curr = beta.store[1,,]
gamma.curr = gamma.store[1,,]
theta.curr = theta.store[1]
phi.curr = phi.store[1,]
phi.beta.curr = phi.store[1]
sigma2.beta.curr = sigma2.beta.store[1]
W.curr = 1/as.matrix(dist(initial.centroids)); diag(W.curr) = 0
K.inv.curr = (diag(1,R) - phi.beta.curr*W.curr)/sigma2.beta.curr
K.curr = solve(K.inv.curr)
# probs = t( (1 + exp(theta.curr*W.curr%*%t((gamma.curr==0)-(gamma.curr==1))))^(-1) )
probs = matrix(1,p,R)

S.list.curr = list()
S.inv.list.curr = list()

for(r in 1:R)
{
  ind = which(clusters.curr==r)
  S.list.curr[[r]] = (exp(-D.train[ind,ind]/phi.curr[r]) + diag(1e-8,length(ind)))
  S.inv.list.curr[[r]] = solve(S.list.curr[[r]])
}

# Run MCMC on this bitch
for(i in 2:I)
{
  # y.curr = y.train
  y.curr = Latent.Draw(X.all=X.train,y.curr=y.curr,z.all=z.train,clusters=clusters.curr,beta.mat=beta.curr,
                       mu.vec=mu.curr,Sigma.inv.list=S.inv.list.curr)
  y.curr = matrix(y.curr,ncol=1)
  
  draw = Sample.phi.beta(mu=mu.curr,beta=beta.curr,gamma=gamma.curr,phi.beta.curr=phi.beta.curr,step=.05,
                         K.inv.curr=K.inv.curr,W=W.curr,sigma2.beta=sigma2.beta.curr)
  phi.beta.curr = draw$phi.beta.draw
  K.inv.curr = draw$K.inv.draw
  K.curr = solve(K.inv.curr)
  
  draw = Sample.sigma2.beta(mu=mu.curr,beta=beta.curr,gamma=gamma.curr,K.inv.curr=K.inv.curr,sigma2.beta.curr=sigma2.beta.curr)
  sigma2.beta.curr = draw$sigma2.beta.draw
  K.inv.curr = draw$K.inv.draw
  K.curr = solve(K.inv.curr)
  
  draw = GP.Phi.Sample.Probit(X.all=X.train,y.all=y.curr,clusters=clusters.curr,D=D.train,step=1,phi.UB=Inf,mu=mu.curr,
                              beta=beta.curr,phi.curr=phi.curr,S.list.curr=S.list.curr,S.inv.list.curr=S.inv.list.curr)
  phi.curr = draw$phi.draw
  S.list.curr = draw$S.list.draw
  S.inv.list.curr =  draw$S.inv.list.draw
  
  draw = Sample.Mu(y.all=y.curr,X.all=X.train,S.inv.list=S.inv.list.curr,K=K.curr,
                   clusters=clusters.curr,beta=beta.curr,mu.curr=mu.curr)
  mu.curr = draw
  
  draw = SSVS.Gibbs(y.all=y.curr,X.all=X.train,clusters=clusters.curr,S.inv.list=S.inv.list.curr,K=K.curr,
                    mu=mu.curr,beta.curr=beta.curr,gamma.curr=gamma.curr,probs=probs)
  gamma.curr = draw$gamma.draw
  beta.curr = draw$beta.draw
  
  draw = Cluster.Assign.CH(X.all=X.train,y.all=y.curr,clusters.curr=clusters.curr,beta=beta.curr,mu=mu.curr,
                           gamma=gamma.curr,S.list=S.list.curr,S.inv.list=S.inv.list.curr,tau2=rep(1,R),g=rep(1e-8,R),
                           phi=phi.curr,D=D.train,coords=coords.train,CHs.curr=CHs.curr,CH.inds.curr=CH.inds.curr)
  S.list.curr = draw$S.list.draw
  S.inv.list.curr = draw$S.inv.list.draw
  clusters.curr = draw$clusters.draw
  CHs.curr = draw$CHs.draw
  CH.inds.curr = draw$CH.inds.draw
  W.curr = 1/as.matrix(dist(Calc.Centroids(clusters.curr,coords.train))); diag(W.curr) = 0
  
  # UHOH - if my cluster assignments change, W changes, in a way that my current value for phi.beta does not allow!
  phi.beta.curr = min(.98*(1/max(eigen(W.curr)$values)),phi.beta.curr)
  K.inv.curr = (diag(1,R) - phi.beta.curr*W.curr)/sigma2.beta.curr
  K.curr = solve(K.inv.curr)
  
  # theta.max = -log(1/.99-1)/max(rowSums(W.curr))
  # bound.prob = 1 - .05
  # rate.seq = seq(.Machine$double.eps,10,length=1000)
  # rate.opt = rate.seq[which.min(abs(pgamma(theta.max,shape=1,rate=rate.seq)-bound.prob))]
  # draw = Sample.theta(theta.curr=theta.curr,gamma=gamma.curr,W=W.curr,a=1,b=rate.opt,theta.UB=100,step.size=.5)
  # theta.curr = draw
  # probs = t( (1 + exp(theta.curr*W.curr%*%t((gamma.curr==0)-(gamma.curr==1))))^(-1) )
  
  for(j in 1:length(test.ind))
  {
     pt = test.ind[j]
     y.preds[i,j] = Predict.ynew(x.new=X.all[pt,],coords.new=coords[pt,],coords=coords.train,clusters=clusters.curr,
                                 S.inv.list=S.inv.list.curr,beta=beta.curr,mu=mu.curr,X=X.train,y=y.curr,phi=phi.curr,
                                 CHs=CHs.curr)
  }

  y.store[i,] = y.curr
  gamma.store[i,,] = gamma.curr
  beta.store[i,,] = beta.curr
  mu.store[i,] = mu.curr
  theta.store[i] = theta.curr
  phi.store[i,] = phi.curr
  phi.beta.store[i] = phi.beta.curr
  sigma2.beta.store[i] = sigma2.beta.curr
  
  if(i%%100==0)
  {
    plot(NA,NA,xlim=c(1,sqrt(n)),ylim=c(1,sqrt(n)),main="True")
    # text(coords[,1],coords[,2],labels=1:n,cex=.75,col=c("black","red")[CH.inds.true+1])
    lapply(CHs.true,function(c) polygon(c,lty=2))
    
    plot(NA,NA,xlim=c(1,sqrt(n)),ylim=c(1,sqrt(n)),main=paste("Iteration",i))
    # text(coords[,1],coords[,2],labels=1:n,col=c("black","red")[CH.inds.curr+1],cex=.75)
    lapply(CHs.curr,function(c) polygon(c,lty=2))
  }
  print(i)
}

LBs = apply(y.preds[-c(1:500),],2,function(x)quantile(x,.025))
UBs = apply(y.preds[-c(1:500),],2,function(x)quantile(x,.975))
Means = colMeans(y.preds[-c(1:500),])
sum((y.test>LBs)&(y.test<UBs))/length(y.test)

par(mfrow=c(3,2))
plot(phi.beta.store,type="l")
abline(h=phi.beta.true,col="chartreuse")
plot(sigma2.beta.store,type="l")
abline(h=sigma2.beta.true,col="chartreuse")
plot(theta.store,type="l")
abline(h=theta.true,col="chartreuse")
for(r in 1:R)
{
  plot(phi.store[,r],type="l")
  abline(h=phi.true[r],col="chartreuse")
}
plot(beta.true,apply(gamma.store[500:(i-1),,],3,function(c) colMeans(c)))

