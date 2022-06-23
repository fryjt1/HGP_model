library(parallel)
library(foreach)
library(doParallel)
library(truncnorm)
library(laGP)
library(mvtnorm)

SSVS.run = function(I,data,Ising=TRUE)
{
  mu.store = matrix(NA,I,R); mu.store[1,] = data$mu
  beta.store = array(NA,c(I,p,R)); beta.store[1,,] = data$beta
  gamma.store = array(NA,c(I,p,R)); gamma.store[1,,] = data$gamma
  theta.store = rep(NA,I); theta.store[1] = data$theta
  
  mu.curr = mu.store[1,]
  beta.curr = beta.store[1,,]
  gamma.curr = gamma.store[1,,]
  theta.curr = theta.store[1]
  K.curr = data$K

  eps = .Machine$double.eps
  theta.max = -log(1/.99-1)/max(rowSums(W))
  bound.prob = 1 - .01
  rate.seq = seq(eps,10,length=1000)
  rate.opt = rate.seq[which.min(abs(pgamma(theta.max,shape=1,rate=rate.seq)-bound.prob))]
  
  for(i in 2:I)
  {
    probs = t( (1 + exp(theta.curr*W%*%t((gamma.curr==0)-(gamma.curr==1))))^(-1) )
    draw = Sample.Mu(y.all=data$y.all,X.all=data$X.all,S.inv.list=data$S.inv.list,K=K.curr,
                     clusters=data$clusters,beta=beta.curr,mu.curr=mu.curr)
    mu.curr = draw
    
    draw = SSVS.Gibbs(y.all=data$y.all,X.all=data$X.all,clusters=data$clusters,S.inv.list=data$S.inv.list,K=K.curr,
                      mu=mu.curr,beta.curr=beta.curr,gamma.curr=gamma.curr,probs=probs)
    gamma.curr = draw$gamma.draw
    beta.curr = draw$beta.draw
    
    if(Ising==TRUE)
    {
      draw = Sample.theta(theta.curr=theta.curr,gamma=gamma.curr,W=W,a=1.5,b=rate.opt,theta.UB=100,step.size=.08)
      theta.curr = draw
    }else{
      theta.curr = 0
    }
    
    theta.store[i] = theta.curr
    beta.store[i,,] = beta.curr
    mu.store[i,] = mu.curr
    gamma.store[i,,] = gamma.curr
    print(i)
  }
  
  burn = round(.1*I,0)
  mu.means = colMeans(mu.store[-c(1:burn),])
  beta.means = apply(beta.store[-c(1:burn),,],3,colMeans)
  gamma.means = apply(gamma.store[-c(1:burn),,],3,colMeans)
  return(gamma.means)
}

B = 100
n.seq = c(625)                  # Changing sample size
p.seq = c(5,15)                    # Changing number of covariates
R.seq = c(2,4,8)                # Changing number of clusters
theta.seq = rep(0,3)    # Changing proportion of nonzero covariates
sigma2.seq = c(1)               # Change noise
beta.seq = c(0,2)             # Changing the mean of nonzero betas
storage = matrix(NA,2*B*length(n.seq)*length(p.seq)*length(R.seq)*
                   length(theta.seq)*length(sigma2.seq)*length(beta.seq),8)
system.time({
count = 1
for(n in n.seq)
{
  # Create a grid for the coordinates
  coords.seq = 1:sqrt(n)
  coords = as.matrix(expand.grid(coords.seq,coords.seq))
  
  for(p in p.seq)
  {
    for(sigma2 in sigma2.seq)
    {
      for(R in R.seq)
      {
        # Use K-means to choose centroids
        fit = kmeans(coords,centers=R)
        centroids = fit$centers
        D.centroids = sqrt(distance(centroids))
        W = 1/D.centroids; diag(W) = 0
        theta.seq = seq(0,-log(1/.95-1)/mean(rowSums(W)),length=length(theta.seq))
        
        for(theta in theta.seq)
        {
          for(beta.mean in beta.seq)
          {
            cl = makeCluster(10)
            registerDoParallel(cl)
            runs = foreach(i=1:B,.combine='rbind',.packages=c("laGP","truncnorm","mvtnorm","foreach")) %dopar% 
            {
              data = Simulate.Data(n=n,p=p,R=R,tau2=rep(sigma2,R),phi=rep(1e-10,R),g=rep(0,R),
                                   beta.mean=beta.mean,sigma2.beta=1,phi.beta=0,centroids=centroids,
                                   var.select=TRUE,theta=theta,pj=.5)
              
              run.Ising = SSVS.run(I=5000,data=data,Ising=TRUE)
              run.ind = SSVS.run(I=5000,data=data,Ising=FALSE)
              AFCCF.Ising = sum(data$gamma*run.Ising + (1-data$gamma)*(1-run.Ising))/(R*p)
              AFCCF.ind = sum(data$gamma*run.ind + (1-data$gamma)*(1-run.ind))/(R*p)
              
              return(c(AFCCF.Ising,AFCCF.ind))
              }
            stopCluster(cl)
            
            storage[count:(count+2*B-1),1:6] = matrix(rep(c(n,p,sigma2,R,round(theta,2),beta.mean),2*B),ncol=6,byrow=TRUE)
            storage[count:(count+2*B-1),7] = rep(c(1,2),B) 
            storage[count:(count+2*B-1),8] = c(t(runs))
            count = count + 2*B

            print(paste("n =",n,"p =",p,"sigma2 =",sigma2,"R =",R,
                          "theta =",round(theta,2),"beta =",beta.mean))
          }
        }
      }
    }
  }
}

})
write.csv(storage, file = "~/Dropbox/Leman Research/SCOMS/Ising-SSVS/Ising_SSVS_sim_bigN.csv")
save.image(file="~/Dropbox/Leman Research/SCOMS/Ising-SSVS/Ising_SSVS_sim_bigN.rdata")
# 
# library(ggplot2)
# colnames(storage) = c("n","p","pj","sigma2","R","phi.beta","beta.mean","method","AFCCF")
# df = as.data.frame(storage)
# 
# p1 = ggplot(aes(y = AFCCF, x=factor(phi.beta), fill=factor(method)), data = df[df$pj==.1,]) +
#   geom_boxplot() + 
#   theme(legend.position="none",panel.background = element_blank(),
#         panel.border = element_rect(colour = "black",fill=NA, size=1))
# print(p1)
# 
# p2 = ggplot(aes(y = AFCCF, x=factor(phi.beta), fill=factor(method)), data = df[df$pj==.5,]) +
#   geom_boxplot() + 
#   theme(legend.position="none",panel.background = element_blank(),
#         panel.border = element_rect(colour = "black",fill=NA, size=1))
# print(p2)
# 
# p3 = ggplot(aes(y = AFCCF, x=factor(phi.beta), fill=factor(method)), data = df[df$pj==.9,]) +
#   geom_boxplot() + 
#   theme(legend.position="none",panel.background = element_blank(),
#         panel.border = element_rect(colour = "black",fill=NA, size=1))
# print(p3)


