library(parallel)
library(foreach)
library(doParallel)

# Right now, if I don't use the GP, I'm fixing the alt. variance to be 1. Maybe I should estimate this...
SSVS.run = function(I,data,GP=TRUE)
{
  phi.beta.store = rep(NA,I); phi.beta.store[1] = data$phi.beta
  sigma2.beta.store = rep(NA,I); sigma2.beta.store[1] = data$sigma2.beta
  mu.store = matrix(NA,I,R); mu.store[1,] = data$mu
  beta.store = array(NA,c(I,p,R)); beta.store[1,,] = data$beta
  gamma.store = array(NA,c(I,p,R)); gamma.store[1,,] = data$gamma
  
  phi.beta.curr = phi.beta.store[1]
  sigma2.beta.curr = sigma2.beta.store[1]
  mu.curr = mu.store[1,]
  beta.curr = beta.store[1,,]
  gamma.curr = gamma.store[1,,]
  K.inv.curr = (diag(1,R) - phi.beta.curr*W)/sigma2.beta.curr
  K.curr = solve(K.inv.curr)
  probs = matrix(data$pj,p,R)
  
  for(i in 2:I)
  {
    draw = Sample.Mu(y.all=data$y.all,X.all=data$X.all,S.inv.list=data$S.inv.list,K=K.curr,
                     clusters=data$clusters,beta=beta.curr,mu.curr=mu.curr)
    mu.curr = draw
    
    draw = SSVS.Gibbs(y.all=data$y.all,X.all=data$X.all,clusters=data$clusters,S.inv.list=data$S.inv.list,K=K.curr,
                      mu=mu.curr,beta.curr=beta.curr,gamma.curr=gamma.curr,probs=probs)
    gamma.curr = draw$gamma.draw
    beta.curr = draw$beta.draw
    
    if(GP==TRUE)
    {
      draw = Sample.phi.beta(mu=mu.curr,beta=beta.curr,gamma=gamma.curr,phi.beta.curr=phi.beta.curr,step=.01,
                             K.inv.curr=K.inv.curr,D.centroids=D.centroids,sigma2.beta=sigma2.beta.curr)
      phi.beta.curr = draw$phi.beta.draw
      K.inv.curr = draw$K.inv.draw
      K.curr = solve(K.inv.curr)
      
      draw = Sample.sigma2.beta(mu=mu.curr,beta=beta.curr,gamma=gamma.curr,K.inv.curr=K.inv.curr,sigma2.beta.curr=sigma2.beta.curr)
      sigma2.beta.curr = draw$sigma2.beta.draw
      K.inv.curr = draw$K.inv.draw
      K.curr = solve(K.inv.curr)
    }else{
      K.curr = diag(1,R)
    }
    phi.beta.store[i] = phi.beta.curr
    sigma2.beta.store[i] = sigma2.beta.curr
    beta.store[i,,] = beta.curr
    mu.store[i,] = mu.curr
    gamma.store[i,,] = gamma.curr
  }
  
  burn = round(.1*I,0)
  mu.means = colMeans(mu.store[-c(1:burn),])
  beta.means = apply(beta.store[-c(1:burn),,],3,colMeans)
  gamma.means = apply(gamma.store[-c(1:burn),,],3,colMeans)
  return(gamma.means)
}

B = 100
n.seq = c(625)          # Changing sample size
p.seq = c(5,15)            # Changing number of covariates
R.seq = c(2,4,8)              # Changing number of clusters
pj.seq = c(.1,.5,.9)    # Changing proportion of nonzero covariates
sigma2.seq = c(1)       # Change noise
phi.seq = rep(0,3)      # I can't put values in here, because it depends on R
beta.seq = c(0,2)       # Changing the mean of nonzero betas
storage = matrix(NA,2*B*length(n.seq)*length(p.seq)*length(R.seq)*
                   length(pj.seq)*length(sigma2.seq)*length(phi.seq)*length(beta.seq),9)

t0 = proc.time()
count = 1
for(n in n.seq)
{
  # Create a grid for the coordinates
  coords.seq = 1:sqrt(n)
  coords = as.matrix(expand.grid(coords.seq,coords.seq))
  
  for(p in p.seq)
  {
    for(pj in pj.seq)
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
          L = eigen(W)$values
          phi.all = seq(1/min(L),1/max(L),length=1000)
          
          # Specify the phi sequence based on W
          phi.seq = rev(seq(0,quantile(phi.all,.975),length=length(phi.seq)))
          for(phi.beta in phi.seq)
          {
            sigma2.beta = 1/solve((diag(1,R) - phi.beta*W))[1]
                                 
            for(beta.mean in beta.seq)
            {
              cl = makeCluster(11)
              registerDoParallel(cl)
              runs = foreach(i=1:B,.combine='rbind',.packages=c("laGP","truncnorm","mvtnorm")) %dopar% 
              {
                data = Simulate.Data(n=n,p=p,R=R,tau2=rep(sigma2,R),phi=rep(1e-10,R),g=rep(0,R),
                                     beta.mean=beta.mean,sigma2.beta=sigma2.beta,phi.beta=phi.beta,centroids=centroids,
                                     var.select=TRUE,theta=0,pj=pj)
                
                run.GP = SSVS.run(I=5000,data=data,GP=TRUE)
                run.ind = SSVS.run(I=5000,data=data,GP=FALSE)
                AFCCF.GP = sum(data$gamma*run.GP + (1-data$gamma)*(1-run.GP))/(R*p)
                AFCCF.ind = sum(data$gamma*run.ind + (1-data$gamma)*(1-run.ind))/(R*p)
                
                return(c(AFCCF.GP,AFCCF.ind))
              }
              stopCluster(cl)
              print(proc.time()-t0)
              t0 = proc.time()
              storage[count:(count+2*B-1),1:7] = matrix(rep(c(n,p,pj,sigma2,R,round(phi.beta,2),beta.mean),2*B),ncol=7,byrow=TRUE)
              storage[count:(count+2*B-1),8] = rep(c(1,2),B) 
              storage[count:(count+2*B-1),9] = c(t(runs))
              count = count + 2*B

              print(paste("n =",n,"p =",p,"pj =",pj,"sigma2 =",sigma2,"R =",R,
                            "phi =",round(phi.beta,2),"beta =",beta.mean))
            }
          }
        }
      }
    }
  }
}

write.csv(storage, file = "~/Dropbox/Leman Research/SCOMS/CAR-SSVS/sim_results_bigN.csv")
save.image(file="~/Dropbox/Leman Research/SCOMS/CAR-SSVS/SSVS_Sim_revise_bigN.rdata")
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
# 
# 
