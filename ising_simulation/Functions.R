library(laGP)
library(truncnorm)
library(doParallel)
library(foreach)
library(mvtnorm)
library(polyclip)
library(lhs)
library(mclust)
library(LaplacesDemon)

# Calculate convex hulls for each cluster
Calc.CHs = function(coords,clusters)
{
  K = max(clusters)
  hulls = list()
  for(i in 1:K)
  {
    subgroup = coords[which(clusters==i),]
    hulls[[i]] = subgroup[chull(subgroup),]
  }
  return(hulls)
}

# Re-calculates convex hulls after a point moves from old to new cluster
Update.CHs = function(old,new,coords,clusters,CHs)
{
   old.ind = which(clusters==old)
   new.ind = which(clusters==new)
   
   subgroup = coords[old.ind,]
   CHs[[old]] = subgroup[chull(subgroup),]
   
   subgroup = coords[new.ind,]
   CHs[[new]] = subgroup[chull(subgroup),]
   
   return(CHs)
}

# Creates a vector indicating if each point is on a convex hull or not
Calc.CH.inds = function(coords,clusters,CHs)
{
  CH.inds = rep(NA,nrow(coords))
  K = max(clusters)
  for(k in 1:K)
  {
    ind = which(clusters==k)
    for(pt in ind)
    {
      CH.inds[pt] = all(coords[pt,]%in%CHs[[k]])
    }
  }
  return(CH.inds)
}

# Re-calculate vector of indicators after point moves from old to new cluster
Update.CH.inds = function(old,new,coords,clusters,CHs,CH.inds)
{
  old.ind = which(clusters==old)
  new.ind = which(clusters==new)
  
  for(pt in old.ind)
  {
    CH.inds[pt] = all(coords[pt,]%in%CHs[[old]])
  }
  for(pt in new.ind)
  {
    CH.inds[pt] = all(coords[pt,]%in%CHs[[new]])
  }
  
  return(CH.inds)
}

# For a point, find the clusters it could possibly move to
Find.Candidate.Clusters = function(coords,clusters.curr,pt)
{
  cluster.ind = rep(NA,max(clusters.curr))
  cluster.ind[clusters.curr[pt]] = 1
  
  loop = 1:max(clusters.curr)
  loop = loop[-clusters.curr[pt]]
  
  for(k in loop)
  {
    clusters.prop = clusters.curr
    clusters.prop[pt] = k
    CHs.prop = Calc.CHs(coords,clusters.prop)
    CHs.prop2 = lapply(CHs.prop,function(z)list(x=z[,1],y=z[,2]))
    
    A <- CHs.prop2[k]
    B <- CHs.prop2[-k]

    cluster.ind[k] = ifelse(length(polyclip(A,B,"intersection"))==0,1,0)
  }

  return(which(cluster.ind==1))
}

# Draw from n-dimensional Ising distribution
Sim.Ising = function(n,theta,W,B)
{
  y = rbinom(n,size=1,prob=.5)
  for(b in 1:B)
  {
    for(i in 1:n)
    {
      weights = W[i,]
      prob = (1+exp(theta*sum(weights*((y==0)-(y==1)))))^(-1)
      y[i] = rbinom(1,1,prob=prob)
    }
  }
  return(y)
}

Simulate.Data = function(n,p,R,tau2=rep(.5,R),phi=rep(1e-10,R),g=rep(1,R),
                         beta.mean=0,sigma2.beta,phi.beta,centroids,var.select=FALSE,theta=0,pj=.5)
{
  clusters = apply(distance(coords,centroids),1,which.min)
  W = 1/sqrt(distance(centroids)); diag(W) = 0
  
  # Simulate gamma j=1,...p and r=1,...,R
  
  if(var.select==TRUE)
  {
    if(theta==0)
    {
      gamma = matrix(rbinom(p*R,1,pj),p,R)
    }else{
    # Use distance between centroids to create Ising weights
    gamma = foreach(i=1:p,.combine='rbind',.export=c("Sim.Ising")) %do% Sim.Ising(n=R,theta=theta,W=W,B=100)
    }
  }else{
    gamma = matrix(1,p,R)
  }
  
  # Simulate beta j=1,...,p and r=1,...,R conditional on gamma
  beta = matrix(0,p,R)
  K = sigma2.beta*solve(diag(1,R) - phi.beta*W)
  
  for(j in 1:p)
  {
    ind = which(gamma[j,]==1)
    if(length(ind)==0)
    {
      beta[j,] = 0
    }else{
      beta[j,ind] = rmvnorm(1,matrix(beta.mean,length(ind),1),K[ind,ind,drop=F])
    }
  }
  
  # Simulate mu r=1,...,R
  mu = matrix(rmvnorm(1,rep(0,R),K),1,R)
  
  # Simulate X with correlated columns
  L = matrix(rnorm(p*5),p,5)
  X.all = scale(rmvnorm(n,rep(0,p),L%*%t(L)+diag(1,p)))
  
  # Simulate y
  y.all = matrix(NA,n,1)
  S.list = list()
  S.inv.list = list()
  D = distance(coords)
  for(r in 1:R)
  {
    ind = which(clusters==r)
    S.list[[r]] = tau2[r]*(exp(-D[ind,ind]/phi[r]) + diag(g[r],length(ind)))
    S.inv.list[[r]] = solve(S.list[[r]])
    y.all[ind,] = mu[,r] + X.all[ind,,drop=F]%*%beta[,r,drop=F] + 
                  c(rmvnorm(1,rep(0,length(ind)),S.list[[r]]))
  }
  return(list(y.all=y.all,X.all=X.all,coords=coords,centroids=centroids,clusters=clusters,sigma2.beta=sigma2.beta,phi.beta=phi.beta,
              W=W,gamma=gamma,pj=pj,tau2=tau2,phi=phi,g=g,beta=beta,mu=mu,K=K,D=D,S.list=S.list,S.inv.list=S.inv.list,theta=theta))
}

# Calculate centroids for each cluster
Calc.Centroids = function(clusters,coords)
{
  R = max(clusters)
  centroids = matrix(NA,R,2)
  
  for(r in 1:R)
  {
    ind = which(clusters==r)
    centroids[r,] = colMeans(coords[ind,])
  }
  return(centroids)
}

# Sample mu for r = 1,...,R
Sample.Mu = function(y.all,X.all,S.inv.list,K,clusters,beta,mu.curr)
{
  R = ncol(beta)
  mu.draw = c(mu.curr)
  
  for(r in 1:R)
  {
    ind = which(clusters==r)
    sig2.mu = K[r,r] - K[r,-r]%*%solve(K[-r,-r])%*%K[-r,r]
    mu.mu = K[r,-r]%*%solve(K[-r,-r])%*%matrix(mu.draw[-r])
    
    t1_S.inv = rowSums(S.inv.list[[r]])
    V = 1/(sum(t1_S.inv) + 1/sig2.mu)
    E = V*(t1_S.inv%*%(y.all[ind,,drop=F]-X.all[ind,,drop=F]%*%beta[,r,drop=F]) + mu.mu/sig2.mu)
    mu.draw[r] = rnorm(1,mean=E,sd=sqrt(V))
  }
  return(mu.draw)
}

# Sample phi for r = 1,...,R
GP.Phi.Sample = function(X.all,y.all,clusters,D,step,phi.UB,mu,beta,phi.curr,
                         tau2,g,S.list.curr,S.inv.list.curr)
{
  R = length(tau2)
  phi.draw = rep(NA,R)
  S.list.draw = list()
  S.inv.list.draw = list()

  for(r in 1:R)
  {
    ind = which(clusters==r)
    y = y.all[ind,,drop=F]
    X = X.all[ind,,drop=F]
    y.curl = y - mu[r] - X%*%beta[,r,drop=F]
    S.curr = S.list.curr[[r]]
    S.curr.inv = S.inv.list.curr[[r]]
    n = nrow(X)
    
    phi.c = phi.curr[r]
    phi.p = rtruncnorm(1,mean=phi.c,sd=step,a=0,b=phi.UB)
    S.prop = tau2[r]*(exp(-D[ind,ind]/phi.p) + diag(g[r],n))
    S.prop.inv = solve(S.prop)
    
    numer = ( -.5*determinant(S.prop,logarithm=TRUE)$modulus-.5*t(y.curl)%*%S.prop.inv%*%(y.curl)
              + Calc.phi_g.prior(D=D[ind,ind],S=S.prop,S.inv=S.prop.inv,tau2=tau2[r],g=g[r],phi=phi.p)
              - log(pnorm(phi.UB,mean=phi.p,sd=step)-pnorm(0,mean=phi.p,sd=step)) )
    denom = ( -.5*determinant(S.curr,logarithm=TRUE)$modulus-.5*t(y.curl)%*%S.curr.inv%*%(y.curl)
              + Calc.phi_g.prior(D=D[ind,ind],S=S.curr,S.inv=S.curr.inv,tau2=tau2[r],g=g[r],phi=phi.c)
              - log(pnorm(phi.UB,mean=phi.c,sd=step)-pnorm(0,mean=phi.c,sd=step))  )
    
    if(log(runif(1))<numer-denom)
    {
      phi.draw[r] = phi.p
      S.list.draw[[r]] = S.prop
      S.inv.list.draw[[r]] = S.prop.inv
    }else
    {
      phi.draw[r] = phi.c
      S.list.draw[[r]] = S.curr
      S.inv.list.draw[[r]] = S.curr.inv
    }
  }
  return(list(S.list.draw=S.list.draw,S.inv.list.draw=S.inv.list.draw,phi.draw=phi.draw))
}

# Sample g for r = 1,...,R
GP.g.Sample = function(X.all,y.all,clusters,D,step,g.UB,mu,beta,phi,
                         tau2,g.curr,S.list.curr,S.inv.list.curr)
{
  R = length(tau2)
  g.draw = rep(NA,R)
  S.list.draw = list()
  S.inv.list.draw = list()
  
  for(r in 1:R)
  {
    ind = which(clusters==r)
    y = y.all[ind,,drop=F]
    X = X.all[ind,,drop=F]
    y.curl = y - mu[r] - X%*%beta[,r,drop=F]
    S.curr = S.list.curr[[r]]
    S.curr.inv = S.inv.list.curr[[r]]
    n = nrow(X)
    
    g.c = g.curr[r]
    g.p = rtruncnorm(1,mean=g.c,sd=step,a=0,b=g.UB)
    S.prop = tau2[r]*(exp(-D[ind,ind]/phi[r]) + diag(g.p,n))
    S.prop.inv = solve(S.prop)
    
    numer = ( -.5*determinant(S.prop,logarithm=TRUE)$modulus-.5*t(y.curl)%*%S.prop.inv%*%(y.curl)
              + Calc.phi_g.prior(D=D[ind,ind],S=S.prop,S.inv=S.prop.inv,tau2=tau2[r],g=g.p,phi=phi[r])
              - log(pnorm(g.UB,mean=g.p,sd=step)-pnorm(0,mean=g.p,sd=step)) )
    denom = ( -.5*determinant(S.curr,logarithm=TRUE)$modulus-.5*t(y.curl)%*%S.curr.inv%*%(y.curl)
             + Calc.phi_g.prior(D=D[ind,ind],S=S.curr,S.inv=S.curr.inv,tau2=tau2[r],g=g.c,phi=phi[r])
              - log(pnorm(g.UB,mean=g.c,sd=step)-pnorm(0,mean=g.c,sd=step))  )
    
    if(log(runif(1))<numer-denom)
    {
      g.draw[r] = g.p
      S.list.draw[[r]] = S.prop
      S.inv.list.draw[[r]] = S.prop.inv
    }else
    {
      g.draw[r] = g.c
      S.list.draw[[r]] = S.curr
      S.inv.list.draw[[r]] = S.curr.inv
    }
  }
  return(list(S.list.draw=S.list.draw,S.inv.list.draw=S.inv.list.draw,g.draw=g.draw))
}

# Sample tau2 for r = 1,...,R
GP.tau2.Sample = function(y.all,X.all,mu,beta,clusters,S.list.curr,S.inv.list.curr,tau2.curr)
{
  R = length(tau2.curr)
  tau2.draw = rep(NA,R)
  S.list.draw = list()
  S.inv.list.draw = list()
  
  for(r in 1:R)
  {
    ind = which(clusters==r)
    y.curl = y.all[ind,] - mu[r] - X.all[ind,,drop=F]%*%beta[,r,drop=F]
    V.inv = tau2.curr[r]*S.inv.list.curr[[r]]
    b = .5*t(y.curl)%*%V.inv%*%y.curl
    a = length(ind)/2
    
    tau2.draw[r] = 1/rgamma(1,shape=a,rate=b)
    S.list.draw[[r]] = (S.list.curr[[r]]/tau2.curr[r])*tau2.draw[r]
    S.inv.list.draw[[r]] = (S.inv.list.curr[[r]]*tau2.curr[r])/tau2.draw[r]
  }
  return(list(S.list.draw=S.list.draw,S.inv.list.draw=S.inv.list.draw,tau2.draw=tau2.draw))
}

# Computes the .5*log(det(...)) using the reference prior for phi,g for a single region r
Calc.phi_g.prior = function(D,S,S.inv,tau2,g,phi)
{
  G.inv = S.inv*tau2
  K = S/tau2 - diag(g,nrow(S))
  dK = K*D/phi^2
  G.inv_dK = G.inv%*%dK
  
  Sig = matrix(NA,3,3)
  
  Sig[1,1] = sum(diag(G.inv_dK))^2
  Sig[1,2] <- Sig[2,1] <- sum(diag(G.inv%*%G.inv_dK))
  Sig[1,3] <- Sig[3,1] <- sqrt(Sig[1,1])
  Sig[2,2] = sum(diag(G.inv%*%G.inv))
  Sig[2,3] <- Sig[3,2] <- sum(diag(G.inv))
  Sig[3,3] = length(nrow(K))
  
  prior = .5*determinant(Sig,logarithm=TRUE)$modulus
  return(prior)
}

# Sample gamma, beta for r = 1,...,R and j = 1,...,p

SSVS.Gibbs = function(y.all,X.all,clusters,S.inv.list,K,mu,beta.curr,gamma.curr,probs)
{
  p = nrow(beta.curr)
  R = ncol(beta.curr)
  beta.draw = beta.curr
  gamma.draw = gamma.curr
  
  for(r in 1:R)
  {
    ind = which(clusters==r)
    y = matrix(y.all[ind,],ncol=1)
    X = X.all[ind,]
    n = nrow(X)
    
    for(j in 1:p)
    {
      ind = which(gamma.draw[j,]==1)
      ind = ind[ind!=r]
      
      if(length(ind)==0)
      {
        sigma2.beta = K[r,r]
        mu.beta = 0
      }else{
        sigma2.beta = K[r,r] - K[r,ind]%*%solve(K[ind,ind])%*%K[ind,r]
        mu.beta = K[r,ind]%*%solve(K[ind,ind])%*%matrix(beta.draw[j,ind],ncol=1)
      }
      
      p.j = probs[j,r]
      x.j = matrix(X[,j],n,1)
      X.curl = matrix(X[,-j],n,p-1)
      y.curl = y - mu[r] - X.curl%*%beta.draw[-j,r,drop=F]
      tx.j_S.inv = t(x.j)%*%S.inv.list[[r]]
      
      V = solve(tx.j_S.inv%*%x.j + 1/sigma2.beta)
      E = V%*%(tx.j_S.inv%*%y.curl + mu.beta/sigma2.beta)
      
      log.BF.inv = dnorm(0,E,sqrt(V),log=TRUE) - dnorm(0,mu.beta,sqrt(sigma2.beta),log=TRUE)
      BF.inv = min(exp(log.BF.inv),1e10) 
      
      prob = ( 1 + ((1-p.j)/p.j)*BF.inv )^(-1)
      gamma.draw[j,r] = rbinom(1,1,prob=prob)
      beta.draw[j,r] = gamma.draw[j,r]*rnorm(1,E,sqrt(V)) + 0 
    }
  }
  
  return(list(gamma.draw=gamma.draw,beta.draw=beta.draw))
}

Log.L = function(y.all,X.all,beta.mat,mu.vec,Sigma.list,clusters)
{
  R = ncol(beta.mat)
  LL = rep(NA,R)
  for(r in 1:R)
  {
    ind = which(clusters==r)
    LL[r] = dmvnorm(y.all[ind,],mean=mu.vec[r]+X.all[ind,]%*%beta.mat[,r],sigma=Sigma.list[[r]],log=TRUE)
  }
  return(sum(LL))
}

Update.Weights = function(coords,clusters)
{
  R = max(clusters)
  centroids = matrix(NA,R,2)
  for(r in 1:R)
  {
    ind = which(clusters==r)
    centroids[r,] = colMeans(coords[ind,])
  }
  weight.mat = exp(-distance(centroids))
  diag(weight.mat) = 0
  weight.mat = sweep(weight.mat,2,colSums(weight.mat),FUN="/"); weight.mat = .5*(weight.mat + t(weight.mat)); 
  weight.mat = sweep(weight.mat,2,colSums(weight.mat),FUN="/")
  return(weight.mat)
}

Update.p = function(gamma.mat.curr,weight.mat,theta)
{
  R = ncol(gamma.mat.curr)
  p = nrow(gamma.mat.curr)
  p.mat = matrix(NA,p,R)
  for(r in 1:R)
  {
    weights = weight.mat[-r,r]
    for(j in 1:p)
    {
      gamma.other = gamma.mat.curr[j,-r]
      p.mat[j,r] = (1+exp(theta*sum(weights*((gamma.other==0)-(gamma.other==1)))))^(-1)
    }
  }
  return(p.mat)
}

Cluster.Assign.CH = function(X.all,y.all,clusters.curr,beta,mu,gamma,S.list,S.inv.list,
                             tau2,g,phi,D,coords,CHs.curr,CH.inds.curr)
{
  n = nrow(X.all)
  R = length(tau2)
  clusters.new = clusters.curr
  CH.inds.draw = CH.inds.curr
  CHs.draw = CHs.curr

  for(i in 1:n)
  {
    draw = clusters.new[i]
    
    # Check to see if pt is on CH and if the cluster won't get too small if moved
    if((CH.inds.draw[i]==1)&&(sum(clusters.new==clusters.new[i])>5))
    {
      my.options = Find.Candidate.Clusters(coords,clusters.new,i)
      if(length(my.options)>1)
      {
        # If a cluster isn't an option, it gets a Pr = 0
        log.weights = rep(NA,R)
        log.weights[which(1:R %in% my.options == FALSE)] = -Inf
        
        # By storing S12 and S22 as I go, I can quickly update later
        S12.list = list() 
        S12_S22.inv.list = list()
        S22.list = list()
        S22.inv.list = list()
        
        # Find log probabilities for each option
        for(j in my.options)
        {
          # Calculation is easier for cluster it is already in
          if(clusters.new[i]==j)
          {
            S = S.list[[j]] # To calculate weight for region r, I need its covariance
            S.inv = S.inv.list[[j]]
            ind = which(clusters.new==j)
            ind.other = ind[ind!=i]
            S.ind = which(ind==i)
            
            V = 1/S.inv[S.ind,S.ind]
            E = ( mu[j] + X.all[i,,drop=F]%*%beta[,j,drop=F] + 
                (-V)*S.inv[S.ind,-S.ind]%*%(y.all[ind.other,]-mu[j]-X.all[ind.other,,drop=F]%*%beta[,j,drop=F]) )
            log.weights[j] = dnorm(y.all[i,],mean=E,sd=sqrt(V),log=TRUE)
          }else{
            # Calculation for other clusters needs rank-one updates
            ind = which(clusters.new==j)
            D.new = D[i,ind]; D.new = D.new[D.new!=0]
            S12 = tau2[j]*exp(-D.new/phi[j])
            S22 = S.list[[j]]
            S22.inv = S.inv.list[[j]]
            
            S12_S22.inv = S12%*%S22.inv
            V = max( (tau2[j]*(1+g[j]) - S12_S22.inv%*%matrix(S12,ncol=1)), 1e-8 )
            E = mu[j] + X.all[i,]%*%beta[,j,drop=F] + S12_S22.inv%*%(y.all[ind,]-mu[j]-X.all[ind,]%*%beta[,j,drop=F])
            log.weights[j] = dnorm(y.all[i,],mean=E,sd=sqrt(V),log=TRUE)           
            
            S12.list[[j]] = S12; S12_S22.inv.list[[j]] = S12_S22.inv
            S22.list[[j]] = S22; S22.inv.list[[j]] = S22.inv
          }
        }
        
        # Based on probabilities, draw a cluster assignment
        m = which.max(log.weights)
        log.Denom = log.weights[m] + log(1+sum(exp(log.weights[-m]-log.weights[m])))
        my.weights = exp(log.weights-log.Denom)
        draw = sample(1:R,size=1,prob=my.weights)
      }
    }
    
    if(draw!=clusters.new[i])
    {
      # Add point to its new covariance matrix
      ind = which(clusters.new==draw)
      S.temp = rbind(0,cbind(0,S.list[[draw]]))
      S.temp[1,1] = tau2[draw]*(1+g[draw])
      S.temp[1,-1] <- S.temp[-1,1] <- S12.list[[draw]]
      o = order(c(i,ind))
      S.temp = S.temp[o,o]
      S.list[[draw]] = S.temp[o,o]
      
      S.inv.temp = matrix(NA,length(ind)+1,length(ind)+1)
      S12_S22.inv = S12_S22.inv.list[[draw]]
      S12 = matrix(S12.list[[draw]],nrow=1); S22 = S22.list[[draw]]
      S22.inv = S22.inv.list[[draw]]
      
      S11 = (tau2[draw]*(1+g[draw]))
      S.inv.temp[1,1] = 1/(S11 - S12_S22.inv%*%matrix(S12,ncol=1))
      S.inv.temp[1,-1] <- S.inv.temp[-1,1] <- (-(S.inv.temp[1,1]))*S12_S22.inv
      S.inv.temp[-1,-1] = 
        ( S22.inv + t(S12_S22.inv)%*%(S12_S22.inv)*S.inv.temp[1,1])
      S.inv.temp = S.inv.temp[o,o]
      S.inv.list[[draw]] = S.inv.temp
      
      # Remove point from its old covariance matrix
      ind = which(clusters.new==clusters.new[i])
      rm.ind = which(ind==i)
      S.temp = S.list[[clusters.new[i]]]
      S.list[[clusters.new[i]]] = S.temp[-rm.ind,-rm.ind]
      
      # Remove point from its old precision matrix
      S.inv.curr = S.inv.list[[clusters.new[i]]]
      S12_S22.inv = matrix(-(1/S.inv.curr[rm.ind,rm.ind])*S.inv.curr[rm.ind,-rm.ind],nrow=1)
      S.inv.temp = S.inv.curr[-rm.ind,-rm.ind] - t(S12_S22.inv)%*%S12_S22.inv*S.inv.curr[rm.ind,rm.ind]
      S.inv.list[[clusters.new[i]]] = S.inv.temp
      
      clusters.new[i] = draw
      
      CHs.draw = Update.CHs(old=clusters.curr[i],new=clusters.new[i],coords=coords,clusters=clusters.new,CHs=CHs.draw)
      CH.inds.draw = Update.CH.inds(old=clusters.curr[i],new=clusters.new[i],coords=coords,
                                    clusters=clusters.new,CHs=CHs.draw,CH.inds=CH.inds.draw)
    }
  }
  return(list(S.list.draw=S.list,S.inv.list.draw=S.inv.list,clusters.draw=clusters.new,CHs.draw=CHs.draw,
              CH.inds.draw=CH.inds.draw))
}

Sample.phi.beta = function(mu,beta,gamma,phi.beta.curr,step,K.inv.curr,D.centroids,sigma2.beta)
{
  R = ncol(beta)
  p = nrow(beta)
  mu = matrix(mu,nrow=1)
  
  W = 1/D.centroids; diag(W) = 0
  L = eigen(W)$values
  L1 = min(L); Ln = max(L)
  phi.LB = 1/L1; phi.UB = 1/Ln
  
  if(sum(gamma)>0)
  {
    beta.sub = beta[which(rowSums(gamma)>0),,drop=F]
    gamma.sub = gamma[which(rowSums(gamma)>0),,drop=F]
    beta.j = lapply(split(beta.sub,seq(nrow(beta.sub))),function(x) x[x!=0])
    
    phi.beta.prop = rtruncnorm(1,a=phi.LB,b=phi.UB,mean=phi.beta.curr,sd=step)
    K.inv.prop = (diag(1,R)-phi.beta.prop*W)/sigma2.beta
    
    K.j.inv.prop = lapply(split(gamma.sub,seq(nrow(gamma.sub))),function(x) K.inv.prop[which(x==1),which(x==1)])
    K.j.inv.curr = lapply(split(gamma.sub,seq(nrow(gamma.sub))),function(x) K.inv.curr[which(x==1),which(x==1)])
    
    log.det.prop = determinant(K.inv.prop,logarithm=TRUE)$modulus + sum(sapply(K.j.inv.prop,function(x) determinant(as.matrix(x),logarithm=TRUE)$modulus))
    log.det.curr = determinant(K.inv.curr,logarithm=TRUE)$modulus + sum(sapply(K.j.inv.curr,function(x) determinant(as.matrix(x),logarithm=TRUE)$modulus))
    
    SS.prop = ( mu%*%K.inv.prop%*%t(mu) +
                  sum(unlist(Map("%*%",Map("%*%",beta.j,K.j.inv.prop),beta.j))) )
    SS.curr = ( mu%*%K.inv.curr%*%t(mu) +
                  sum(unlist(Map("%*%",Map("%*%",beta.j,K.j.inv.curr),beta.j))) )
  }else{
    phi.beta.prop = rtruncnorm(1,a=phi.LB,b=phi.UB,mean=phi.beta.curr,sd=step)
    K.inv.prop = (diag(1,R)-phi.beta.prop*W)/sigma2.beta

    log.det.prop = determinant(K.inv.prop,logarithm=TRUE)$modulus
    log.det.curr = determinant(K.inv.curr,logarithm=TRUE)$modulus 
    
    SS.prop =  mu%*%K.inv.prop%*%t(mu) 
    SS.curr =  mu%*%K.inv.curr%*%t(mu) 
  }
  
  numer = ( .5*log.det.prop - .5*SS.prop
            + log(sqrt( sum((L/(1-phi.beta.prop*L))^2) - (1/n)*sum(L/(1-phi.beta.prop*L))^2 ))
            - log(pnorm(phi.UB,mean=phi.beta.prop,sd=step)-pnorm(phi.LB,mean=phi.beta.prop,sd=step)) )
  denom = ( .5*log.det.curr - .5*SS.curr
            + log(sqrt( sum((L/(1-phi.beta.curr*L))^2) - (1/n)*sum(L/(1-phi.beta.curr*L))^2 ))
            - log(pnorm(phi.UB,mean=phi.beta.curr,sd=step)-pnorm(phi.LB,mean=phi.beta.curr,sd=step))  )
  
  if(log(runif(1))<numer-denom)
  {
    phi.beta.draw = phi.beta.prop
    K.inv.draw = K.inv.prop
  }else
  {
    phi.beta.draw = phi.beta.curr
    K.inv.draw = K.inv.curr
  }
  return(list(phi.beta.draw = phi.beta.draw, K.inv.draw = K.inv.draw))
}

Sample.sigma2.beta = function(mu,beta,gamma,K.inv.curr,sigma2.beta.curr)
{
  R = ncol(beta)
  mu = matrix(mu,nrow=1)

  if(sum(gamma)>0)
  {
    beta.sub = beta[which(rowSums(gamma)>0),,drop=F]
    gamma.sub = gamma[which(rowSums(gamma)>0),,drop=F]
    beta.j = lapply(split(beta.sub,seq(nrow(beta.sub))),function(x) x[x!=0])
    K.inv.curl = K.inv.curr*sigma2.beta.curr
    K.j.inv.curl = lapply(split(gamma.sub,seq(nrow(gamma.sub))),function(x) K.inv.curl[which(x==1),which(x==1)])
    
    SS = ( mu%*%K.inv.curl%*%t(mu) +
             sum(unlist(Map("%*%",Map("%*%",beta.j,K.j.inv.curl),beta.j))) )
  }else{
    K.inv.curl = K.inv.curr*sigma2.beta.curr
    SS = mu%*%K.inv.curl%*%t(mu) 
  }

  a.star = .5*(R + sum(gamma))
  b.star = .5*SS
  sigma2.beta.draw = 1/rgamma(1,shape=a.star,rate=b.star)
  K.inv.draw = K.inv.curl/sigma2.beta.draw

  return(list(sigma2.beta.draw=sigma2.beta.draw,K.inv.draw=K.inv.draw))
}

Sample.theta = function(theta.curr,gamma,W,a,b,theta.UB,step.size)
{
  theta.prop = rtruncnorm(1,a=0,b=theta.UB,mean=theta.curr,sd=step.size)
  p.prop = t( (1 + exp(theta.prop*W%*%t((gamma==0)-(gamma==1))))^(-1) )
  p.curr = t( (1 + exp(theta.curr*W%*%t((gamma==0)-(gamma==1))))^(-1) )
  
  log.numer = ( sum(log(p.prop[gamma==1]))+sum(log((1-p.prop)[gamma==0])) - (a-1)*log(theta.prop) - b*theta.prop
              - log(pnorm(theta.UB,theta.prop,step.size)-pnorm(0,theta.prop,step.size)) )
  log.denom = ( sum(log(p.curr[gamma==1]))+sum(log((1-p.curr)[gamma==0])) - (a-1)*log(theta.curr) - b*theta.curr
              - log(pnorm(theta.UB,theta.curr,step.size)-pnorm(0,theta.curr,step.size)) )
  log.ratio = log.numer - log.denom
  
  if(log(runif(1))<log.ratio)
  {
    theta.draw = theta.prop
  }else{
    theta.draw = theta.curr
  }
  return(theta.draw)
}

Cluster.Assign.Probit = function(X.all,y.all,clusters.curr,beta.mat,mu.vec,gamma.mat,Sigma.list,Sigma.inv.list,
                                 tau2.vec,g.vec,phi.vec,D,z.all)
{
  n = nrow(X.all)
  R = length(tau2.vec)
  clusters.new = clusters.curr
  y.new = y.all
  z = c(z.all)
  LBs = rep(0,n); LBs[z==0] = -Inf
  UBs = rep(0,n); UBs[z==1] = Inf
  
  for(i in 1:n)
  {
    log.weights = rep(NA,R)
    S12.list = list() # By storing S12 and S22 as I go, I can quickly update later
    S12_S22.inv.list = list()
    S22.list = list()
    S22.inv.list = list()
    
    for(j in 1:R)
    {
      if(clusters.new[i]==j)
      {
        S = Sigma.list[[j]] # To calculate weight for region r, I need its covariance
        S.inv = Sigma.inv.list[[j]]
        ind = which(clusters.new==j)
        ind.other = ind[ind!=i]
        S.ind = which(ind==i)
        
        V = 1/S.inv[S.ind,S.ind]
        E = mu.vec[j] + X.all[i,]%*%beta.mat[,j] + (-V)*S.inv[S.ind,-S.ind]%*%(y.new[ind.other,]-mu.vec[j]-X.all[ind.other,]%*%beta.mat[,j])
        log.weights[j] = dnorm(y.new[i,],mean=E,sd=sqrt(V),log=TRUE)
        y.temp = rtruncnorm(1,a=LBs[i],b=UBs[i],mean=E,sd=sqrt(V))
      }else{
        ind = which(clusters.new==j)
        D.new = D[i,ind]; D.new = D.new[D.new!=0]
        S12 = tau2.vec[j]*exp(-D.new/phi.vec[j] + diag(1e-8,length(ind)))
        S22 = Sigma.list[[j]]
        S22.inv = Sigma.inv.list[[j]]
        
        S12_S22.inv = S12%*%S22.inv
        V = (tau2.vec[j]*(1+g.vec[j]) - S12_S22.inv%*%matrix(S12,ncol=1))
        V = max(V,1e-6)
        E = mu.vec[j] + X.all[i,]%*%beta.mat[,j] + S12_S22.inv%*%(y.new[ind,]-mu.vec[j]-X.all[ind,]%*%beta.mat[,j])
        log.weights[j] = dnorm(y.new[i,],mean=E,sd=sqrt(V),log=TRUE)
        
        S12.list[[j]] = S12; S12_S22.inv.list[[j]] = S12_S22.inv
        S22.list[[j]] = S22; S22.inv.list[[j]] = S22.inv
      }
    }
    
    y.new[i,] = y.temp
    m = which.max(log.weights)
    log.Denom = log.weights[m] + log(1+sum(exp(log.weights[-m]-log.weights[m])))
    my.weights = exp(log.weights-log.Denom)
    
    draw = sample(1:R,size=1,prob=my.weights)
    draw = clusters.new[i] # THIS IS A HACK CHANGE IT BACK
    if(draw!=clusters.new[i])
    {
      # Add point to its new covariance matrix
      ind = which(clusters.new==draw)
      Sigma.temp = rbind(0,cbind(0,Sigma.list[[draw]]))
      Sigma.temp[1,1] = tau2.vec[draw]*(1+g.vec[draw])
      Sigma.temp[1,-1] <- Sigma.temp[-1,1] <- S12.list[[draw]]
      o = order(c(i,ind))
      Sigma.temp = Sigma.temp[o,o]
      Sigma.list[[draw]] = Sigma.temp[o,o]
      
      Sigma.inv.temp = matrix(NA,length(ind)+1,length(ind)+1)
      S12_S22.inv = S12_S22.inv.list[[draw]]
      S12 = matrix(S12.list[[draw]],nrow=1); S22 = S22.list[[draw]]
      S22.inv = S22.inv.list[[draw]]
      
      S11 = (tau2.vec[draw]*(1+g.vec[draw]))
      Sigma.inv.temp[1,1] = 1/(S11 - S12_S22.inv%*%matrix(S12,ncol=1))
      Sigma.inv.temp[1,-1] <- Sigma.inv.temp[-1,1] <- (-(Sigma.inv.temp[1,1]))*S12_S22.inv
      Sigma.inv.temp[-1,-1] =
        ( S22.inv + t(S12_S22.inv)%*%(S12_S22.inv)*Sigma.inv.temp[1,1])
      Sigma.inv.temp = Sigma.inv.temp[o,o]
      Sigma.inv.list[[draw]] = Sigma.inv.temp
      
      # Remove point from its old covariance matrix
      ind = which(clusters.new==clusters.new[i])
      rm.ind = which(ind==i)
      Sigma.temp = Sigma.list[[clusters.new[i]]]
      Sigma.list[[clusters.new[i]]] = Sigma.temp[-rm.ind,-rm.ind]
      
      # Remove point from its old precision matrix
      Sigma.inv.curr = Sigma.inv.list[[clusters.new[i]]]
      S12_S22.inv = matrix(-(1/Sigma.inv.curr[rm.ind,rm.ind])*Sigma.inv.curr[rm.ind,-rm.ind],nrow=1)
      Sigma.inv.temp = Sigma.inv.curr[-rm.ind,-rm.ind] - t(S12_S22.inv)%*%S12_S22.inv*Sigma.inv.curr[rm.ind,rm.ind]
      Sigma.inv.list[[clusters.new[i]]] = Sigma.inv.temp
      
      clusters.new[i] = draw
    }
  }
  return(list(y.new=y.new,Sigma.list.draw=Sigma.list,Sigma.inv.list.draw=Sigma.inv.list,clusters.draw=clusters.new))
}

Latent.Draw = function(X.all,y.curr,z.all,clusters,beta.mat,mu.vec,Sigma.inv.list)
{
  y.draw = y.curr
  n = nrow(X.all)
  LBs <- UBs <- rep(0,n)
  UBs[z.all==1] = Inf
  LBs[z.all==0] = -Inf
  
  for(i in 1:n)
  {
    r = clusters[i]
    ind = which(clusters==r) # which elements of 1,...,n are in the cluster
    ind.other = ind[ind!=i] # which elements are in the cluster, not counting itself
    S.ind = which(ind==i) # which member of cluster i is
    tau2.y = 1/Sigma.inv.list[[r]][S.ind,S.ind]
    mu.y = mu.vec[r] + X.all[i,]%*%beta.mat[,r] + (-tau2.y)*Sigma.inv.list[[r]][S.ind,-S.ind]%*%(y.draw[ind.other,] - mu.vec[r] - X.all[ind.other,]%*%beta.mat[,r])
    y.draw[i,] = rtruncnorm(1,mean=mu.y,sd=sqrt(tau2.y),a=LBs[i],b=UBs[i])
  }
  return(y.draw=y.draw)
}

Find.Candidate.Clusters.newpt = function(coords,coords.new,clusters) 
{
  cluster.ind = rep(NA,max(clusters))
  
  for(k in 1:max(clusters))
  {
    clusters.prop = c(k,clusters)
    coords.prop = rbind(coords.new,coords)
    hulls.prop = Calc.Convex.Hulls(grid=coords.prop,groups=clusters.prop)
    
    base = Polygons(list(Polygon(as.matrix(hulls.prop[[k]]))),"base")
    sel = Polygons(lapply(hulls.prop[-k],function(c)Polygon(c)),"sel")
    shape = SpatialPolygons(list(base, sel))
    ind = if(is.null(gIntersection(shape["base"], shape["sel"]))){FALSE}else{TRUE}
    
    cluster.ind[k] = 1-ind
  }
  
  return(which(cluster.ind==1))
}

Cross.Validate.CH = function(X.all,y.all,clusters,beta.mat,mu.vec,Sigma.list,Sigma.inv.list,
                             tau2.vec,g.vec,phi.vec,D,coords,X.new,coords.new)
{
  n = nrow(X.all)
  R = length(tau2.vec)
  clusters.assign = rep(NA,nrow(X.new))
  y.pred = matrix(NA,nrow=nrow(X.new))
  
  for(i in 1:nrow(X.new))
  {
    my.options = Find.Candidate.Clusters.newpt(coords=coords,coords.new=coords.new[i,],clusters=clusters)
    log.weights = rep(NA,R)
    log.weights[which(1:R %in% my.options == FALSE)] = -Inf
    
    if(length(my.options)==0)
    {
      draw = clusters[which.min(distance(matrix(coords.new[i,],nrow=1),coords))]
      
    }else{
      if(length(my.options)==1){draw = my.options[1]}else{
        
        for(j in my.options)
        {
          ind = which(clusters==j)
          D.new = distance(matrix(coords.new[i,],nrow=1),coords[ind,])
          S12 = tau2.vec[j]*exp(-D.new/phi.vec[j])
          S22 = Sigma.list[[j]]
          S22.inv = Sigma.inv.list[[j]]
          
          S12_S22.inv = S12%*%S22.inv
          V = max( (tau2.vec[j]*(1+g.vec[j]) - S12_S22.inv%*%matrix(S12,ncol=1)), 1e-8 )
          E = mu.vec[j] + X.all[i,]%*%beta.mat[,j] + S12_S22.inv%*%(y.all[ind,]-mu.vec[j]-X.all[ind,]%*%beta.mat[,j])
          log.weights[j] = dnorm(y.all[i,],mean=E,sd=sqrt(V),log=TRUE)
        }
        m = which.max(log.weights)
        log.Denom = log.weights[m] + log(1+sum(exp(log.weights[-m]-log.weights[m])))
        my.weights = exp(log.weights-log.Denom)
        
        draw = sample(1:R,size=1,prob=my.weights)
      }
    }
    y.pred[i,] = rnorm(1,mean=mu.vec[draw] + X.new[i,]%*%beta.mat[,draw],sd=1)
  }
  
  return(y.pred)
}
