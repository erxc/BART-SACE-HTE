
rm(list = ls())

gendata = function(){
  ## study parameters
  N = 3000 # sample size
  # N = 1000
  
  sigma2 = 2 # outcome variance
  
  ## cluster ids and two covariates
  x1 = rnorm(N)
  x2 = rnorm(N)
  x3 = rnorm(N)
  x4 = rnorm(N)
  x5 = rnorm(N)
  x6 = rnorm(N)
  X = cbind(x1, x2, x3, x4, x5, x6)
  
  ## strata model parameters
  # generate latent variable for G group status
  p00 = pnorm(-1+0.8*sin(pi*x1*x2)-0.6*(x3-0.5)**2+0.2*exp(-abs(x4))+0.8*tanh(x5*x6))
  G00 = rbinom(N,1,p00)
  
  G10 = rep(0,N)
  p10 = pnorm(-0.5+0.7*cos(pi*x1[G00==0]*x2[G00==0])+0.5*(x3[G00==0]-0.5)**2-
                0.7*exp(-abs(x4[G00==0]))+0.5*tanh(x5[G00==0]*x6[G00==0]))
  G10[G00==0] = rbinom(sum(G00==0),1,p10)
  G11 = rep(1,N) - G10 - G00
  
  # outcome model parameters
  d = rbinom(N,1,1/2) # individual treatment
  
  resid = rnorm(N, 0, sqrt(sigma2)) # outcome variance
  
  Ytmp = rep(NA, N) # outcome (NA are death)
  
  X111 = X[(G11==1)&(d==1),]
  X101 = X[(G10==1)&(d==1),]
  X110 = X[(G11==1)&(d==0),]
  
  Ytmp[(G11==1)&(d==1)] = 0.8 + 0.7*cos(pi*X111[,1]*X111[,2]) + 
    0.4*(X111[,3]-0.5)**3 + 0.7*tanh(X111[,4]*X111[,5]*X111[,6])
  
  Ytmp[(G10==1)&(d==1)] = 0.3 + 0.2*sin(pi*X101[,1]*X101[,2]) +
    0.2*(X101[,3]-0.5)**3 + 0.2*tanh(X101[,4]*X101[,5]*X101[,6])
  
  Ytmp[(G11==1)&(d==0)] = -0.4 + 0.6*sin(pi*X110[,1]*X110[,2]) +
    0.1*(X110[,3]-0.5)**3 + 0.5*tanh(X110[,4]*X110[,5]*X110[,6])
  
  Y = Ytmp + resid
  
  G = 1*G10 + 2*G11
  
  # true TE (SACE)
  XG2 = X[G==2,]
  
  TE = mean(0.8 + 0.7*cos(pi*XG2[,1]*XG2[,2]) + 
              0.4*(XG2[,3]-0.5)**3 + 0.7*tanh(XG2[,4]*XG2[,5]*XG2[,6]) - 
              (-0.4 + 0.6*sin(pi*XG2[,1]*XG2[,2]) +
                 0.1*(XG2[,3]-0.5)**3 + 0.5*tanh(XG2[,4]*XG2[,5]*XG2[,6])))
  
  # For G, the true status of 0, 1 or 2 will not be known, 
  # and death in treatment group and survival in control group have known G status   
  return(data.frame(cbind(Y, x1, x2, x3, x4, x5, x6, d, id=1:N, G, TE)))
}

## Combined Gibbs sampler
combined.Gibbs = function(df=df, S){ # S is the number of iterations
  
  # model information
  N = nrow(df)
  Nst = length(which(!is.na(df$Y)))
  X = as.matrix(df[,c('x1','x2','x3','x4','x5','x6')])
  Y = df$Y
  D = df$d
  J = 200
  w = 2
  
  ## initial values
  # outcome model
  library(BayesTree, quietly=T)
  m.oc.11.1 = rep(NA,N)
  m.oc.fit = bart(x.train=X[!is.na(Y),], y.train=Y[!is.na(Y)], x.test=X[is.na(Y),], ndpost=200, verbose=F)
  m.oc.11.1[!is.na(Y)] = m.oc.fit$yhat.train.mean
  m.oc.11.1[is.na(Y)] = m.oc.fit$yhat.test.mean
  m.oc.11.0 = m.oc.11.1
  m.oc.10.1 = m.oc.11.1
  sigma2 = m.oc.fit$sigest
  G = df$G # from the truth
  G0 = df$G # the true strata membership
  
  # strata model
  library(geepack, quietly=T)
  lmm1 = geeglm((G==0) ~ x1+x2+x3+x4+x5+x6-1, family=binomial, id=id, data=df)
  m.Z = lmm1$fitted.values[,1]
  
  lmm2 = geeglm((G==1) ~ x1+x2+x3+x4+x5+x6-1, family=binomial, id=id, data=df[G!=0,])
  m.W = (X %*% lmm2$coefficients)[,1]
  
  ## priors 
  # outcome model
  cc = 0.001  # prior variance terms
  dd = 0.001
  
  library(truncnorm, quietly=T)
  # every individual has a Z
  Z = ifelse (G<1, rtruncnorm(1,a=0,mean=0,sd=sqrt(1)), rtruncnorm(1,b=0,mean=0,sd=sqrt(1)))
  # every individual in strata 10 & 11 has a W
  W = rep(NA,N)
  W[G==1] = rtruncnorm(sum(G==1), a=0, mean=0, sd=sqrt(1))
  W[G==2] = rtruncnorm(sum(G==2), b=0, mean=0, sd=sqrt(1))
  W[G==0] = NA
  
  ## store chain S
  SIGMA_111 = rep(0,S)
  SIGMA_110 = rep(0,S)
  SIGMA_101 = rep(0,S)
  TE = rep(0,S) # sample TE by estimating the potential outcomes of G=11.
  PEHE = array(NA, c(length(which(G0==2)), S))
  IE_RB_num = rep(0,S) # relative bias of ISCE, numerator
  IE_RB_den = rep(0,S) # relative bias of ISCE, denominator
  IE_RM = array(NA, c(length(which(G0==2)), S))
  
  M.oc.11.1 = matrix(NA, nrow=N, ncol=S) # BART means for G=11 & D=1
  M.oc.11.0 = matrix(NA, nrow=N, ncol=S) # BART means for G=11 & D=0
  M.oc.10.1 = matrix(NA, nrow=N, ncol=S) # BART means for G=10 & D=1
  
  GF = matrix(NA, nrow=3, ncol=S)
  GV = matrix(NA, nrow=N, ncol=S)
  M.Z = matrix(NA, nrow=N, ncol=S) # BART means for Z
  M.W = matrix(NA, nrow=N, ncol=S) # BART means for W
  
  library(dbarts, quietly=T)
  # 1000: ntrees = 50/50
  # 3000: ntrees = 75/75
  db.control.1 = dbartsControl(updateState=F, verbose=F, n.burn=0L, n.samples=1L, 
                               n.trees=50L, n.thin=1L, n.chains=1L)
  db.control.2 = dbartsControl(updateState=T, verbose=F, n.burn=0L, n.samples=1L, 
                               n.trees=50L, n.thin=1L, n.chains=1L)
  
  for(i in 1:S){
  
    # update outcome BART means, m.oc.11.1, m.oc.11.0, m.oc.10.1
    # use df[(G==g)&(D==d),] to train, and df[!((G==g)&(D==d)),] to test and predict
    # have to save all BART means, because G may change after each iteration
    if (sum((G==2)&(D==1)) > ncol(X)+1) {
      m.oc.11.1.train = df[(G==2)&(D==1), c('Y','x1','x2','x3','x4','x5','x6')]
      m.oc.11.1.test = df[!((G==2)&(D==1)), c('x1','x2','x3','x4','x5','x6')]
      m.oc.11.1.sampler = dbarts(Y ~ x1+x2+x3+x4+x5+x6, m.oc.11.1.train, m.oc.11.1.test, 
                                 node.prior = normal(2), control = db.control.1)
      m.oc.11.1.samples = m.oc.11.1.sampler$run()
      m.oc.11.1 = rep(NA,N)
      m.oc.11.1[(G==2)&(D==1)] = m.oc.11.1.samples$train[,1]
      m.oc.11.1[!((G==2)&(D==1))] = m.oc.11.1.samples$test[,1]
      
      rate_111 = dd + 0.5*(sum((Y[(G==2)&(D==1)]-m.oc.11.1[(G==2)&(D==1)])**2))
      sigma2_111 = rgamma(1, shape=(cc+sum((G==2)&(D==1))/2), rate=rate_111)**(-1)
    }
    M.oc.11.1[,i] = m.oc.11.1
    SIGMA_111[i] = sigma2_111
    
    m.oc.11.0.train = df[(G==2)&(D==0), c('Y','x1','x2','x3','x4','x5','x6')]
    m.oc.11.0.test = df[!((G==2)&(D==0)), c('x1','x2','x3','x4','x5','x6')]
    m.oc.11.0.sampler = dbarts(Y ~ x1+x2+x3+x4+x5+x6, m.oc.11.0.train, m.oc.11.0.test, 
                               node.prior = normal(2), control = db.control.1)
    m.oc.11.0.samples = m.oc.11.0.sampler$run()
    m.oc.11.0 = rep(NA,N)
    m.oc.11.0[(G==2)&(D==0)] = m.oc.11.0.samples$train[,1]
    m.oc.11.0[!((G==2)&(D==0))] = m.oc.11.0.samples$test[,1]
    
    rate_110 = dd + 0.5*(sum((Y[(G==2)&(D==0)]-m.oc.11.0[(G==2)&(D==0)])**2))
    sigma2_110 = rgamma(1, shape=(cc+sum((G==2)&(D==0))/2), rate=rate_110)**(-1)
    M.oc.11.0[,i] = m.oc.11.0
    SIGMA_110[i] = sigma2_110
    
    if (sum((G==1)&(D==1)) > ncol(X)+1) {
      m.oc.10.1.train = df[(G==1)&(D==1), c('Y','x1','x2','x3','x4','x5','x6')]
      m.oc.10.1.test = df[!((G==1)&(D==1)), c('x1','x2','x3','x4','x5','x6')]
      m.oc.10.1.sampler = dbarts(Y ~ x1+x2+x3+x4+x5+x6, m.oc.10.1.train, m.oc.10.1.test, 
                                 node.prior = normal(2), control = db.control.1)
      m.oc.10.1.samples = m.oc.10.1.sampler$run()
      m.oc.10.1 = rep(NA,N)
      m.oc.10.1[(G==1)&(D==1)] = m.oc.10.1.samples$train[,1]
      m.oc.10.1[!((G==1)&(D==1))] = m.oc.10.1.samples$test[,1]
      
      rate_101 = dd + 0.5*(sum((Y[(G==1)&(D==1)]-m.oc.10.1[(G==1)&(D==1)])**2))
      sigma2_101 = rgamma(1, shape=(cc+sum((G==1)&(D==1))/2), rate=rate_101)**(-1)
    }
    M.oc.10.1[,i] = m.oc.10.1
    SIGMA_101[i] = sigma2_101
    
    # potential outcomes
    Y1 = m.oc.11.1[G==2]
    Y0 = m.oc.11.0[G==2]
    te = mean(Y1-Y0)  
    TE[i] = te 
    
    # assessment metrics
    tgt.IE = which(G0==2)
    XG2 = df[tgt.IE,c('x1','x2','x3','x4','x5','x6')]
    IE0 = (0.8 + 0.7*cos(pi*XG2[,1]*XG2[,2]) + 
             0.4*(XG2[,3]-0.5)**3 + 0.7*tanh(XG2[,4]*XG2[,5]*XG2[,6])) - 
      (-0.4 + 0.6*sin(pi*XG2[,1]*XG2[,2]) +
         0.1*(XG2[,3]-0.5)**3 + 0.5*tanh(XG2[,4]*XG2[,5]*XG2[,6]))
    IE.est = m.oc.11.1[tgt.IE]-m.oc.11.0[tgt.IE]
    PEHE[,i] = IE.est-IE0
    IE_RB_num[i] = mean(IE.est-IE0)
    IE_RB_den[i] = mean(IE0)
    IE_RM[,i] = IE.est-IE0
    
    # update outcome BART means, m.Z, m.W
    m.Z.train = cbind(Z=Z, df[, c('x1','x2','x3','x4','x5','x6')])
    m.Z.sampler = dbarts(Z ~ x1+x2+x3+x4+x5+x6, m.Z.train, control = db.control.2)
    m.Z.samples = m.Z.sampler$run()
    m.Z = m.Z.samples$train[,1]
    M.Z[,i] = m.Z

    m.W.train = cbind(W=W[G!=0], df[G!=0, c('x1','x2','x3','x4','x5','x6')])
    m.W.test = df[G==0, c('x1','x2','x3','x4','x5','x6')]
    m.W.sampler = dbarts(W ~ x1+x2+x3+x4+x5+x6, m.W.train, m.W.test, control = db.control.2)
    m.W.samples = m.W.sampler$run()
    m.W = rep(NA, N)
    m.W[G!=0] = m.W.samples$train[,1]
    m.W[G==0] = m.W.samples$test[,1]
    M.W[,i] = m.W

    # update G
    G = rep(NA,N)
    G[is.na(Y)&(D==1)] = 0 # D = 1, S(1) = 0
    G[(!is.na(Y))&(D==0)] = 2 # D = 0, S(0) = 1
    # D = 0, S(0) = 0
    tgt.D0S0 = intersect(which(is.na(Y)), which(D==0))
    p00.Z = pnorm(m.Z[tgt.D0S0])
    p10.Z = (1-pnorm(m.Z[tgt.D0S0])) * pnorm(m.W[tgt.D0S0])
    G[tgt.D0S0] = rep(1,length(tgt.D0S0)) - rbinom(length(tgt.D0S0),1,p00.Z/(p00.Z+p10.Z))
    # D = 1, S(1) = 1
    tgt.D1S1 = intersect(which(!is.na(Y)), which(D==1))
    p10.W = pnorm(m.W[tgt.D1S1]) * dnorm(Y[tgt.D1S1], m.oc.10.1[tgt.D1S1], sqrt(sigma2_101))
    p11.W = (1-pnorm(m.W[tgt.D1S1])) * dnorm(Y[tgt.D1S1], m.oc.11.1[tgt.D1S1], sqrt(sigma2_111))
    G[tgt.D1S1] = rep(2,length(tgt.D1S1)) - rbinom(length(tgt.D1S1),1,p10.W/(p10.W+p11.W))
    GV[,i] = G
    GF[,i] = c(sum(G==0),sum(G==1),sum(G==2))/N

    # update Z & W
    Z = rep(NA,N)
    Z[G==0] = rtruncnorm(sum(G==0), a=0, mean=(m.Z[G==0]), sd=1)
    Z[G!=0] = rtruncnorm(sum(G!=0), b=0, mean=(m.Z[G!=0]), sd=1)
    W = rep(NA,N)
    if (sum(G==1) > 0) {
      W[G==1] = rtruncnorm(sum(G==1), a=0, mean=(m.W[G==1]), sd=1)
    }
    if (sum(G==2) > 0) {
      W[G==2] = rtruncnorm(sum(G==2), b=0, mean=(m.W[G==2]), sd=1)
    }
    W[G==0] = NA

    if (i%%100==0) {
      cat("\r", paste("sampling iteration", i, sep=" "))
      flush.console()
    }
  }   
  return(list(SIGMA_111=SIGMA_111, SIGMA_110=SIGMA_110, SIGMA_101=SIGMA_101, TE=TE, GV=GV, GF=GF,
              M.oc.11.1=M.oc.11.1, M.oc.11.0=M.oc.11.0, M.oc.10.1=M.oc.10.1,
              M.Z=M.Z, M.W=M.W, PEHE=PEHE, IE_RB_num=IE_RB_num, IE_RB_den=IE_RB_den, IE_RM=IE_RM))
}

n_rep = 1
sigma_111_rep = matrix(NA, nrow=n_rep, ncol=3)
sigma_110_rep = matrix(NA, nrow=n_rep, ncol=3)
sigma_101_rep = matrix(NA, nrow=n_rep, ncol=3)
TE_rep = matrix(NA, nrow=n_rep, ncol=4)
IE_rep = matrix(NA, nrow=n_rep, ncol=4)
GF_rep = array(NA, dim=c(n_rep,4,3))
n_burn = 3000
n_iter = 3000

# test code
set.seed(12345)
for (ii in 1:n_rep) {
  df = gendata()
  S1 = combined.Gibbs(df=df, S=n_burn+n_iter)
  sigma_111_rep[ii,1] = mean(S1$SIGMA_111[-(1:n_burn)])
  sigma_111_rep[ii,2] = quantile(S1$SIGMA_111[-(1:n_burn)], probs=0.025)
  sigma_111_rep[ii,3] = quantile(S1$SIGMA_111[-(1:n_burn)], probs=0.975)
  
  sigma_110_rep[ii,1] = mean(S1$SIGMA_110[-(1:n_burn)])
  sigma_110_rep[ii,2] = quantile(S1$SIGMA_110[-(1:n_burn)], probs=0.025)
  sigma_110_rep[ii,3] = quantile(S1$SIGMA_110[-(1:n_burn)], probs=0.975)
  
  sigma_101_rep[ii,1] = mean(S1$SIGMA_101[-(1:n_burn)])
  sigma_101_rep[ii,2] = quantile(S1$SIGMA_101[-(1:n_burn)], probs=0.025)
  sigma_101_rep[ii,3] = quantile(S1$SIGMA_101[-(1:n_burn)], probs=0.975)
  
  TE_rep[ii,1] = df$TE[1]
  TE_rep[ii,2] = mean(S1$TE[-(1:n_burn)])
  TE_rep[ii,3] = quantile(S1$TE[-(1:n_burn)], probs=0.025)
  TE_rep[ii,4] = quantile(S1$TE[-(1:n_burn)], probs=0.975)  
  
  IE_rep[ii,1] = sqrt(mean((apply(S1$PEHE[,-(1:n_burn)], 1, mean))**2)) # PEHE
  IE_rep[ii,2] = mean(S1$IE_RB_num[-(1:n_burn)]) # IE_RB_num
  IE_rep[ii,3] = mean(S1$IE_RB_den[-(1:n_burn)]) # IE_RB_den
  IE_rep[ii,4] = mean((apply(S1$IE_RM[,-(1:n_burn)], 1, mean))**2) # IE_RM
  
  GF_rep[ii,1,] = table(df$G)/nrow(df)
  GF_rep[ii,2,] = apply(S1$GF[,-(1:n_burn)], 1, mean)
  GF_rep[ii,3,] = apply(S1$GF[,-(1:n_burn)], 1, function(x) quantile(x, probs=0.025))
  GF_rep[ii,4,] = apply(S1$GF[,-(1:n_burn)], 1, function(x) quantile(x, probs=0.975))
  
  if (ii%%1==0) {
    cat("\r", paste("repetition", n_rep, "iteration", ii, sep=" "))
    flush.console()
  }
}

# sigma
mean(sigma_111_rep[,1]) # mean estimate
mean(sigma_111_rep[,1]-2)/2 # relative bias
sqrt(mean((sigma_111_rep[,1]-2)**2)) # RMSE
sum((sigma_111_rep[,2] <= 2)*(sigma_111_rep[,3] >= 2))/n_rep # coverage percentage

mean(sigma_110_rep[,1]) # mean estimate
mean(sigma_110_rep[,1]-2)/2 # relative bias
sqrt(mean((sigma_110_rep[,1]-2)**2)) # RMSE
sum((sigma_110_rep[,2] <= 2)*(sigma_110_rep[,3] >= 2))/n_rep # coverage percentage

mean(sigma_101_rep[,1]) # mean estimate
mean(sigma_101_rep[,1]-2)/2 # relative bias
sqrt(mean((sigma_101_rep[,1]-2)**2)) # RMSE
sum((sigma_101_rep[,2] <= 2)*(sigma_101_rep[,3] >= 2))/n_rep # coverage percentage

# SACE
mean(TE_rep[,1]) # truth
mean(TE_rep[,2]) # mean estimate
(mean(TE_rep[,2])-mean(TE_rep[,1]))/mean(TE_rep[,1]) # relative bias
sqrt(mean((TE_rep[,2]-TE_rep[,1])**2)) # RSME
sum((TE_rep[,3] <= TE_rep[,1])*(TE_rep[,4] >= TE_rep[,1]))/n_rep # coverage percentage

# ISCE
mean(IE_rep[,1]) # PEHE
mean(IE_rep[,2])/mean(IE_rep[,3]) # relative bias
sqrt(mean(IE_rep[,4])) # RMSE

# GF
apply(GF_rep[,1,], 2, mean)
apply(GF_rep[,2,], 2, mean)
apply((GF_rep[,2,]-GF_rep[,1,])/GF_rep[,1,], 2, mean)
apply((GF_rep[,3,] <= GF_rep[,1,])*(GF_rep[,4,] >= GF_rep[,1,]), 2, sum)/n_rep

