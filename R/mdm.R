#' Calculate the log predictive likelihood for a specified set of parents and a fixed delta.
#'
#' @param Yt the vector of observed time series, length T
#' @param Ft the matrix of covariates, dim = number of thetas (p) x number of time points (T), usually a row of 1s to represent an intercept and the time series of the parent nodes.
#' @param delta discount factor (scalar).
#' @param m0 the value of the prior mean at time t=0, scalar, and assuming the mean is the same for all nodes. The default is zero. (theta_{0} | D_{0}, phi) ~ N(m_{0},C*_{0} x phi^-1), D_{0} denotes the set of initial information.
#' @param CS0 controls the scaling of the prior variance matrix C*_{0} at time t=0. The default is 3, giving a non-informative prior for C*_{0}, 3 x (p x p) identity matrix.
#' @param n0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. n0 has to be higher than 0.
#' @param d0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. 
#'
#' @return
#' mt the vector or matrix of the posterior mean (location parameter), dim = p x T.
#' Ct the posterior scale matrix, dim = p x p x T, C_{t} = C*_{t} x S_{t}, where S_{t} is a point estimate for the observation variance phi^-1. 
#' CSt the posterior scale matrix, dim = p x p x T, C_{t} = C*_{t} x S_{t}, where S_{t} is a point estimate for the observation variance phi^-1.
#' Rt the prior scale matrix, dim = p x p x T. R_{t} = R*_{t} x S_{t-1}, where S_{t-1} is a point estimate for the observation variance phi^-1 at the previous time point.
#' RSt the prior scale matrix, dim = p x p x T. R_{t} = R*_{t} x S_{t-1}, where S_{t-1} is a point estimate for the observation variance phi^-1 at the previous time point.
#' nt and dt the vectors of the hyperparameters for the precision phi with length T.
#' S the vector of the point estimate for the observation variance phi^-1 with length T.
#' ft the vector of the one-step forecast location parameter with length T.
#' Qt the vector of the one-step forecast scale parameter with length T.
#' ets the vector of the standardised forecast residuals with length T, defined as Y_{t} - f_{t} / sqrt (Q_{t}).
#' lpl the vector of the Log Predictive Likelihood with length T.
#' 
#' @export
dlm.lpl <- function(Yt, Ft, delta, m0 = 0, CS0 = 3, n0 = 0.001, d0 = 0.001){
  
  CS0 = CS0*diag(nrow(Ft))
  m0 = rep(m0,nrow(Ft))
  Nt = length(Yt)+1 # the length of the time series + t0
  p = nrow(Ft)      # the number of parents and one for an intercept (i.e. the number of thetas)
  
  Y = numeric(Nt)
  Y[2:Nt] = Yt
  
  F1 = array(0, dim=c(p,Nt))
  F1[,2:Nt] = Ft
  
  # Set up allocation matrices, including the priors
  mt = array(0,dim=c(p,Nt))
  mt[,1]=m0
  
  Ct = array(0,dim=c(p,p,Nt))
  CSt = array(0,dim=c(p,p,Nt))
  CSt[,,1] = CS0
  
  Rt = array(0,dim=c(p,p,Nt))
  RSt = array(0,dim=c(p,p,Nt))
  
  nt = numeric(Nt) 
  nt[1]=n0
  
  dt = numeric(Nt)
  dt[1]=d0
  
  S = numeric(Nt)
  S[1]=dt[1]/nt[1]
  
  ft = numeric(Nt)
  Qt = numeric(Nt)
  ets = numeric(Nt)
  lpl = numeric(Nt)
  
  # Updating
  
  for (t in 2:Nt){
    
    # Posterior at {t-1}: (theta_{t-1}|D_{t-1}) ~ T_{n_{t-1}}[m_{t-1}, C_{t-1} = C*_{t-1} x d_{t-1}/n_{t-1}]
    # Prior at {t}: (theta_{t}|D_{t-1}) ~ T_{n_{t-1}}[m_{t-1}, R_{t}]
    # D_{t-1} = D_{0},Y_{1},...,Y_{t-1} D_{0} is the initial information set
    
    # R*_{t} = C*_{t-1}/delta
    RSt[,,t] = CSt[,,(t-1)] / delta
    Rt[,,t] = RSt[,,t] * S[(t-1)] 
    # One-step forecast: (Y_{t}|D_{t-1}) ~ T_{n_{t-1}}[f_{t}, Q_{t}]
    ft[t] = t(F1[,t]) %*% mt[,(t-1)]
    QSt = as.vector(1 + t(F1[,t]) %*% RSt[,,t] %*% F1[,t])
    Qt[t] = QSt * S[(t-1)]
    et = Y[t] - ft[t]
    ets[t] = et / sqrt(Qt[t])
    
    # Posterior at t: (theta_{t}|D_{t}) ~ T_{n_{t}}[m_{t}, C_{t}]
    # D_{t} = D_{0},Y_{1},...,Y_{t}
    At = (RSt[,,t] %*% F1[,t])/QSt
    mt[,t] = mt[,(t-1)] + (At*et)
    
    nt[t] = nt[(t-1)] + 1
    dt[t] = dt[(t-1)] + (et^2)/QSt
    S[t]=dt[t]/nt[t] 
    
    CSt[,,t] = RSt[,,t] - (At %*% t(At))*QSt
    Ct[,,t] = S[t]*CSt[,,t]
    
    # Log Predictive Likelihood 
    lpl[t] = lgamma((nt[(t-1)]+1)/2)-lgamma(nt[(t-1)]/2)-0.5*log(pi*nt[(t-1)]*Qt[t])-((nt[(t-1)]+1)/2)*log(1+(1/nt[(t-1)])*et^2/Qt[t])
  }
  
  mt = mt[,2:Nt]; Ct = Ct[,,2:Nt]; CSt = CSt[,,2:Nt]; Rt = Rt[,,2:Nt]; RSt = RSt[,,2:Nt]
  nt = nt[2:Nt]; dt = dt[2:Nt]; S = S[2:Nt]; ft = ft[2:Nt]; Qt = Qt[2:Nt]; ets = ets[2:Nt]; lpl = lpl[2:Nt]
  
  output <- list(mt=mt,Ct=Ct,CSt=CSt,Rt=Rt,RSt=RSt,nt=nt,dt=dt,S=S,ft=ft,Qt=Qt,ets=ets,lpl=lpl)
  return(output)
}

#' A function to generate all the possible models. 
#'
#' @param Nn number of nodes; the number of columns of the dataset can be used.
#' @param node the node of interest (i.e., the node to find parents for).
#'
#' @return
#' output.model = a matrix with dimensions (Nn-1) x number of models, where number of models = 2^(Nn-1).
#' 
#' @export
model.generator<-function(Nn,node){
  
  # Create the model 'no parents' (the first column of the matrix is all zeros)
  empt=rep(0,(Nn-1)) 
  
  for (k in 1:(Nn-1)) {
    
    # Calculate all combinations when number of parents = k
    #m=combn(c(1:Nn)[-node],k)
    if (Nn==2 & node==1) {
      model = matrix(c(0,2),1,2)
    } else { 
      m=combn(c(1:Nn)[-node],k) 
      
      # Expand the array so that unconnected edges are represented by zeros  
      empt.new=array(0,dim=c((Nn-1),ncol(m)))
      empt.new[1:k,]=m
      
      # Bind the matrices together; the next set of models are added to this matrix
      model=cbind(empt,empt.new)
      empt=model
    } 
  }
  
  colnames(model)=NULL
  output.model<-model
  
  return(output.model)
  
}

#' A function for an exhaustive search, calculates the optimum value of the discount factor.
#'
#' @param Data  Dataset with dimension number of time points T x Number of nodes Nn.
#' @param node  node of interest.
#' @param nbf   Log Predictive Likelihood will sum from (and including) this time point. 
#' @param delta a vector of potential values for the discount factor.
#' @param cpp boolean true (default): fast C++ implementation, false: native R code.
#' @param m0 the value of the prior mean at time t=0, scalar, and assuming the mean is the same for all nodes. The default is zero. (theta_{0} | D_{0}, phi) ~ N(m_{0},C*_{0} x phi^-1), D_{0} denotes the set of initial information.
#' @param CS0 controls the scaling of the prior variance matrix C*_{0} at time t=0. The default is 3, giving a non-informative prior for C*_{0}, 3 x (p x p) identity matrix.
#' @param n0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. n0 has to be higher than 0.
#' @param d0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. 
#'
#' @return
#' model.store a matrix with the model, LPL and chosen discount factor for all possible models.
#' runtime an estimate of the run time of the function, using proc.time().
#' @export
exhaustive.search <- function(Data, node, nbf=15, delta=seq(0.5,1,0.01), cpp=TRUE, m0 = 0, CS0 = 3, n0 = 0.001, d0 = 0.001) {
  
  ptm=proc.time()  
  
  Nn=ncol(Data) # the number of nodes
  Nm=2^(Nn-1)   # the number of models per node
  
  M=model.generator(Nn,node) # Generate all the possible models
  models=rbind(1:Nm,M) # Label each model with a 'model number'
  
  
  Yt=Data[,node]    # the time series of the node we wish to find parents for
  Nt=length(Yt)     # the number of time points
  nd=length(delta)  # the number of deltas
  
  # Create empty arrays for the lpl scores and the optimum deltas
  lpldet=array(NA,c(Nm,length(delta)))
  lplmax=rep(NA,Nm)
  DF.hat=rep(NA,Nm)
  
  # Now create Ft. 
  for (z in 1:Nm) {
    pars=models[(2:Nn),z] 
    pars=pars[pars!=0]
    Ft=array(1,dim=c(Nt,length(pars)+1))
    if (ncol(Ft)>1) {
      Ft[,2:ncol(Ft)]=Data[,pars] # selects parents 
    }  
    
    # Calculate the log predictive likelihood, for each value of delta
    for (j in 1:nd) {
      if (cpp) {
        # new C++ implementation
        lpl=c(dlmLplCpp(Yt, t(Ft), delta[j], m0, CS0, n0, d0))
        lpldet[z,j]=sum(lpl[nbf:Nt])
      } else {
        # original native R
        a=dlm.lpl(Yt, t(Ft), delta=delta[j], m0=m0, CS0=CS0, n0=n0, d0=d0)
        lpldet[z,j]=sum(a$lpl[nbf:Nt])
      }
    }
    
    lplmax[z]=max(lpldet[z,],na.rm=TRUE)
    DF.hat[z]=delta[lpldet[z,]==max(lpldet[z,],na.rm=TRUE)] # add which here for index
  }
  
  # Output model.store
  model.store=rbind(models,lplmax,DF.hat)
  rownames(model.store)=NULL
  
  runtime=(proc.time()-ptm)
  
  output<-list(model.store=model.store,runtime=runtime)    
  return(output)
}

#' Mean centers timeseries in a 2D array timeseries x nodes,
#' i.e. each timeseries of each node has mean of zero.
#'
#' @param X 2D array with dimensions timeseries x nodes.
#'
#' @return M 2D array.
#' @export
center <- function(X) {
  d = dim(X)
  M = matrix(NA, d[1], d[2])
  
  for (i in 1:d[2]) {
    M[,i]=scale(X[,i], center = T, scale = F)
  }
  
  return(M)
}

#' Estimate subject's full network: runs exhaustive search on very node.
#'
#' @param X array with dimensions timeseries x nodes.
#' @param id subject ID. If set, results are saved to a txt file.
#' @param nbf  Log Predictive Likelihood will sum from (and including) this time point. 
#' @param delta a vector of potential values for the discount factor.
#' @param cpp boolean true (default): fast C++ implementation, false: native R code.
#' @param m0 the value of the prior mean at time t=0, scalar, and assuming the mean is the same for all nodes. The default is zero. (theta_{0} | D_{0}, phi) ~ N(m_{0},C*_{0} x phi^-1), D_{0} denotes the set of initial information.
#' @param CS0 controls the scaling of the prior variance matrix C*_{0} at time t=0. The default is 3, giving a non-informative prior for C*_{0}, 3 x (p x p) identity matrix.
#' @param n0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. n0 has to be higher than 0.
#' @param d0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. 
#'
#' @return store list with results.
#' @export
subject <- function(X, id=NULL, nbf=15, delta=seq(0.5,1,0.01), cpp=TRUE, m0 = 0, CS0 = 3, n0 = 0.001, d0 = 0.001) {
  N=ncol(X)  # nodes
  M=2^(N-1)  # rows/models
  models = array(rep(0,(N+2)*M*N),dim=c(N+2,M,N))
  
  for (n in 1:N) {
    tmp=exhaustive.search(X, n, nbf=nbf, delta=delta, cpp=cpp, m0=m0, CS0=CS0, n0=n0, d0=d0)
    models[,,n]=tmp$model.store
    if (!is.null(id)) {
      write(t(models[,,n]), file=sprintf("%s_node_%03d.txt", id, n), ncolumns = M)
    }
  }
  
  store=list()
  store$models=models
  store$winner=getWinner(models,N)
  store$adj=getAdjacency(store$winner,N)
  store$thr=getThreshAdj(store$adj, store$models, store$winner)
  
  return(store)
}

#' Runs exhaustive search on a single node and saves results in txt file.
#'
#' @param X array with dimensions timeseries x nodes.
#' @param n node number.
#' @param id subject ID. If set, results are saved to a txt file.
#' @param nbf  Log Predictive Likelihood will sum from (and including) this time point. 
#' @param delta a vector of potential values for the discount factor.#'
#' @param cpp boolean true (default): fast C++ implementation, false: native R code.
#' @param m0 the value of the prior mean at time t=0, scalar, and assuming the mean is the same for all nodes. The default is zero. (theta_{0} | D_{0}, phi) ~ N(m_{0},C*_{0} x phi^-1), D_{0} denotes the set of initial information.
#' @param CS0 controls the scaling of the prior variance matrix C*_{0} at time t=0. The default is 3, giving a non-informative prior for C*_{0}, 3 x (p x p) identity matrix.
#' @param n0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. n0 has to be higher than 0.
#' @param d0 prior hypermarameter of precision phi ~ G(n_{0}/2; d_{0}/2). The default is a non-informative prior, with n0 = d0 = 0.001. 
#' 
#' @return store list with results.
#' @export
node <- function(X, n, id=NULL, nbf=15, delta=seq(0.5,1,0.01), cpp=TRUE, m0 = 0, CS0 = 3, n0 = 0.001, d0 = 0.001) {
  N=ncol(X)  # nodes
  M=2^(N-1)  # rows/models
  
  store=exhaustive.search(X, n, nbf=nbf, delta=delta, cpp=cpp, m0=m0, CS0=CS0, n0=n0, d0=d0)
  if (!is.null(id)) {
    write(t(store$model.store), file=sprintf("%s_node_%03d.txt", id, n), ncolumns = M)
  }
  return(store)
}

#' Reads single subject's network from txt files.
#' @param path path.
#' @param id identifier to select all subjects' nodes, e.g. pattern containing subject ID and session number.
#' @param nodes number of nodes.
#'
#' @return store list with results.
#' @export
read.subject <- function(path, id, nodes) {
  models = array(0,dim=c(nodes+2,2^(nodes-1),nodes))
  for (n in 1:nodes) {
    #file=sprintf("%s_node_%03d.txt", id, n)
    #models[,,n] = as.matrix(read.table(file)) # quite slow
    file=list.files(path, pattern=glob2rx(sprintf("%s*_node_%03d.txt", id, n)))
    models[,,n] = as.matrix(fread(file.path(path,file))) # faster, from package "data.table"
  }
  store=list()
  store$models=models
  store$winner=getWinner(models,nodes)
  store$adj=getAdjacency(store$winner,nodes)
  store$thr=getThreshAdj(store$adj, store$models, store$winner)
  
  return(store)
}

#' Get winner network by maximazing log predictive likelihood (LPL)
#' from a set of models.
#'
#' @param models 2D matrix, or 3D models x node.
#' @param nodes number of nodes.
#'
#' @return winner array with highest scored model(s).
#' @export
getWinner <- function(models, nodes) {
  
  dims=length(dim(models))
  
  if (dims==2) {
    winner = models[,which.max(models[nodes+1,])]
  } else if (dims==3) {
    winner = array(0, dim=c(nodes+2,nodes))
    for (n in 1:nodes) {
      winner[,n]=models[,which.max(models[nodes+1,,n]),n]
    }
  }
  return(winner)
}

#' Get adjacency and associated likelihoods (LPL) and disount factros (df) of winning models.
#'
#' @param winner, 2D matrix.
#' @param nodes number of nodes.
#'
#' @return adj, 2D adjacency matrix.
#' @export
getAdjacency <- function(winner, nodes) {
  
  am = array(rep(0,nodes*nodes),dim=c(nodes,nodes))
  lpl = df = array(rep(NA,nodes*nodes),dim=c(nodes,nodes))
  for (n in 1:nodes) {
    p = winner[2:nodes,n]  # parents
    p = p[p>0]
    am[p,n] = 1
    lpl[p,n] = winner[nodes+1,n]
    df[p,n]  = winner[nodes+2,n]
  }
  return(list(am=am, lpl=lpl, df=df))
}

#' Plots network as adjacency matrix.
#'
#' @param adj 2D adjacency matrix.
#' @param title title.
#' @param label label for colormap.
#' @param hasColMap FALSE turns off color map, default is NULL (on).
#' @param lim vector with min and max value for color scaling.
#'
#' @export
gplotMat <- function(adj, title=NULL, label=NULL, hasColMap=NULL, lim=c(0, 1)) {
  x = melt(adj)
  names(x)[1] = "Parent"
  names(x)[2] = "Child"
  
  ggplot(x, aes_string(x = "Child", y = "Parent", fill = "value")) +
    geom_tile(color = "gray60") +

    scale_fill_gradient2(
      low = "white",
      high = "red",
      mid = "orange",
      midpoint = sum(lim)/2,
      limit = lim,
      space = "Lab",
      name = label) + 

    theme(#axis.ticks.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          text = element_text(size=12),
          axis.text =  element_text(size=12),
          plot.title = element_text(size=12)
          ) +
    scale_y_reverse() + ggtitle(title) + guides(fill=hasColMap)
}

#' Performes a binomial test with FDR correction for network edges in an adjacency matrix.
#'
#' @param adj adjacency matrix, nodes x nodes x subj, or nodes x nodes x runs x subj.
#' @param alter type of binomial test, "two.sided" (default), "less", or "greater"
#'
#' @return store list with results.
#' @export
binom.nettest <- function(adj, alter="two.sided", fdr=0.05) {
  
  mydim=dim(adj)
  M = sum(adj) # total edges over all N subjects, all R(R-1) edges
  
  if (length(mydim) == 3) { # without runs
    N=mydim[3] # No. of subjects
    N_Comp=mydim[1]
    adj_ = apply(adj, c(1,2), sum)
    
  } else if (length(mydim) == 4) { # with runs
    N=mydim[4] # No. of subjects
    N_runs=mydim[3]
    N_Comp=mydim[1]
    N=N*N_runs # adjust N by for no. of runs
    adj_ = apply(adj, c(1,2), sum) # sum acrosss subjects and runs
  }
  
  # binom test for every edge occurance
  p0 = M/N/N_Comp/(N_Comp-1) # H0 edge probability
  
  p = array(NA, dim=c(N_Comp,N_Comp))
  for (i in 1:N_Comp) {
    for (j in 1:N_Comp) {
      tmp=binom.test(adj_[i,j],N,p=p0, alternative=alter)
      p[i,j]=tmp$p.value
    }
  }
  
  # FDR
  p_fdr=matrix(p.adjust(p, method = "fdr"),N_Comp,N_Comp)
  adj_fdr=adj_
  adj_fdr[p_fdr>=fdr]=NA
  
  store=list()
  store$p0=p0
  store$p=p
  store$p_fdr=p_fdr
  store$adj=adj_/N
  store$adj_fdr=adj_fdr/N # significant proportions
  
  return(store)
}

#' Reshapes a 2D concatenated time series into 3D according to no. of subjects and volumes.
#'
#' @param ts a 2D time series volumes x nodes.
#' @param N No. of subjects.
#' @param V No. of volumes.
#' 
#' @return M 3D matrix, time series x nodes x subjects.
#' @export
reshapeTs <- function(ts, N, V) {
  NC = ncol(ts)
  M = array(NA, dim=c(V,NC,N))
  for (i in 1:N) {
    idx = ((i-1)*V+1):(V*i)
    M[,,i] = ts[idx,]
  }
  return(M)
}

#' Correlation of time series.
#'
#' @param ts a 3D time series time series x nodes x subjects.
#' 
#' @return M correlation matrix.
#' @export
corTs <- function(ts) {
  d=dim(ts)
  N=d[3] # No. subjects
  N_nodes=d[2]
  R=array(NA, dim=c(N_nodes,N_nodes,N))
  for (s in 1:N) {
    R[,,s]=cor(ts[,,s])
  }
  M = apply(R, c(1,2), mean)
  return(M)
}

#' Get specific parent model from all models.
#'
#' @param models a 2D model matrix.
#' @param parents a vector with parent nodes.
#' 
#' @return mod specific parent model.
#' @export
getModel <- function(models, parents) {
  Nn = nrow(models) - 2 # No. of nodes
  Nm = ncol(models) # No. of models
  parents = c(parents, rep(0, (Nn - 1) - length(parents))) # add fill zeros
  for (m in 1:Nm) {
    if (sum(models[2:Nn,m] == parents) == Nn - 1) {
      mod = models[,m]
    }
  }
  return(mod)
}

#' A group is a list containing restructured data from subejcts for easier group analysis.
#'
#' @param subj a list of subjects.
#' 
#' @return group a list.
#' @export
mdm.group <- function(subj) {
  Nn=ncol(subj[[1]]$adj$am)
  N=length(subj)
  
  am = lpl = df = tam = tbi = array(NA, dim=c(Nn,Nn,N))
  tlpls = array(NA, dim=c(Nn,Nn,2,N))
  for (s in 1:N) {
    am[,,s]  = subj[[s]]$adj$am
    lpl[,,s] = subj[[s]]$adj$lpl
    df[,,s]  = subj[[s]]$adj$df
    
    # thresholded measures
    tam[,,s]  = subj[[s]]$thr$am
    tbi[,,s]  = subj[[s]]$thr$bi
    tlpls[,,,s]= subj[[s]]$thr$lpls
  }
  
  group=list(am=am,lpl=lpl,df=df,tam=tam,tbi=tbi,tlpls=tlpls)
  return(group)
}

#' A group is a list containing restructured data from subejcts for easier group analysis.
#'
#' @param subj a list of subjects.
#' 
#' @return group a list.
#' @export
patel.group <- function(subj) {
  Nn=ncol(subj[[1]]$kappa)
  N=length(subj)
  
  kappa = tkappa = tau = ttau = net = tnet = array(NA, dim=c(Nn,Nn,N))
  for (s in 1:N) {
    kappa[,,s]  = subj[[s]]$kappa
    tkappa[,,s] = subj[[s]]$tkappa
    tau[,,s]    = subj[[s]]$tau
    ttau[,,s]   = subj[[s]]$ttau
    net[,,s]    = subj[[s]]$net
    tnet[,,s]   = subj[[s]]$tnet
  }
  
  group=list(kappa=kappa,tkappa=tkappa,tau=tau,ttau=ttau,net=net,tnet=tnet)
  return(group)
}

#' Get thresholded adjacency network.
#'
#' @param adj list with network adjacency from getAdjacency().
#' @param models matrix 3D with full model estimates.
#' @param winner matrix 2D with winning models.
#' 
#' @return thr list with thresholded network adjacency.
#' @export
getThreshAdj <- function(adj, models, winner) {
  
  Nn = ncol(adj$am)
  # determine bidirectional edges
  bi = array(0, dim=c(Nn, Nn))
  for (i in 1:Nn) {
    for (j in 1:Nn) {
      if (adj$am[i,j] == 1 & adj$am[j,i] == 1) {
        bi[i,j] = 1
      }
    }
  }
  
  # Calculate models
  B=bi*upper.tri(bi)
  lpls=array(NA, dim=c(Nn, Nn, 2))
  for (i in 1:Nn) {
    for (j in 1:Nn) {
      # if bidirectional, calculate 3 models
      # A+B:i<>j, A:i>j, B:i<j
      # A+B: LPLj + LPLi
      # A  : LPLj + LPLi-j
      # B  : LPLi + LPLj-i
      if (B[i,j] == 1) {
        
        # bidirectional LPL
        lpls[i,j,1] = lpls[j,i,1] = adj$lpl[i,j] + adj$lpl[j,i]
        
        # uni i->j
        p = winner[,i][2:Nn]
        p = p[p != j & p!= 0] # remove node j
        lpls[i,j,2] = adj$lpl[i,j] + getModel(models[,,i], p)[Nn+1]
        
        # uni j->i
        p = winner[,j][2:Nn]
        p = p[p != i & p!= 0] # remove node i
        lpls[j,i,2] = adj$lpl[j,i] + getModel(models[,,j], p)[Nn+1]
      }
    }
  }
  
  # am matrix
  am=adj$am
  BF=20 # bayes factor threshold
  for (i in 1:Nn) {
    for (j in 1:Nn) {
      if (B[i,j] == 1) {
        if (lpls[i,j,1] - BF <= max(lpls[i,j,2], lpls[j,i,2]) ) {
          if (lpls[i,j,2] > lpls[j,i,2]) {
            am[i,j] = 1; am[j,i] = 0
          } else {
            am[i,j] = 0; am[j,i] = 1
          }
        }
      }
    }
  }
    
  thr=list()
  thr$bi=bi # bidirectional edges
  thr$lpls=lpls # lpls
  thr$am=am # adjacency matrix (thresholded)
  
  return(thr)
}

#' Performance of estimates, such as sensitivity, specificity, and more.
#'
#' @param x estimated binary network matrix.
#' @param xtrue, true binary network matrix.
#' 
#' @return perf vector.
#' @export
perf <- function(x, xtrue) {
  
  # in case xtrue continas NA instead of 0
  xtrue[is.na(xtrue)]=0
  
  d = dim(x)
  Nn=d[1]
  if (length(d) == 3) {
    N=d[3]
  } else if (length(d) == 2) {
    x=array(x, dim=c(Nn,Nn,1))
    N=1
  }
  
  perf=array(NA,dim=c(N,8))
  
  for (i in 1:N) {
    # see https://en.wikipedia.org/wiki/Sensitivity_and_specificity
    TP = sum(x[,,i] & xtrue)
    FP = sum((x[,,i] - xtrue) == 1)
    FN = sum((xtrue - x[,,i]) == 1)
    TN = sum(!x[,,i] & !xtrue) - ncol(x[,,i])
    
    tpr = TP/(TP+FN) # 1
    spc = TN/(TN+FP) # 2
    ppv = TP/(TP+FP) # 3
    npv = TN/(TN+FN) # 4
    fpr = FP/(FP+TN) # 5
    fnr = FN/(TP+FN) # 6
    fdr = FP/(TP+FP) # 7
    acc = (TP+TN)/(TP+FP+FN+TN) # 8
    
    perf[i,]=c(tpr,spc,ppv,npv,fpr,fnr,fdr,acc)
  }
  
  return(perf)
}
  
#' Scaling data. Zero centers and scales the nodes (SD=1).
#'
#' @param X time x node 2D matrix, or 3D with subjects as the 3rd dimension.
#'
#' @return S centered and scaled matrix.
#' @export
scaleTs <- function(X) {
  D=dim(X)
  if (length(D)==2) {
    tmp=array(NA, dim=c(D[1],D[2],1))
    tmp[,,1]=X
    X=tmp
    D=dim(X)
  }
  
  S=array(NA, dim=c(D[1],D[2],D[3]))
  
  # center data
  for (s in 1:D[3]) {
    colm=colMeans(X[,,s])
    S[,,s]=X[,,s]-t(array(rep(colm, D[1]), dim=c(D[2],D[1])))
  }
  
  # scale data
  for (s in 1:D[3]) {
    SD=sqrt(mean(apply(S[,,s], 2, var)))
    S[,,s]=S[,,s]/SD
  }
  
  if (D[3]==1) {
    S=S[,,1]
  }
  return(S)
}

#' Patel.
#'
#' @param X time x node 2D matrix.
#' @param lower percentile cuttoff.
#' @param upper percentile cuttoff for 0-1 scaling.
#' @param bin threshold for conversion to binary values.
#' @param TK significance threshold for connection strength kappa.
#' @param TT significance threshold for direction tau.
#'
#' @return PT list with strengths kappa, direction tau, and net structure.
#' @export
patel <- function(X, lower=0.1, upper=0.9, bin=0.75, TK=0, TT=0) {
  
  nt=nrow(X)
  nn=ncol(X)
  
  # scale data into 0,1 interval
  X10 = apply(X, 2, quantile, lower) # cutoff 0.1 percentile
  X90 = apply(X, 2, quantile, upper) # cutoff 0.9 percentile
  a = sweep(sweep(X, 2, X10), 2, X90-X10, FUN="/") # center, scale data
  X2 = apply(a, c(1,2), function(v) max(min(v,1),0)) # keep between 0,1
  
  # binarize
  X2=(X2>bin)*1 # convert to double
  # Joint activation probability of timeseries a, b
  # See Table 2, Patel et al. 2006
  theta1 = crossprod(X2)/nt        # a=1, b=1 -- a and b active
  theta2 = crossprod(X2,1-X2)/nt   # a=1, b=0 -- a active, b not
  theta3 = crossprod(1-X2,X2)/nt   # a=0, b=1
  theta4 = crossprod(1-X2,1-X2)/nt # a=0, b=0
  
  # directionality tau [-1, 1]
  tau = matrix(0, ncol(X2), ncol(X2))
  inds = theta2 >= theta3
  tau[inds] = 1 - (theta1[inds] + theta3[inds])/(theta1[inds] + theta2[inds])
  tau[!inds] = (theta1[!inds] + theta2[!inds])/(theta1[!inds] + theta3[!inds]) - 1
  tau=-tau
  tau[as.logical(diag(nn))]=NA
  # tau(a,b) positive, a is ascendant to b (a is parent)
  
  # functional connectivity kappa [-1, 1]
  E=(theta1+theta2)*(theta1+theta3)
  max_theta1=min(theta1+theta2,theta1+theta3)
  min_theta1=max(0,2*theta1+theta2+theta3-1)
  inds = theta1>=E
  D = matrix(0, ncol(X2), ncol(X2))
  D[inds]=0.5+(theta1[inds]-E[inds])/(2*(max_theta1-E[inds]))
  D[!inds]=0.5-(theta1[!inds]-E[!inds])/(2*(E[!inds]-min_theta1))
  
  kappa=(theta1-E)/(D*(max_theta1-E) + (1-D)*(E-min_theta1))
  kappa[as.logical(diag(nn))]=NA
  
  # theresholding
  tkappa = kappa
  tkappa[kappa >= TK[1] & kappa <= TK[2]] = 0
  ttau = tau
  ttau[tau >= TT[1] & tau <= TT[2]] = 0
  
  # binary nets: we focus on positive associations
  net = (tkappa > 0)*(tau > 0)  # only kappa thresholded
  net[diag(nn)==1]=0 # remove NA
  tnet = (tkappa > 0)*(ttau > 0) # both tau and kappa thresholded
  tnet[diag(nn)==1]=0 # remove NA
  
  PT=list(kappa=kappa, tau=tau, tkappa=tkappa, ttau=ttau, net=net, tnet=tnet)
  return(PT)
}

#' Permutation test for Patel's kappa. Creates a distribution of values
#' kappa under the null hypothesis.
#'
#' @param X time x node x subjects 3D matrix.
#' @param alpha sign. level
#'
#' @return stat lower and upper significance thresholds.
#' @export
perm.test <- function(X, alpha=0.05) {
  
  low = alpha/2      # two sided test
  up  = 1 - alpha/2
  
  N  = dim(X)[3] # Nr. of subjects
  Nn = dim(X)[2] # Nr. of nodes
  
  # shuffle across subjects with fixed nodes
  ka = array(NA,dim=c(Nn,Nn,N)) # kappa null distribution
  ta = array(NA,dim=c(Nn,Nn,N)) # kappa null distribution
  X_= array(NA, dim=dim(X))
  for (s in 1:N) {
    for (n in 1:Nn) {
      X_[,n,s]=X[,n,sample(N,1)] # draw a random subject (with repetition)
    }
    p=patel(X_[,,s])
    ka[,,s]=p$kappa
    ta[,,s]=p$tau
  }
  
  # two sided sign. test
  stat = list()
  stat$kappa = quantile(ka[!is.na(ka)], probs=c(low, up)) # two sided
  stat$tau   = quantile(ta[!is.na(ta)], probs=c(low, up)) # two sided
  
  return(stat)
}

rmna <- function(M) {
  
  M[is.na(M)] = 0
  return(M)
  
}

rmdiag <- function(M) {
  
  M[as.logical(diag(nrow(M)))]=0
  return(M)
}