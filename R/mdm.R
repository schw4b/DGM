#' Calculate the filtering distribution for a specified set of parents and a fixed delta.
#'
#' @param Yt time series of the node of interest (dim = Nt).
#' @param Ft the time series of the parents (number of parents = p), and 1 for the intercept (dim = p x Nt).
#' @param delta discount factor (scalar).
#' @param m0 prior means at time t=0 with length p. The default is non-informative prior, with zero mean.
#' @param CS0 squared matrix of prior variance. The default is non-informative prior, with prior variance equal to 3 times the observed variance.
#' @param n0 prior hypermarameters of precision phi ~ G(n0/2; d0/2). The default is non-informative priors, with value of 0.001. n0 has to be higher than 0.
#' @param d0 prior hypermarameters of precision phi ~ G(n0/2; d0/2). The default is non-informative priors, with value of 0.001. n0 has to be higher than 0.
#'
#' @return
#' mt = the filtered posterior mean, dim = p x T.
#' Ct = the filtered posterior variance, dim = p x p x T.
#' CSt.
#' Rt = the prior variance, dim = p X p X T.
#' RSt
#' nt and dt = the prior hypermarameters of precision phi with length T.
#' ft = the one-step forecast mean with length T.
#' Qt = the one-step forecast variance with length T.
#' ets = the standardised residuals with length T.
#' lpl = Log Predictive Likelihood with length T.
#' 
#' @export
dlm.filt.rh <- function(Yt, Ft, delta, m0 = numeric(nrow(Ft)), CS0 = 3*diag(nrow(Ft)), n0 = 0.001, d0 = 0.001){
  
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
  Ct[,,1] = CS0
  CSt = array(0,dim=c(p,p,Nt))
  
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
  
  # Filtering
  
  for (t in 2:Nt){
    
    # Posterior at {t-1}: (theta_{t-1}|y_{t-1}) ~ T_{n_{t-1}}[m_{t-1}, C_{t-1} = C*_{t-1} x d_{t-1}/n_{t-1}]
    # Prior at {t}: (theta_{t}|y_{t-1}) ~ T_{n_{t-1}}[m_{t}, R_{t}]
    
    # RSt ~ C*_{t-1}/delta
    RSt[,,t] = Ct[,,(t-1)] / (S[(t-1)]*delta)
    Rt[,,t] = RSt[,,t] * (S[(t-1)]) 
    # One-step forecast: (Y_{t}|y_{t-1}) ~ T_{n_{t-1}}[f_{t}, Q_{t}]
    ft[t] = t(F1[,t]) %*% mt[,(t-1)]
    QSt = as.vector(1 + t(F1[,t]) %*% RSt[,,t] %*% F1[,t])
    Qt[t] = QSt * S[(t-1)]
    et = Y[t] - ft[t]
    ets[t] = et / sqrt(Qt[t])
    
    # Posterior at t: (theta_{t}|y_{t}) ~ T_{n_{t}}[m_{t}, C_{t}]
    At = (RSt[,,t] %*% F1[,t])/QSt
    mt[,t] = mt[,(t-1)] + (At*et)
    
    nt[t] = nt[(t-1)] + 1
    dt[t] = dt[(t-1)] + (et^2)/QSt
    S[t]=dt[t]/nt[t] 
    
    CSt[,,t] = RSt[,,t] - (At %*% t(At))*QSt
    Ct[,,t] = S[t]*CSt[,,t]
    
    # Log Predictive Likelihood (degrees of freedom = nt[(t-1)], not nt[t])
    lpl[t] = lgamma((nt[(t-1)]+1)/2)-lgamma(nt[(t-1)]/2)-0.5*log(pi*nt[(t-1)]*Qt[t])-((nt[(t-1)]+1)/2)*log(1+(1/nt[(t-1)])*et^2/Qt[t])
  }
  
  mt = mt[,2:Nt]; Ct = Ct[,,2:Nt]; CSt = CSt[,,2:Nt]; Rt = Rt[,,2:Nt]; RSt = RSt[,,2:Nt]
  nt = nt[2:Nt]; dt = dt[2:Nt]; S = S[2:Nt]; ft = ft[2:Nt]; Qt = Qt[2:Nt]; ets = ets[2:Nt]; lpl = lpl[2:Nt]
  
  filt.output <- list(mt=mt,Ct=Ct,CSt=CSt,Rt=Rt,RSt=RSt,nt=nt,dt=dt,S=S,ft=ft,Qt=Qt,ets=ets,lpl=lpl)
  return(filt.output)
}

#' A function to generate all the possible models. 
#'
#' @param Nn number of nodes; the number of columns of the dataset can be used.
#' @param node the node of interest (i.e., the node to find parents for).
#'
#' @return
#' output.model = a matrix with dimensions (nn-1) x number of models, where number of models = 2^(nn-1).
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
#' @param Data  Dataset with dimension number of time points Nt x Number of nodes Nn.
#' @param node  node of interest.
#' @param nbf   Log Predictive Likelihood will be calculated from this time point. 
#' @param delta a vector of potential values for the discount factor.
#' @param cpp boolean true (default): fast C++ implementation, false: native R code.
#'
#' @return
#' model.store = a matrix with the model, LPL and chosen discount factor for all possible models.
#' 
#' @export
exhaustive.search <- function(Data,node,nbf=15,delta=seq(0.5,1,0.01),cpp=TRUE) {
  
  ptm=proc.time()  
  
  Nn=ncol(Data) # the number of nodes
  Nm=2^(Nn-1)   # the number of models per node
  
  M=model.generator(Nn,node) # Generate all the possible models
  #M=cbind(M,M[,1]) # Add in a model to represent the AR alternative
  models=rbind(1:Nm,M) # Label each model with a 'model number'
  
  
  # Find the Log Predicitive Likelihood and associated discount factor for the zero parent model
  Yt=Data[,node]    # the time series of the node we wish to find parents for
  Nt=length(Yt)     # the number of time points
  nd=length(delta)  # the number of deltas
  
  # Create empty arrays for the lpl scores and the optimum deltas
  lpldet=array(NA,c(Nm,length(delta)))
  lplmax=rep(NA,Nm)
  DF.hat=rep(NA,Nm)
  
  # Now create Ft. 
  for (z in 1:Nm) {
    par=models[(2:Nn),z] # par is distinguished from pars, which are the selected parents at each stage
    par=par[par!=0]
    Ft=array(1,dim=c(Nt,length(par)+1))
    if (ncol(Ft)>1) {
      Ft[,2:ncol(Ft)]=Data[,par] # selects parents
    }  
    
    # Calculate the log predictive likelihood, for each value of delta, for the specified models
    for (j in 1:nd) {
      if (cpp) {
        # new C++ implementation
        lpl=c(dlmFiltCpp(Yt,t(Ft),delta[j]))
        lpldet[z,j]=sum(lpl[nbf:Nt])
      } else {
        # original native R
        a=dlm.filt.rh(Yt, t(Ft), delta=delta[j])
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

#' Search a subject's network; runs exhaustive search on very node.
#'
#' @param X 3D array with dimensions timeseries x nodes x subjects.
#' @param id subject ID. If set, results are saved to a txt file.
#'
#' @return store list with results.
#' @export
subject <- function(X, id=NULL) {
  N=ncol(X)  # nodes
  M=2^(N-1)  # rows/models
  models = array(rep(0,(N+2)*M*N),dim=c(N+2,M,N))
  
  for (n in 1:N) {
    tmp=exhaustive.search(X,n)
    models[,,n]=tmp$model.store
  }
  
  if (!is.null(id)) {
    for (n in 1:N) {
      write(t(models[,,n]), file=sprintf("%s_node_%03d.txt", id, n), ncolumns = M)
    }
  }
  store=list()
  store$models=models
  store$winner=getWinner(models,N)
  store$adj=getAdjacencyMatrix(store$winner,N)
  
  return(store)
}

#' Reads single subject's network from text files.
#'
#' @param id subject ID.
#' @param nodes number of nodes.
#'
#' @return store list with results.
#' @export
read.subject <- function(id, nodes) {
  
  models = array(0,dim=c(nodes+2,2^(nodes-1),nodes))
  for (n in 1:nodes) {
    file=sprintf("%s_node_%03d.txt", id, n)
    models[,,n] = as.matrix(read.table(file))
  }
  store=list()
  store$models=models
  store$winner=getWinner(models,nodes)
  store$adj=getAdjacencyMatrix(store$winner,nodes)
  
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

#' Get adjacency matrix from winning models.
#'
#' @param winner, 2D matrix.
#' @param nodes number of nodes.
#'
#' @return adj, 2D adjacency matrix.
#' @export
getAdjacencyMatrix <- function(winner, nodes) {
  
  adj = array(rep(0,nodes*nodes),dim=c(nodes,nodes))
  for (n in 1:nodes) {
    p = winner[2:nodes,n]  # parents
    p = p[p>0]
    adj[p,n] = 1
  }
  return(adj)
}

#' Plots network as graph.
#'
#' @param adj, 2D adjacency matrix.
#'
#' @export
plotNet <- function(adj) {
  plot.igraph(graph.adjacency(adj, mode="directed", weighted=T, diag=F))
}
