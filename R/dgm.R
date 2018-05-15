#' Specify the priors. Without inputs, defaults will be used.
#' 
#' @param m0 the value of the prior mean at time \code{t=0}, scalar (assumed to be the same
#'  for all nodes). The default is zero.
#' @param CS0 controls the scaling of the prior variance matrix \code{C*_{0}} at time 
#'  \code{t=0}. The default is 3, giving a non-informative prior for \code{C*_{0}, 3 x (p x p)}
#'  identity matrix. \code{p} is the number of thetas.
#' @param n0 prior hyperparameter of precision \code{phi ~ G(n_{0}/2; d_{0}/2)}. The default
#'  is a non-informative prior, with \code{n0 = d0 = 0.001}. \code{n0} has to be higher than 0.
#' @param d0 prior hyperparameter of precision \code{phi ~ G(n_{0}/2; d_{0}/2)}. The default
#'  is a non-informative prior, with \code{n0 = d0 = 0.001}.
#' 
#' @details At time \code{t=0}, \code{(theta_{0} | D_{0}, phi) ~ N(m_{0},C*_{0} x phi^{-1})},
#'  where \code{D_{0}} denotes the set of initial information.
#' 
#' @return \code{priors} a list with the prior hyperparameters. Relevant to \code{\link{dlm.lpl},
#'  \link{exhaustive.search}, \link{node}, \link{subject}}.
#'  
#' @references West, M. & Harrison, J., 1997. Bayesian Forecasting and Dynamic Models. Springer New York.
#' 
#' @examples
#' pr=priors.spec()
#' pr=priors.spec(n0=0.002)
priors.spec <- function(m0 = 0, CS0 = 3, n0 = 0.001, d0 = 0.001) {
  
  priors = list(m0 = m0, CS0 = CS0, n0 = n0, d0 = d0)
  return(priors)
  
}

#' Calculate the log predictive likelihood for a specified set of parents and a fixed delta.
#'
#' @param Yt the vector of observed time series, length \code{T}.
#' @param Ft the matrix of covariates, dim = number of thetas (\code{p}) x number of time
#'  points (\code{T}), usually a row of 1s to represent an intercept and the time series of
#'  the parent nodes.
#' @param delta discount factor (scalar).
#' @param priors list with prior hyperparameters.
#' 
#' @return
#' \item{mt}{the vector or matrix of the posterior mean (location parameter), dim = \code{p x T}.}
#' \item{Ct}{and \code{CSt} the posterior scale matrix \code{C_{t}} is \code{C_{t} = C*_{t} x S_{t}},
#'  with dim = \code{p x p x T}, where \code{S_{t}} is a point estimate for the observation variance
#'  \code{phi^{-1}}}
#' \item{Rt}{and \code{RSt} the prior scale matrix \code{R_{t}} is \code{R_{t} = R*_{t} x S_{t-1}},
#'  with dim = \code{p x p x T}, where \code{S_{t-1}} is a point estimate for the observation
#'  variance \code{phi^{-1}} at the previous time point.}
#' \item{nt}{and \code{dt} the vectors of the updated hyperparameters for the precision \code{phi}
#'  with length \code{T}.}
#' \item{S}{the vector of the point estimate for the observation variance \code{phi^{-1}} with
#'  length \code{T}.}
#' \item{ft}{the vector of the one-step forecast location parameter with length \code{T}.}
#' \item{Qt}{the vector of the one-step forecast scale parameter with length \code{T}.}
#' \item{ets}{the vector of the standardised forecast residuals with length \code{T},
#'  \eqn{\newline} defined as \code{(Y_{t} - f_{t}) / sqrt (Q_{t})}.}
#' \item{lpl}{the vector of the Log Predictive Likelihood with length \code{T}.}
#' 
#' @references West, M. & Harrison, J., 1997. Bayesian Forecasting and Dynamic Models. Springer New York.
#' 
#' @examples
#' data("utestdata")
#' Yt = myts[,1]
#' Ft = t(cbind(1,myts[,2:5]))
#' m = dlm.lpl(Yt, Ft, 0.7)
#' 
#' 
dlm.lpl <- function(Yt, Ft, delta, priors = priors.spec() ) {
  
  m0 = priors$m0
  CS0 = priors$CS0
  n0 = priors$n0
  d0 = priors$d0
  
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
  
  for (i in 2:Nt){
    
    # Posterior at {t-1}: (theta_{t-1}|D_{t-1}) ~ T_{n_{t-1}}[m_{t-1}, C_{t-1} = C*_{t-1} x d_{t-1}/n_{t-1}]
    # Prior at {t}: (theta_{t}|D_{t-1}) ~ T_{n_{t-1}}[m_{t-1}, R_{t}]
    # D_{t-1} = D_{0},Y_{1},...,Y_{t-1} D_{0} is the initial information set
    
    # R*_{t} = C*_{t-1}/delta
    RSt[,,i] = CSt[,,(i-1)] / delta
    Rt[,,i] = RSt[,,i] * S[(i-1)] 
    # One-step forecast: (Y_{t}|D_{t-1}) ~ T_{n_{t-1}}[f_{t}, Q_{t}]
    ft[i] = t(F1[,i]) %*% mt[,(i-1)]
    QSt = as.vector(1 + t(F1[,i]) %*% RSt[,,i] %*% F1[,i])
    Qt[i] = QSt * S[(i-1)]
    et = Y[i] - ft[i]
    ets[i] = et / sqrt(Qt[i])
    
    # Posterior at t: (theta_{t}|D_{t}) ~ T_{n_{t}}[m_{t}, C_{t}]
    # D_{t} = D_{0},Y_{1},...,Y_{t}
    At = (RSt[,,i] %*% F1[,i])/QSt
    mt[,i] = mt[,(i-1)] + (At*et)
    
    nt[i] = nt[(i-1)] + 1
    dt[i] = dt[(i-1)] + (et^2)/QSt
    S[i]=dt[i]/nt[i] 
    
    CSt[,,i] = RSt[,,i] - (At %*% t(At))*QSt
    Ct[,,i] = S[i]*CSt[,,i]
    
    # Log Predictive Likelihood 
    lpl[i] = lgamma((nt[(i-1)]+1)/2)-lgamma(nt[(i-1)]/2)-0.5*log(pi*nt[(i-1)]*Qt[i])-((nt[(i-1)]+1)/2)*log(1+(1/nt[(i-1)])*et^2/Qt[i])
  }
  
  mt = mt[,2:Nt]; Ct = Ct[,,2:Nt]; CSt = CSt[,,2:Nt]; Rt = Rt[,,2:Nt]; RSt = RSt[,,2:Nt]
  nt = nt[2:Nt]; dt = dt[2:Nt]; S = S[2:Nt]; ft = ft[2:Nt]; Qt = Qt[2:Nt]; ets = ets[2:Nt]; lpl = lpl[2:Nt]
  
  output <- list(mt=mt,Ct=Ct,CSt=CSt,Rt=Rt,RSt=RSt,nt=nt,dt=dt,S=S,ft=ft,Qt=Qt,ets=ets,lpl=lpl)
  return(output)
}

#' A function to generate all the possible models. 
#'
#' @param Nn number of nodes; the number of columns of the dataset can be used.
#' @param node The node to find parents for.
#'
#' @return
#' output.model = a matrix with dimensions (Nn-1) x number of models, where number of models = 2^(Nn-1).
#' 
#' @examples 
#' m=model.generator(5,1)
#' 
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
#' @param node  The node to find parents for.
#' @param nbf   Log Predictive Likelihood will sum from (and including) this time point. 
#' @param delta a vector of potential values for the discount factor.
#' @param cpp boolean true (default): fast C++ implementation, false: native R code.
#' @param priors list with prior hyperparameters.
#'
#' @return
#' model.store a matrix with the model, LPL and chosen discount factor for all possible models.
#' runtime an estimate of the run time of the function, using proc.time().
#' 
#' @examples
#' data("utestdata")
#' result=exhaustive.search(myts,3)
#' 
exhaustive.search <- function(Data, node, nbf=15, delta=seq(0.5,1,0.01), cpp=TRUE, priors=priors.spec() ) {
  
  ptm=proc.time()
  
  m0 = priors$m0
  CS0 = priors$CS0
  n0 = priors$n0
  d0 = priors$d0
  
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
        a=dlm.lpl(Yt, t(Ft), delta=delta[j], priors=priors)
        lpldet[z,j]=sum(a$lpl[nbf:Nt])
      }
    }
    
    if (sum(is.na(lpldet[z,])) == length(lpldet[z,])) {
      lplmax[z] = -.Machine$double.xmax
    } else {
      lplmax[z]=max(lpldet[z,],na.rm=TRUE)
      DF.hat[z] = delta[which.max(lpldet[z,])]
    }
  }
  
  # Output model.store
  model.store=rbind(models,lplmax,DF.hat)
  rownames(model.store)=NULL
  
  runtime=(proc.time()-ptm)
  
  return(list(model.store=model.store,runtime=runtime))
}

#' Mean centers timeseries in a 2D array timeseries x nodes,
#' i.e. each timeseries of each node has mean of zero.
#'
#' @param X 2D array with dimensions timeseries x nodes.
#'
#' @return M 2D array.
#' 
#' @examples
#' data("utestdata")
#' myts=center(myts)
#' 
#' 
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
#' @param priors list with prior hyperparameters.
#' @param path a path where results are written.
#' @param method ether exhaustive, foward, backward, or both.
#'
#' @return store list with results.
#' 
#' @examples
#' data("utestdata")
#' # select only 3-nodes to speed-up this example
#' sub=subject(myts[,1:3]) 
#' sub=subject(myts[,1:3], method="both")
#' 
subject <- function(X, id=NULL, nbf=15, delta=seq(0.5,1,0.01), cpp=TRUE,
                    priors = priors.spec(), path = getwd(), method = "exhaustive") {
  
  if (!is.element(method, c("both", "exhaustive", "forward", "backward"))) {
    stop("Method must be either exhaustive, forward, backward or both")
  }
  
  Nn=ncol(X)  # nodes
  models=list()
  
  for (n in 1:Nn) {
    
    if (method == "both") {
      tmp=stepwise.combine(X, n, nbf=nbf, delta=delta, priors=priors)
      models[[n]] = tmp$model.store
    } else if (method == "forward") {
      tmp=stepwise.forward(X, n, nbf=nbf, delta=delta, priors=priors)
      models[[n]] = tmp$model.store
    } else if (method == "backward") {
      tmp=stepwise.backward(X, n, nbf=nbf, delta=delta, priors=priors)
      models[[n]] = tmp$model.store
    } else if (method == "exhaustive") {
      tmp=exhaustive.search(X, n, nbf=nbf, delta=delta, cpp=cpp, priors=priors)
      models[[n]] = tmp$model.store
    }
    
    if (!is.null(id)) {
      write(t(models[[n]]), file=file.path(path, sprintf("%s_node_%03d.txt", id, n)), 
            ncolumns = ncol(tmp$model.store))
    }
  }
  
  store=list()
  store$models=models
  store$winner=getWinner(models,Nn)
  store$adj=getAdjacency(store$winner,Nn)
  
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
#' @param priors list with prior hyperparameters.
#' @param path a path where results are written.
#' @param method can be exhaustive (default), forward, backward, or both.
#' 
#' @return store list with results.
#' 
#' @examples
#' \donttest{
#' data("utestdata")
#' m=node(myts, 3, id="SUB001_5nodes")
#' }
#' 
node <- function(X, n, id=NULL, nbf=15, delta=seq(0.5,1,0.01), cpp=TRUE, priors=priors.spec(),
                 path=getwd(), method = "exhaustive") {
  
  if (!is.element(method, c("both", "exhaustive", "forward", "backward"))) {
    stop("Method must be either exhaustive, forward, backward or both")
  }
  
  if (method == "both") {
    store=stepwise.combine(X, n, nbf=nbf, delta=delta, priors=priors)
  } else if (method == "forward") {
    store=stepwise.forward(X, n, nbf=nbf, delta=delta, priors=priors)
  } else if (method == "backward") {
    store=stepwise.backward(X, n, nbf=nbf, delta=delta, priors=priors)
  } else if (method == "exhaustive") {
    store=exhaustive.search(X, n, nbf=nbf, delta=delta, cpp=cpp, priors=priors)
  }
  
  if (!is.null(id)) {
    write(t(store$model.store), file=file.path(path, sprintf("%s_node_%03d.txt", id, n)), 
          ncolumns = ncol(store$model.store))
  }
  return(store)
}

#' Reads single subject's network from txt files.
#' @param path path.
#' @param id identifier to select all subjects' nodes, e.g. pattern containing subject ID and session number.
#' @param nodes number of nodes.
#' @param modelStore can be set to false to save memory.
#'
#' @return store list with results.
#' 
#' @examples
#' \donttest{
#' read.subject(path='~/myData', id='ID00012', nodes=5)
#' }
#' 
#' 
read.subject <- function(path, id, nodes, modelStore=TRUE) {
  
  models = list()
  for (n in 1:nodes) {
    #file=sprintf("%s_node_%03d.txt", id, n)
    #models[,,n] = as.matrix(read.table(file)) # quite slow
    file=list.files(path, pattern=glob2rx(sprintf("%s*_node_%03d.txt", id, n)))
    # we have to use a list as model dimension is variable in stepwise.
    models[[n]] = as.matrix(fread(file.path(path,file))) # faster, from package "data.table"
  }
  store=list()
  if (modelStore) {
    store$models=models
  }
  store$winner=getWinner(models,nodes)
  store$adj=getAdjacency(store$winner,nodes)
  
  return(store)
}

#' Get winner network by maximazing log predictive likelihood (LPL)
#' from a set of models.
#'
#' @param models 2D matrix, or 3D models x node.
#' @param nodes number of nodes.
#'
#' @return winner array with highest scored model(s).
#' 
getWinner <- function(models, nodes) {
  
  if (is.matrix(models)) {
    winner = models[,which.max(models[nodes+1,])]
    
  } else if (is.list(models)) {
    winner = array(0, dim=c(nodes+2,nodes))
    for (n in 1:nodes) {
      winner[,n]=models[[n]][,which.max(models[[n]][nodes+1,])]
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
#' 
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
#' @param colMapLabel label for colormap.
#' @param hasColMap FALSE turns off color map, default is NULL (on).
#' @param lim vector with min and max value, data outside this range will be removed.
#' @param gradient gradient colors.
#' @param nodeLabels node labels.
#' @param axisTextSize text size of the y and x tick labels.
#' @param xAngle orientation of the x tick labels.
#' @param titleTextSize text size of the title.
#' @param barWidth width of the colorbar.
#' @param textSize width of the colorbar.
#' 
#' @examples
#' # Generate some sample binary 5-node network structures for N=20, then compute
#' # proportion at each edge
#' N=20
#' x = array(rbinom(n=5*5*N, size=1, prob=0.30), dim=c(5,5,N))
#' A = apply(x, c(1,2), mean)
#' \donttest{
#' gplotMat(A, title = "network", colMapLabel = '%', barWidth = 0.3)
#' }
#' 
gplotMat <- function(adj, title=NULL, colMapLabel=NULL, hasColMap=NULL, lim=c(0, 1),
                     gradient=c("white", "orange", "red"), nodeLabels=waiver(), axisTextSize=12,
                     xAngle=0, titleTextSize=12, barWidth = 1, textSize=12) {
  colnames(adj)=NULL
  rownames(adj)=NULL
  
  x = melt(adj)
  names(x)[1] = "Parent"
  names(x)[2] = "Child"
  
  # handle scales in case custom labeling is set
  if (is.list(nodeLabels)) {
    x_scale = scale_x_continuous()
    y_scale = scale_y_reverse()
  } else {
    x_scale = scale_x_continuous(breaks = 1:ncol(adj), labels = nodeLabels)
    y_scale = scale_y_reverse(breaks = 1:ncol(adj), labels = nodeLabels)
  }
  
  if (is.null(hasColMap)) {
    myguides = guides(fill=guide_colorbar(barwidth = barWidth))
  } else {
    myguides = guides(fill=hasColMap)
  }
  
  ggplot(x, aes_string(x = "Child", y = "Parent", fill = "value")) +
    geom_tile(color = "gray60") +
    
    scale_fill_gradient2(
      na.value = "transparent",
      low  = gradient[1],
      mid  = gradient[2],
      high = gradient[3],
      midpoint = sum(lim)/2,
      limit = lim,
      space = "Lab",
      name = colMapLabel) + 
    
    theme(#axis.ticks.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      text = element_text(size=textSize),
      plot.title = element_text(size=titleTextSize),
      axis.text.x = element_text(size=axisTextSize,angle=xAngle),
      axis.text.y = element_text(size=axisTextSize),
      panel.background = element_blank()
      #panel.grid.major = element_line(colour="black", size = (1.5)),
      #panel.grid.minor = element_line(size = (0.2), colour="grey")
    ) +
    
    x_scale + y_scale + ggtitle(title) + myguides
}

#' Performes a binomial test with FDR correction for network edge occurrence.
#'
#' @param adj adjacency matrix, nodes x nodes x subj, or nodes x nodes x runs x subj.
#' @param alter type of binomial test, "two.sided" (default), "less", or "greater"
#' @param fdr false discovery rate (FDR) control, default is 0.05.
#'
#' @return store list with results.
#'
#' @examples
#' # Generate some sample binary 5-node network structures for N=20, then perform
#' # significance testing.
#' N=20
#' x = rmdiag(array(rbinom(n=5*5*N, size=1, prob=0.10), dim=c(5,5,N)))
#' x[1,2,2:N]=1; x[2,3,seq(1,N,2)]=1 # add some consitent edges
#' A = apply(x, c(1,2), mean)
#' l = binom.nettest(x)
#'
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
      if (i !=j) {
        tmp=binom.test(adj_[i,j],N,p=p0, alternative=alter)
        p[i,j]=tmp$p.value
      }
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
#' 
#' @examples 
#' # Let's say subjects are concatenated in a 2D matrix
#' # (samples x nodes), with each having 200 samples.
#' # generate some sample data
#' N=20
#' Nn=5
#' x = array(rnorm(200*N*Nn), dim=c(200*N,Nn))
#' ts = reshapeTs(x,N,200)
#' 
reshapeTs <- function(ts, N, V) {
  NC = ncol(ts)
  M = array(NA, dim=c(V,NC,N))
  for (i in 1:N) {
    idx = ((i-1)*V+1):(V*i)
    M[,,i] = ts[idx,]
  }
  return(M)
}

#' Mean correlation of time series across subjects.
#'
#' @param ts a 3D time series time series x nodes x subjects.
#' 
#' @return M correlation matrix.
#' 
#' @examples
#' # create some sample data with 200 samples,
#' # 5 nodes, and 2 subjects
#' ts = array(rnorm(200*5*2), dim=c(200,5,2))
#' M = corTs(ts)
#' 
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

#' Extract specific parent model with assocated df and ME from complete model space.
#'
#' @param models a 2D model matrix.
#' @param parents a vector with parent nodes.
#' 
#' @return mod specific parent model.
#' 
#' @examples
#' data("utestdata")
#' r=exhaustive.search(myts,3)
#' # get model with parents 1, 2, and 4.
#' m=getModel(r$model.store,c(1,2,4))
#' 
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

#' Get model number from a set of parents.
#'
#' @param models a 2D model matrix.
#' @param parents a vector with parent nodes.
#' 
#' @return nr model number.
getModelNr <- function(models, parents) {
  Nn = nrow(models) + 1 # No. of nodes
  Nm = ncol(models) # No. of models
  
  parents = c(parents, rep(0, Nn-length(parents)-1)) # add fill zeros
  for (m in 1:Nm) {
    if (all(models[,m] == parents)) {
      nr = m
      break;
    }
  }
  
  return(nr)
}

#' A group is a list containing restructured data from subejcts for easier group analysis.
#'
#' @param subj a list of subjects.
#' 
#' @return group a list.
#' 
#' @examples
#' # create some sample data with 200 samples,
#' # 3 nodes, and 2 subjects
#' ts = array(rnorm(200*3*2), dim=c(200,3,2))
#' mysubs=list()
#' mysubs[[1]]=subject(ts[,,1])
#' mysubs[[2]]=subject(ts[,,2])
#' g=dgm.group(mysubs)
#' 
dgm.group <- function(subj) {
  Nn=ncol(subj[[1]]$adj$am)
  N=length(subj)
  
  am = lpl = df = tam = tbi = array(NA, dim=c(Nn,Nn,N))
  df_ = array(NA, dim=c(N,Nn))
  tlpls = array(NA, dim=c(Nn,Nn,2,N))
  winner = array(NA, dim=c(Nn+2,Nn,N))
  
  for (s in 1:N) {
    am[,,s]  = subj[[s]]$adj$am
    lpl[,,s] = subj[[s]]$adj$lpl
    df[,,s]  = subj[[s]]$adj$df
    df_[s,]  = subj[[s]]$winner[nrow(subj[[s]]$winner),]
    winner[,,s]  = subj[[s]]$winner
    
    # pruning
    if (!is.null(subj[[s]]$thr)) {
      tam[,,s]  = subj[[s]]$thr$am
      tbi[,,s]  = subj[[s]]$thr$bi
      tlpls[,,,s]= subj[[s]]$thr$lpls
    }
  }
  
  group=list(am=am,lpl=lpl,df=df,tam=tam,tbi=tbi,tlpls=tlpls,
             df_=df_,winner=winner)
  
  return(group)
}

#' A group is a list containing restructured data from subejcts for easier group analysis.
#'
#' @param subj a list of subjects.
#' 
#' @return group a list.
#' 
#' @examples
#' # create some sample data with 200 samples,
#' # 3 nodes, and 2 subjects
#' ts = array(rnorm(200*3*2), dim=c(200,3,2))
#' mysubs=list()
#' mysubs[[1]]=patel(ts[,,1])
#' mysubs[[2]]=patel(ts[,,2])
#' g=patel.group(mysubs)
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

#' Get pruned adjacency network.
#'
#' @param adj list with network adjacency from getAdjacency().
#' @param models list of models.
#' @param winner matrix 2D with winning models.
#' @param e bayes factor for network pruning.
#' 
#' @return thr list with pruned network adjacency.
#'
#' @examples
#' data("utestdata")
#' # select only 3-nodes to speed-up this example
#' sub=subject(myts[,1:3])
#' p=pruning(sub$adj, sub$models, sub$winner)
pruning <- function(adj, models, winner, e = 20) {
  
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
        
        # uni i->j, LPL without parent j
        p = winner[,i][2:Nn]
        p = p[p != j & p!= 0] # remove node j
        lpls[i,j,2] = adj$lpl[i,j] + getModel(models[[i]], p)[Nn+1]
        
        # uni j->i, LPL without parent i
        p = winner[,j][2:Nn]
        p = p[p != i & p!= 0] # remove node i
        lpls[j,i,2] = adj$lpl[j,i] + getModel(models[[j]], p)[Nn+1]
      }
    }
  }
  
  # am matrix
  am=adj$am
  for (i in 1:Nn) {
    for (j in 1:Nn) {
      if (B[i,j] == 1) {
        if (lpls[i,j,1] - e <= max(lpls[i,j,2], lpls[j,i,2]) ) {
          # if unidirectional lpl is larger than bidirectional with a
          # bayes factor penatly, take the simpler unidirectional model.
          if (lpls[i,j,2] > lpls[j,i,2]) {
            am[i,j] = 1; am[j,i] = 0
          } else if (lpls[i,j,2] < lpls[j,i,2])  {
            am[i,j] = 0; am[j,i] = 1
          }
        }
      }
    }
  }
  
  thr=list()
  thr$bi=bi # bidirectional edges
  thr$lpls=lpls # lpls
  thr$am=am # adjacency matrix (pruned)
  
  return(thr)
}

#' Performance of estimates, such as sensitivity, specificity, and more.
#'
#' @param x estimated binary network matrix.
#' @param true, true binary network matrix.
#' 
#' @return p list with results.
#' 
#' @examples
#' trueNet=matrix(c(0,0,0,1,0,0,0,1,0),3,3)
#' am=matrix(c(0,0,0,1,0,1,0,1,0),3,3)
#' p=perf(am, trueNet)
perf <- function(x, true) {
  
  d = dim(x)
  Nn=d[1]
  if (length(d) == 3) {
    N=d[3]
  } else if (length(d) == 2) {
    x=array(x, dim=c(Nn,Nn,1))
    N=1
  }
  
  subj=array(NA,dim=c(N,8))
  cases=array(NA,dim=c(N,4))
  
  for (i in 1:N) {
    TP = sum(x[,,i] & true)
    FP = sum((x[,,i] - true) == 1)
    FN = sum((true - x[,,i]) == 1)
    TN = sum(!x[,,i] & !true) - ncol(x[,,i])
    
    cases[i,]=c(TP,FP,FN,TN)
    
    # see https://en.wikipedia.org/wiki/Sensitivity_and_specificity
    tpr = TP/(TP+FN) # 1
    spc = TN/(TN+FP) # 2
    ppv = TP/(TP+FP) # 3
    npv = TN/(TN+FN) # 4
    fpr = FP/(FP+TN) # 5
    fnr = FN/(TP+FN) # 6
    fdr = FP/(TP+FP) # 7
    acc = (TP+TN)/(TP+FP+FN+TN) # 8
    
    subj[i,]=c(tpr,spc,ppv,npv,fpr,fnr,fdr,acc)
  }
  colnames(cases) = c("TP", "FP", "FN", "TN")
  colnames(subj) = c("tpr", "spc", "ppv", "npv", "fpr", "fnr", "fdr", "acc")
  
  p = list()
  p$subj=subj
  p$cases=cases
  
  p$tpr = sum(cases[,1])/(sum(cases[,1]) + sum(cases[,3]))
  p$spc = sum(cases[,4])/(sum(cases[,4]) + sum(cases[,2]))
  p$acc = (sum(cases[,1]) + sum(cases[,4]))/(sum(cases[,1]) + sum(cases[,2]) + sum(cases[,3]) + sum(cases[,4]))
  p$ppv = sum(cases[,1])/(sum(cases[,1]) + sum(cases[,2]))
  return(p)
}

#' Scaling data. Zero centers and scales the nodes (SD=1).
#'
#' @param X time x node 2D matrix, or 3D with subjects as the 3rd dimension.
#'
#' @return S centered and scaled matrix.
#' 
#' @examples
#' # create some sample data
#' ts = array(rnorm(200*5, mean=5, sd=10), dim=c(200,5))
#' ts = scaleTs(ts)
#' 
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
#' 
#' @examples
#' # Generate some sample data
#' x=array(rnorm(200*5), dim=c(200,5))
#' p=patel(x)
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
  tau = matrix(0, nn, nn)
  inds = theta2 >= theta3
  tau[inds] = 1 - (theta1[inds] + theta3[inds])/(theta1[inds] + theta2[inds])
  tau[!inds] = (theta1[!inds] + theta2[!inds])/(theta1[!inds] + theta3[!inds]) - 1
  #tau=-tau # inverse 
  tau[as.logical(diag(nn))]=NA
  # tau(a,b) positive, a is ascendant to b (a is parent)
  
  # functional connectivity kappa [-1, 1]
  E=(theta1+theta2)*(theta1+theta3)
  max_theta1=pmin(theta1+theta2,theta1+theta3)
  min_theta1=pmax(array(0,dim=c(nn,nn)),2*theta1+theta2+theta3-1)
  inds = theta1>=E
  D = matrix(0, nn, nn)
  D[inds]=0.5+(theta1[inds]-E[inds])/(2*(max_theta1[inds]-E[inds]))
  D[!inds]=0.5-(theta1[!inds]-E[!inds])/(2*(E[!inds]-min_theta1[!inds]))
  
  kappa=(theta1-E)/(D*(max_theta1-E) + (1-D)*(E-min_theta1))
  kappa[as.logical(diag(nn))]=NA
  
  # thresholding
  tkappa = kappa
  tkappa[kappa >= TK[1] & kappa <= TK[2]] = 0 # this is a tow-sided test, kappa ranges from -1 to 1
  ttau = tau
  ttau[tau >= TT[1] & tau <= TT[2]] = 0
  
  # binary nets: we focus on positive associations
  net = (tkappa > 0)*(tau < 0)  # only kappa thresholded
  net[diag(nn)==1]=0 # remove NA
  tnet = (tkappa > 0)*(ttau < 0) # both tau and kappa thresholded
  tnet[diag(nn)==1]=0 # remove NA
  
  PT=list(kappa=kappa, tau=tau, tkappa=tkappa, ttau=ttau, net=net, tnet=tnet)
  return(PT)
}

#' Randomization test for Patel's kappa. Creates a distribution of values
#' kappa under the null hypothesis.
#'
#' @param X time x node x subjects 3D matrix.
#' @param alpha sign. level
#' @param K number of randomizations, default is 1000.
#'
#' @return stat lower and upper significance thresholds.
#' 
#' @examples
#' # create some sample data with 200 samples,
#' # 3 nodes, and 2 subjects
#' ts = array(rnorm(200*3*5), dim=c(200,3,5))
#' mysubs=list()
#' mysubs[[1]]=patel(ts[,,1])
#' mysubs[[2]]=patel(ts[,,2])
#' mysubs[[3]]=patel(ts[,,3])
#' mysubs[[4]]=patel(ts[,,4])
#' mysubs[[5]]=patel(ts[,,5])
#' g=patel.group(mysubs)
#' r=rand.test(rmdiag(g$kappa), K=100)
rand.test <- function(X, alpha=0.05, K=1000) {
  
  low = alpha/2      # two sided test
  up  = 1 - alpha/2
  
  N  = dim(X)[3] # Nr. of subjects
  Nn = dim(X)[2] # Nr. of nodes
  Nt = dim(X)[1] # Nr. of nodes
  
  ka = array(NA,dim=c(Nn,Nn,K)) # kappa null distribution
  ta = array(NA,dim=c(Nn,Nn,K)) # tau null distribution
  X_= array(NA, dim=c(Nt,Nn))
  for (k in 1:K) {
    s = sample(N,Nn, replace = F)
    # shuffle across subjects with fixed nodes (because of node variance)
    for (n in 1:Nn) {
      X_[,n]=X[,n,s[n]]
    }
    p=patel(X_)
    ka[,,k]=p$kappa
    ta[,,k]=p$tau
  }
  
  # two sided sign. test
  stat = list()
  stat$kappa = quantile(ka[!is.na(ka)], probs=c(low, up)) # two sided
  stat$tau   = quantile(ta[!is.na(ta)], probs=c(low, up)) # two sided
  
  return(stat)
}

#' Removes NAs from matrix.
#' 
#' @param M Matrix
#'
#' @return matrix with NAs removed.
#' 
#' @examples 
#' M=array(NA, dim=c(3,3))
#' M[1,2]=0.9
#' M=rmna(M)
rmna <- function(M) {
  
  M[is.na(M)] = 0
  return(M)
  
}

#' Removes diagonal of NA's from matrix.
#' 
#' @param M Matrix
#'
#' @return matrix with diagonal of 0's.
#' 
#' @examples 
#' M=array(rnorm(3*3), dim=c(3,3))
#' M[as.logical(diag(3))] = NA
#' M=rmna(M)
rmdiag <- function(M) {
  
  M[as.logical(diag(nrow(M)))]=0
  return(M)
}

#' Turns asymetric network into an symmetric network. Helper function to
#' determine the detection of a connection while ignoring directionality.
#' 
#' @param M 3D matrix nodes x nodes x subjects
#'
#' @return 3D matrix nodes x nodes x subjects
#' 
#' @examples 
#' M=array(NA, dim=c(3,3,2))
#' M[,,1]=matrix(c(0,0,0,1,0,0,0,1,0),3,3)
#' M[,,2]=matrix(c(0,0,0,1,0,0,0,0,0),3,3)
#' M_=symmetric(M)
symmetric <- function(M) {
  
  d = dim(M)
  R = M
  for (i in 1:d[3]) {
    R[,,i] = pmax(M[,,i], t(M[,,i]))
  }
  return(R)
}

#' Stepise forward non-exhaustive greedy search, calculates the optimum value of the discount factor.
#'
#' @param Data  Dataset with dimension number of time points \code{T} x number of nodes \code{Nn}.
#' @param node  The node to find parents for.
#' @param nbf   The Log Predictive Likelihood will sum from (and including) this time point. 
#' @param delta A vector of values for the discount factor.
#' @param max.break If \code{TRUE}, the code will break if adding / removing parents does not
#' improve the LPL. If \code{FALSE}, the code will continue to the zero parent / all parent model.
#' Default is \code{TRUE}.
#' @param priors List with prior hyperparameters.
#'
#' @return
#' model.store The parents, LPL and chosen discount factor for the subset of models scored using this method.
#' 
stepwise.forward <- function(Data, node, nbf=15, delta=seq(0.5,1,0.01), 
                             max.break=TRUE, priors=priors.spec()){
  ptm=proc.time()
  
  # Define the prior hyperparameters
  m0 = priors$m0
  CS0 = priors$CS0
  n0 = priors$n0
  d0 = priors$d0
  
  Nn = ncol(Data) # the number of nodes
  Nm = 2^(Nn-1)   # the number of models (per node)
  
  # Begin at model 1, the no parent (intercept only) model
  model.store = array(0,dim=c((Nn+2),1))
  
  # Find the Log Predicitive Likelihood and discount factor for the zero parent model
  Yt = Data[,node]    # the time series of the node we wish to find parents for
  Nt = length(Yt)     # the number of time points
  nd = length(delta)  # the number of deltas
  
  # Create Ft. For the zero parent model, this is simply a column of ones, representing an intercept.
  Ft=rep(1,Nt)      
  lpl.delta=rep(NA,nd) 
  
  for (j in 1:nd){
    lpl = dlmLplCpp(Yt,t(Ft),delta=delta[j],m0_=m0,CS0_=CS0,n0=n0,d0=d0)
    lpl.delta[j]=sum(lpl[nbf:Nt])}
  
  lpl.max = max(lpl.delta,na.rm=TRUE)
  w = which(lpl.delta==lpl.max)
  DF.hat = delta[w]
  
  # Store the model number, model, LPL and discount factor
  model.store[(Nn+1),] = lpl.max
  model.store[(Nn+2),] = DF.hat
  
  # The parents in the model
  pars = numeric(0)
  
  for (N in 1:(Nn-1)){
    
    # Find all the models with the correct length and containing the previously selected parents
    pars_add = c(1:Nn)[-c(node,pars)]
    
    ms.new = array(0,dim=c((Nn+2),(Nn-N)))
    if (length(pars>0)){ms.new[2:(length(pars)+1),] = pars}
    ms.new[(length(pars)+2),] = pars_add
    
    # Find the LPL and discount factor of this subset of models models
    nms=ncol(ms.new) # How many models are being considered?
    
    # Create empty arrays for the LPL scores and the deltas
    lpl.delta=array(NA,c(nms,length(delta)))
    lpl.max=rep(NA,nms)
    DF.hat=rep(NA,nms)
    
    # Now create Ft.
    for (i in 1:nms){
      pars = ms.new[(2:Nn),i] 
      pars = pars[pars!=0]
      Ft=array(1,dim=c(Nt,length(pars)+1))
      if (ncol(Ft)>1){Ft[,2:ncol(Ft)] = Data[,pars]}
      
      # Calculate the Log Predictive Likelihood for each value of delta, for the specified models
      for (j in 1:nd){
        lpl = dlmLplCpp(Yt,t(Ft),delta=delta[j],m0_=m0,CS0_=CS0,n0=n0,d0=d0)
        lpl.delta[i,j]=sum(lpl[nbf:Nt])}
      
      lpl.max[i] = max(lpl.delta[i,],na.rm=TRUE)
      w = which(lpl.delta[i,]==lpl.max[i])
      DF.hat[i] = delta[w]}
    
    ms.new[(Nn+1),] = lpl.max
    ms.new[(Nn+2),] = DF.hat
    
    # Find the highest LPL so far
    max.score = max(model.store[(Nn+1),]) 
    logBF = lpl.max - max.score
    W = which(logBF > 0)
    
    # Update model.store
    model.store = cbind(model.store,ms.new)
    
    if (length(W)==0 & max.break==TRUE){break} 
    else{W.max = which.max(logBF)
    
    # Update the model parents
    pars = ms.new[(2:Nn),W.max]
    pars = pars[pars!=0]}}
  
  model.store[1,] = c(1:ncol(model.store)) # attach a model number
  
  runtime=(proc.time()-ptm)
  
  return(list(model.store = model.store, runtime=runtime))
}

#' Stepise backward non-exhaustive greedy search, calculates the optimum value of the discount factor.
#'
#' @param Data  Dataset with dimension number of time points \code{T} x number of nodes \code{Nn}.
#' @param node  The node to find parents for.
#' @param nbf   The Log Predictive Likelihood will sum from (and including) this time point. 
#' @param delta A vector of values for the discount factor.
#' @param max.break If \code{TRUE}, the code will break if adding / removing parents does not
#' improve the LPL. If \code{FALSE}, the code will continue to the zero parent / all parent model.
#' Default is \code{TRUE}.
#' @param priors List with prior hyperparameters.
#'
#' @return
#' model.store The parents, LPL and chosen discount factor for the subset of models scored using this method.
#' 
stepwise.backward <- function(Data, node, nbf=15, delta=seq(0.5,1,0.01), 
                              max.break=TRUE, priors=priors.spec()){
  ptm=proc.time()
  
  # Define the prior hyperparameters
  m0 = priors$m0
  CS0 = priors$CS0
  n0 = priors$n0
  d0 = priors$d0
  
  Nn = ncol(Data) # the number of nodes
  Nm = 2^(Nn-1)   # the number of models (per node)
  
  # Begin at the all parent model
  model.store = array(0,dim=c((Nn+2),1)) 
  model.store[(2:Nn),1] = c(1:Nn)[-node]
  
  # Find the Log Predicitive Likelihood and discount factor for the all parent model
  Yt = Data[,node]    # the time series of the node we wish to find parents for
  Nt = length(Yt)     # the number of time points
  nd = length(delta)  # the number of deltas
  
  # The parents in the model
  pars = model.store[2:Nn,1]
  pars = pars[pars!=0]
  
  Ft = array(1,dim=c(Nt,Nn)) 
  Ft[,2:Nn] = Data[,pars]
  
  lpl.delta=rep(NA,nd) 
  
  for (j in 1:nd){
    lpl = dlmLplCpp(Yt,t(Ft),delta=delta[j],m0_=m0,CS0_=CS0,n0=n0,d0=d0)
    lpl.delta[j]=sum(lpl[nbf:Nt])}
  
  lpl.max = max(lpl.delta,na.rm=TRUE)
  w = which(lpl.delta==lpl.max)
  DF.hat = delta[w]
  
  model.store[(Nn+1),] = lpl.max
  model.store[(Nn+2),] = DF.hat
  
  for (N in 1:(Nn-1)){
    
    # Find all the models with the correct length and missing the previously removed parents
    pars_add = combn(pars,(Nn-(N+1)))
    
    ms.new = array(0,dim=c((Nn+2),(Nn-N))) 
    if (length(pars_add>0)){ms.new[2:(nrow(pars_add)+1),] = pars_add}
    
    # Find the LPL and discount factor of these models
    nms=ncol(ms.new) # How many models are being considered?
    
    # Create empty arrays for the lpl scores and the deltas
    lpl.delta=array(NA,c(nms,length(delta)))
    lpl.max=rep(NA,nms)
    DF.hat=rep(NA,nms)
    
    # Now create Ft. 
    for (i in 1:nms){
      pars=ms.new[(2:Nn),i]
      pars=pars[pars!=0]
      
      Ft=array(1,dim=c(Nt,length(pars)+1))
      if (ncol(Ft)>1){Ft[,2:ncol(Ft)]=Data[,pars]}
      
      # Calculate the log predictive likelihood for each value of delta, for the specified models
      for (j in 1:nd){
        lpl = dlmLplCpp(Yt,t(Ft),delta=delta[j],m0_=m0,CS0_=CS0,n0=n0,d0=d0)
        lpl.delta[i,j]=sum(lpl[nbf:Nt])}
      
      lpl.max[i] = max(lpl.delta[i,],na.rm=TRUE)
      w = which(lpl.delta[i,]==lpl.max[i])
      DF.hat[i] = delta[w]}
    
    ms.new[(Nn+1),] = lpl.max
    ms.new[(Nn+2),] = DF.hat
    
    max.score = max(model.store[(Nn+1),]) # Find the highest LPL calculated so far
    logBF = lpl.max - max.score
    W = which(logBF > 0)  
    
    # Update model.store
    model.store = cbind(model.store,ms.new)
    
    if (length(W)==0 & max.break==TRUE){break}
    else{W.max = which.max(logBF)
    
    # Update the model parents 
    pars = ms.new[(2:Nn),W.max]
    pars = pars[pars!=0]}}
  
  model.store[1,] = c(1:ncol(model.store)) # attach a model number
  
  runtime=(proc.time()-ptm)
  
  return(list(model.store = model.store, runtime=runtime))
}

#' Stepise combine
#'
#' @param Data  Dataset with dimension number of time points \code{T} x number of nodes \code{Nn}.
#' @param node  The node to find parents for.
#' @param nbf   The Log Predictive Likelihood will sum from (and including) this time point. 
#' @param delta A vector of values for the discount factor.
#' @param max.break If \code{TRUE}, the code will break if adding / removing parents does not
#' improve the LPL. If \code{FALSE}, the code will continue to the zero parent / all parent model.
#' Default is \code{TRUE}.
#' @param priors List with prior hyperparameters.
#'
#' @return
#' model.store The parents, LPL and chosen discount factor for the subset of models scored using this method.
#' 
stepwise.combine <- function(Data, node, nbf=15, delta=seq(0.5,1,0.01), 
                             max.break=TRUE, priors=priors.spec()) {
  
  ptm=proc.time()
  
  fw=stepwise.forward(Data, node, nbf=nbf, delta=delta, priors=priors, max.break=max.break)
  bw=stepwise.backward(Data, node, nbf=nbf, delta=delta, priors=priors, max.break=max.break)
  model.store = mergeModels(fw$model.store,bw$model.store)
  
  runtime=(proc.time()-ptm)
  
  return(list(model.store=model.store, runtime=runtime))
}

#' Calculate the location and scale parameters for the time-varying coefficients 
#' given all the observations. West, M. & Harrison, J., 1997. Bayesian Forecasting
#' and Dynamic Models. Springer New York.
#' 
#' @param mt the vector or matrix of the posterior mean (location parameter), dim = \code{p x T}, 
#' where \code{p} is the number of thetas (at any time \code{t}) and \code{T} is the number of time points
#' @param CSt the posterior scale matrix with dim = \code{p x p x T} (unscaled by the observation variance)
#' @param RSt the prior scale matrix with dim = \code{p x p x T} (unscaled by the observation variance)
#' @param nt vector of the updated hyperparameters for the precision \code{phi} with length \code{T}
#' @param dt vector of the updated hyperparameters for the precision \code{phi} with length \code{T}
#' 
#' @return
#' smt = the location parameter of the retrospective distribution with dimension \code{p x T}
#' sCt = the scale matrix of the retrospective distribution with dimension \code{p x p x T} 
#' 
dlm.retro <- function(mt, CSt, RSt, nt, dt) {
  
  # Convert vectors to matrices
  if (is.vector(mt)){
    mt = array(mt, dim=c(1,length(mt)))
    CSt = array(CSt, dim=c(1,1,length(CSt)))
    RSt = array(RSt, dim=c(1,1,length(RSt)))}
  
  p = nrow(mt) # the number of thetas at any time t
  Nt = ncol(mt) # the number of time points
  smt = array(NA, dim=c(p,Nt))
  sCSt = array(NA, dim=c(p,p,Nt))
  
  # Values at the last time point
  smt[,Nt] = mt[,Nt]
  sCSt[,,Nt] = CSt[,,Nt] 
  
  # For other time points
  for (i in (Nt-1):1){
    #inv.sR = solvecov(RSt[,,(i+1)], cmax = 1e+10)$inv # overcomes the limitation when it cannot invert the matrix package:fpc
    inv.sR = solve(RSt[,,(i+1)])
    B = CSt[,,i] %*% inv.sR
    smt[,i] = mt[, i] + B %*% (smt[,(i+1)] - mt[,i])
    sCSt[,,i] = CSt[,,i] + B %*% (sCSt[,,(i+1)] - RSt[,,(i+1)]) %*% t(B)}
  
  # Multiply by the observation variance
  sCt = sCSt * (dt[Nt] / nt[Nt])
  
  result = list(smt=smt, sCt=sCt)
  return(result)
}


#' Checks results and returns job number for incomplete nodes.
#' @param path path to results.
#' @param ids subjects ids.
#' @param Nr Number of runs.
#' @param Nn Number of nodes.
#'
#' @return jobs job numbers
#' 
getIncompleteNodes <- function(path, ids, Nr, Nn) {
  
  f=list.files(path, pattern=glob2rx('*.txt'))
  
  # get info flag
  info = strsplit(f[1], "_")[[1]][6]
  
  idx = rep(NA, length(ids)*Nr*Nn)
  c=1
  for (i in 1:length(ids))  {
    for (r in 1:Nr)  {
      for (n in 1:Nn)  {
        s = sprintf("%s_Run_%03d_Comp_%03d_%s_node_%03d.txt", ids[i], r, n, info, n)
        idx[c] = s %in% f
        c=c+1
      }
    }
  }
  jobs=which(!idx)
  return(jobs)
}

#' Merges forward and backward model store.
#' @param fw forward model.
#' @param bw backward model.
#'
#' @return m model store.
#' 
mergeModels <- function(fw, bw) {
  Nn = nrow(fw)-2
  
  a = apply(fw[2:(Nn),], 2, sort)
  b = apply(bw[2:(Nn),], 2, sort)
  
  idx = array(NA, ncol(b))
  for (i in 1:ncol(b)) {
    idx[i] = any(apply(a, 2, function(x, want) isTRUE(all.equal(x, want)), b[,i]))
  }
  
  model.store=cbind(fw, bw[,!idx])
  return(model.store)
}

#' Removes reciprocal connections in the lower diagnoal of the network matrix.
#' @param M adjacency matrix
#'
#' @return M adjacency matrix without reciprocal connections.
#' 
rmRecipLow <- function(M) {
  # get edges, just upper diagonal
  edges = which((M*upper.tri(M))==1, arr.ind = T)
  
  for (i in 1:nrow(edges)) {
    # if edge symmatric
    if (M[edges[i,1], edges[i,2]] == 1 &&
        M[edges[i,2], edges[i,1]] == 1) {
      M[edges[i,2], edges[i,1]] = 0 # remove edge in lower diag.
    }
  }
  
  return(M)
}

#' Quick diagnostics on delta.
#' @param path path to results files.
#' @param id subject identifier.
#' @param nodes number of nodes.

#' @return x array node model's delta
#' 
diag.delta <- function(path, id, nodes) {
  
  s = read.subject(path=path, id=id, nodes=nodes)
  x=t(as.matrix(s$winner[nodes+2,], 1, nodes))
  rownames(x) = 'delta'
  colnames(x) = seq(1,nodes)
  
  return(x)
}

#' Threshold correlation matrix to match a given number of edges.
#' @param R correlation matrix.
#' @param n number of edges.

#' @return A thresholded matrix.
#' 
cor2adj <- function(R, n) {
  
  # if uneven change to even
  if (n %% 2) {
    n = n + sample(c(-1,1),1)
  }

  l = length(R)
  v=quantile(abs(R), probs = 1-n/l)
  A = R
  A[R < v & R > -v] = 0
  
  return(A)
}

#' Comparing two population proportions on the network with FDR correction.
#' @param x1 network matrix with successes in group 1.
#' @param n1 sample size group 1.
#' @param x2 network matrix with successes in group 2.
#' @param n2 sample size group 2.
#' @param alpha alpha level for uncorrected test.
#' @param fdr alpha level for FDR.
#' 
#' @return store List with test statistics and p-values.
#' 
prop.nettest <- function(x1, n1, x2, n2, alpha=0.05, fdr=0.05) {
  
  # # Verified test with example from
  # # https://onlinecourses.science.psu.edu/stat500/node/55
  # x1 =  52; n1 =  69
  # x2 = 120; n2 = 131

  p1 = x1/n1
  p2 = x2/n2
  
  p = (x1+x2)/(n1+n2)
  z = (p1 - p2)/sqrt(p*(1-p)*(1/n1 + 1/n2))
  pval = 2*pnorm(-abs(z)) # get two-sided p-value
  
  # FDR
  pval_fdr=array(p.adjust(pval, method = "fdr"), dim=dim(pval))
  
  z_fdr=z
  z_fdr[pval_fdr>=fdr]=NA
  
  z_uncorr=z
  z_uncorr[pval>=alpha]=NA
  
  store=list()
  store$z=z
  store$pval=pval
  store$pval_fdr=pval_fdr
  store$z_fdr=z_fdr
  store$z_uncorr=z_uncorr
  
  return(store)
}


#' Comparing connectivity strenght of two groups with FDR correction.
#' @param s matrix with Nn x Nn x N.
#' @param g group assignment, vector of type factor of size N.
#' @param fdr FDR alpha level.
#' @param alpha alpha level for uncorrected test.
#' 
#' @return store List with test statistics and p-values.
#' 
ttest.nettest <- function(s, g, alpha=0.05, fdr=0.05) {
  
  Nn = dim(s)[1]
  N  = dim(s)[3]
  
  t = array(NA, c(Nn, Nn))
  t.pval = array(NA, c(Nn, Nn))
  
  for (i in 1:Nn) {
    for  (j in 1:Nn){
      if (!all(is.na(s[i,j,]))) { # if not NaN
        tmp=t.test(s[i,j,] ~ g, var.equal=TRUE)
        #tmp=wilcox.test(s[i,j,] ~ g)
        t[i,j] = tmp$statistic
        t.pval[i,j] = tmp$p.value
      }
    }
  }
  t.df = tmp$parameter
  
  # FDR
  t.pval_fdr = t.pval
  t.pval_fdr[!is.na(t.pval)] = p.adjust(t.pval[!(is.na(t.pval))], method = "fdr")
  
  t_fdr=t
  t_fdr[t.pval_fdr>=fdr]=NA
  
  t_uncorr=t
  t_uncorr[t.pval>=alpha]=NA
  
  store=list()
  store$t=t
  store$t.pval=t.pval
  
  store$t_fdr=t_fdr
  store$t.pval_fdr=t.pval_fdr
  
  store$t_uncorr=t_uncorr
  
  return(store)
}