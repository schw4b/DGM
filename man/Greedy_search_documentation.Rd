\name{Greedy_search_documentation}
\alias{Greedy_search_documentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to speed up the model search}
\description{Functions for forward selection and backward elimination.}
\usage{stepwise.forward(Data,node,nbf=15,delta=seq(0.5,1,0.01),max.break=TRUE,priors=priors.spec())
stepwise.backward(Data,node,nbf=15,delta=seq(0.5,1,0.01),max.break=TRUE,priors=priors.spec())
stepwise.combine(forward_matrix,backward_matrix)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Data}{Dataset with dimension number of time points \code{T} x number of nodes \code{Nn}.}
\item{node}{The node to find parents for.}
\item{nbf}{The Log Predictive Likelihood will sum from (and including) this time point.}
\item{delta}{A vector of values for the discount factor.}
\item{max.break}{If \code{TRUE}, the code will break if adding / removing parents does not improve the LPL. If \code{FALSE}, the code will continue to the zero parent / all parent model.}
\item{priors.spec}{List with prior hyperparameters.}
\item{forward_matrix}{The winning sets of parents using a Forward Selection model search. A matrix with dimension \code{Nn+2} x \code{Nn}, rows \code{1:Nn} are the parents (ones and zeros), rows \code{(Nn+1):(Nn+2)} are the LPL and discount factor.}
\item{backward_matrix}{The winning sets of parents using a Backward Elimination model search. A matrix with dimension \code{Nn+2} x \code{Nn}, rows \code{1:Nn} are the parents (ones and zeros), rows \code{(Nn+1):(Nn+2)} are the LPL and discount factor.}
}

\details{}
\value{
\item{model.store}{The parents, LPL and chosen discount factor for the subset of models scored using this method.}
\item{step_combine_matrix}{The adjacency network, LPLs and discount factors when the Forward Selection and Backward Elimination model searches are combined.}
}

\references{}
\author{Ruth Harbord}
\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
