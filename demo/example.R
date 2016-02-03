setwd("~/workspace/mdm") # set path to your working directory
source("R/mdm.R")

load("data/dts.Rda")
ts.plot(dts, type="l", col=c("red", "blue"))
cor(dts)^2

# add third node
rn1 = rnorm(230, 0, 0.2)
rn2 = rnorm(230, 0, 0.2)
rn3 = rnorm(230, 0, 0.2)
rn4 = rnorm(230, 0, 0.2)

d = cbind(dts,rn1)
d = cbind(dts,rn1,rn2,rn3,rn4)
cor(d)^2 > .3 & cor(d)^2 < 1

idx = 1:30
ts.plot(d[idx,1:3], type="l", col=c("red", "blue", 'black'))

dim(d)
mymodel = exhaustive.search(d,2)

# for burster

run=as.numeric(Sys.getenv('SGE_TASK_ID')) # Use the task ID to run the code in parallel

load("Data_NS.Rda") # A list object, as the number of time points differs between subjects
Data=Data_NS

Ns=dim(Data)[1]
Nn=dim(Data)[3]

subj=as.vector(t(matrix(rep(1:Ns,Nn),Ns,Nn)))
nodes=rep(1:Nn,Ns)

model.set=exhaustive.search(Data=Data[subj[run],,],node=nodes[run],nbf=15,delta=seq(0.5,1,0.01))
write.table(model.set$model.store,file=paste("NS_redone_subj",subj[run],"_n_",nodes[run],".txt",sep=""))
write.table(model.set$runtime[3],file=paste("NS_redone_subj",subj[run],"_n_",nodes[run],"_run_time.txt",sep=""))
