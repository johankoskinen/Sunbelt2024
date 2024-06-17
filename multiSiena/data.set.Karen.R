#### === functions for creating a set of dummy multilevel SAOM objects =====
data.set.karen.set.up <- function(n,M,seed=123,nbrNodes=NULL,nmain=NULL,nwarm=NULL){
require(multiSiena)
require(sna)
average.degree <- 3
net.density <- average.degree/(n-1)
net1 <- rgraph(n,tprob =net.density)
net2 <- net1
ties.list <- which( net1==1, arr.ind = TRUE)
diag(net1) <- 99
none.tie.list <-  which( net1==0, arr.ind = TRUE)
diag(net1) <- 0
net1[ ties.list[1,1], ties.list[1,2]] <- 0

net1[ none.tie.list[1,1], none.tie.list[1,2]] <- 1

behave.DV <- round(matrix( runif(n*2)*5,n,2))

mynet1 <- sienaDependent(array(c(net1, net2), dim=c(n,n,2) ),allowOnly = FALSE)# 
mybeh <- sienaDependent( behave.DV, type = "behavior" ,allowOnly = FALSE)# w
mydata <- sienaDataCreate(mynet1,mybeh)
myeffects <- getEffects( mydata )
myeffects <- setEffect( myeffects, name = "mynet1",density,initialValue = -1.4163)# 
myeffects <- setEffect( myeffects,recip,initialValue = 1.1383 )# set
myeffects <- includeEffects( myeffects, transTrip, transRecTrip )
myeffects <- setEffect( myeffects, egoX, 
                        interaction1 = "mybeh" ,initialValue = -.5 )# social selection effects
myeffects <- setEffect( myeffects,  altX, 
                        interaction1 = "mybeh" ,initialValue = .25 )# social selection effects
myeffects <- setEffect( myeffects, egoXaltX,
                        interaction1 = "mybeh",initialValue = .5 )# social selection homophilye
myeffects <- setEffect( myeffects, avAlt, name="mybeh",
                        interaction1 = "mynet1",   initialValue = 1.2 ) # a positive influence effect
myeffects <- setEffect(myeffects, name = "mynet1", Rate, type="rate",
                       initialValue =  4.2525)
myeffects <- setEffect(myeffects, name = "mybeh", Rate, type="rate",
                       initialValue = 3)

sim_model  <-  sienaAlgorithmCreate( 
  projname = 'sim_model',
  cond = FALSE, 
  useStdInits = FALSE, nsub = 0 ,
  n3 = M,# run M simulations
  simOnly = TRUE)   # by 


sim_ans <- siena07( sim_model,#simulation settings
                    data = mydata,# our starting data
                    effects = myeffects,# the effects and paramter values we have set for model
                    returnDeps = TRUE,# this is to actually return the simulated networks and behaviours
                    batch=TRUE )

simnets <- vector('list',M)

simbeh <- matrix(NA,n,M)


for (i in c(1:M)){
  
  emptyNetwork <- network::network.initialize(n, 
                                              bipartite = NULL)
  matrixNetwork <- sparseMatrixExtraction(i, sim_ans$f, sim_ans$sims, 
                                          1, 'Data1', "mynet1")
  sparseMatrixNetwork <- as(matrixNetwork, "dgTMatrix")
  simnets[[i]] <- network::network.edgelist(cbind(sparseMatrixNetwork@i + 
                                                    1, sparseMatrixNetwork@j +1, 1), emptyNetwork)
  
  simbeh[,i] <- sim_ans$sims[[i]][[1]][[2]][[1]]
}





# eval(parse(text= paste("x",k,"<- 5",sep=''))
network.list <- vector('list',M)

for (i in c(1:M))
{
  net2 <- as.matrix.network(simnets[[i]])
  G  <- sienaNet(array( c(net1, net2 ), dim = c(n, n, 2)))
  att.g <- cbind( behave.DV[,1],simbeh[,i])
  beh <- sienaDependent(att.g, type = "behavior",allowOnly=FALSE )
  eval(parse(text= paste("grp.",i,"  <- sienaDataCreate(net = G,beh = beh)",sep='') ) )
  eval(parse(text= paste("network.list[[",i,"]]"," <- grp.",i,sep='') ) )
}


my.Karen <- sienaGroupCreate(network.list)

## == define model
GroupEffects <- getEffects(my.Karen)
GroupEffects <- includeEffects( GroupEffects, transTrip)
GroupEffects <- includeEffects( GroupEffects, simX, interaction1 = "beh" )
GroupsModel <- sienaAlgorithmCreate(projname = 'Karen',seed=seed)

## === non-Bayesian estimation was used to check stability:

#model.1 <- siena07(GroupsModel, data = FourGroups, effects = GroupEffects)
#model.1

# === set random effects:
GroupEffects <- setEffect( GroupEffects, density,random=TRUE)
GroupEffects <- setEffect( GroupEffects, recip,random=TRUE)
GroupEffects <- setEffect( GroupEffects, transTrip,random=TRUE)
GroupEffects <- setEffect(GroupEffects, simX, interaction1 = "beh" ,random=TRUE)
GroupEffects <- setEffect(GroupEffects,linear,name ='beh'  ,random=TRUE)

## === set up sienaBayes ===
p <- 7 # 2x Rates plus 5 effects
Mu <- rep(0,p)
Mu[2] <- -0.2  # outdegree
Mu[3] <- 0.5   # reciprocity
Mu[5] <- 0.4 # beh similarity


Sig <- matrix(0,p,p)
diag(Sig) <- 0.01


#groupModel.e <- sienaBayes(GroupsModel, data = my.Karen,
                           # initgainGlobal=0.1, initgainGroupwise = 0.001,
                           # effects = GroupEffects, priorMu = Mu, priorSigma = Sig,
                           # priorKappa = 0.01,
                           # nwarm=nwarm, nmain=nmain, nrunMHBatches=40,
                           # nbrNodes=nbrNodes, silentstart=FALSE)
Karen <- list(GroupsModel=GroupsModel,my.Karen=my.Karen,effects = GroupEffects, priorMu = Mu, priorSigma = Sig,
              priorKappa = 0.01,
              nwarm=nwarm, nmain=nmain, nrunMHBatches=40,
              nbrNodes=nbrNodes)
Karen 
}