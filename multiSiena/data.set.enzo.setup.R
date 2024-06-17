# ==== Dataset Enzo ====== #
create.enzo <- function(nbrNodes=NULL,nmain=NULL,nwarm=NULL,seed=123,n=5,mu.deg =-0.2,mu.rec =0.5, mu.sim= 0.4 ){
#install.packages('RUnit')
#install.packages("RSienaTest", repos="http://R-Forge.R-project.org",type = "source")

require(multiSiena)
# n <- 5# number of nodes
w1 <- matrix(0,n,n)# universal empty matrix wave 1
w1[1,2] <- 1
w1[2,1] <- 1
w1[3,4] <- 1
w1[4,5] <- 1
w2 <- w1#  universal empty matrix wave 2
w2[3,4] <- 0
att <- matrix(0,n,2)#  universal empty attribute matrix waves 1 and 2
att[1:2,1] <- 1
coord <- matrix(c(1,1,
1,2,
2,2,
2,1,
.5,.5),n,2,byrow=TRUE)
# === group 1 ========
w1.g <- w1
w2.g <- w2
att.g <- att
# reciprochal tie
w2.g[5,4] <- 1
# transitive tie
w2.g[4,3] <- 1
w2.g[3,5] <- 1
# attribute
att.g[c(1,3),2] <- 1
G  <- sienaNet(array( c(w1.g, w2.g), dim = c(n, n, 2)))
beh <- sienaDependent(att.g, type = "behavior",allowOnly=FALSE )
grp.1  <- sienaDataCreate(net = G,beh = beh)
# === group 2 ========
w1.g <- w1
w2.g <- w2
att.g <- att
# reciprochal tie
w2.g[5,4] <- 1
w2.g[5,4] <- 1
# transitive tie
w2.g[4,3] <- 1
w2.g[2,4] <- 1
w2.g[2,5] <- 1
# attribute
att.g[c(1,2,4),2] <- 1
G  <- sienaNet(array( c(w1.g, w2.g), dim = c(n, n, 2)))
beh <- sienaDependent(att.g, type = "behavior" ,allowOnly=FALSE)
grp.2  <- sienaDataCreate(net = G,beh = beh)# === group 3 ======== group 3
w1.g <- w1
w2.g <- w2
att.g <- att
w1.g[2,3] <- 1
# remove
w2.g[4,5] <- 0
# reciprochal tie
w2.g[5,4] <- 1
w2.g[5,4] <- 1
# transitive tie
w2.g[4,3] <- 1
w2.g[1,3] <- 1
w2.g[2,3] <- 1
w2.g[2,4] <- 1
# attribute
att.g[c(1,2,4),2] <- 1
G  <- sienaNet(array( c(w1.g, w2.g), dim = c(n, n, 2)))
beh <- sienaDependent(att.g, type = "behavior",allowOnly=FALSE )
grp.3  <- sienaDataCreate(net = G,beh = beh)
# === group 4 ========
w1.g <- w1
w2.g <- w2
att.g <- att
# remove
w2.g[4,5] <- 0
# reciprochal tie
w2.g[3,4] <- 1
w2.g[4,3] <- 1
# transitive tie
w2.g[4,2] <- 1
w2.g[2,3] <- 1
w2.g[1,3] <- 1
w2.g[3,1] <- 1
# attribute
att.g[5,1] <- 1
att.g[c(1,2),2] <- 1
G  <- sienaNet(array( c(w1.g, w2.g), dim = c(n, n, 2)))
beh <- sienaDependent(att.g, type = "behavior",allowOnly=FALSE )
grp.4  <- sienaDataCreate(net = G,beh = beh)
# === group 5 ========
w1.g <- w1
w2.g <- w2
att.g <- att
# remove
w2.g[4,5] <- 0
w1.g[2,1] <- 0
# reciprochal tie
w2.g[3,4] <- 1
w2.g[4,3] <- 1
w2.g[3,2] <- 1
# transitive tie
w2.g[4,2] <- 1
w2.g[2,3] <- 1
w2.g[1,3] <- 1
w2.g[3,1] <- 1
# attribute
att.g[5,1] <- 1
att.g[c(1,2),2] <- 1
G  <- sienaNet(array( c(w1.g, w2.g), dim = c(n, n, 2)))
beh <- sienaDependent(att.g, type = "behavior",allowOnly=FALSE )
grp.5  <- sienaDataCreate(net = G,beh = beh)

# === group 6 ========
w1.g <- w1
w2.g <- w2
att.g <- att

w2.g[2,4] <- 1
# reciprochal tie
w2.g[5,4] <- 1
# transitive tie
w2.g[4,3] <- 1
w2.g[3,5] <- 1
w2.g[1,4] <- 1
# attribute
att.g[c(1,3),2] <- 1
G  <- sienaNet(array( c(w1.g, w2.g), dim = c(n, n, 2)))
beh <- sienaDependent(att.g, type = "behavior",allowOnly=FALSE )
grp.6  <- sienaDataCreate(net = G,beh = beh)


# === all together now ===
SixGroups <- sienaGroupCreate(list(grp.1,grp.2,grp.3,grp.4,grp.5,grp.6))


## == define model
GroupEffects <- getEffects(SixGroups)
GroupEffects <- includeEffects( GroupEffects, transTrip)
GroupEffects <- includeEffects( GroupEffects, simX, interaction1 = "beh" )
GroupsModel <- sienaAlgorithmCreate(projname = 'Enzo',seed=seed)

## === non-Bayesian estimation was used to check stability:

#model.1 <- siena07(GroupsModel, data = FourGroups, effects = GroupEffects)
#model.1

# === set random effects:
GroupEffects <- setEffect( GroupEffects, density,random=TRUE)
GroupEffects <- setEffect( GroupEffects, recip,random=TRUE)
GroupEffects <- setEffect( GroupEffects, transTrip,random=TRUE)
GroupEffects <- setEffect(GroupEffects, simX, interaction1 = "beh" ,random=TRUE)
GroupEffects <- setEffect(GroupEffects,linear,name ='beh'  ,random=TRUE)

## === run sienaBayes ===


p <- 7 # 2x Rates plus 5 effects
Mu <- rep(0,p)
Mu[2] <- mu.deg  # outdegree
Mu[3] <- mu.rec  # reciprocity
Mu[5] <- mu.sim # beh similarity


Sig <- matrix(0,p,p)
diag(Sig) <- 0.01


#groupModel.e <- sienaBayes(GroupsModel, data = FourGroups,
#				initgainGlobal=0.1, initgainGroupwise = 0.001,
#                effects = GroupEffects, priorMu = Mu, priorSigma = Sig,
#				priorKappa = 0.01,
#				nwarm=nwarm, nmain=nmain, nrunMHBatches=40,
#                nbrNodes=nbrNodes, silentstart=FALSE)

#save(list=ls(),file='Enzo.RData')
SixGroups
}


