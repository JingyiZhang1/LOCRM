######################################################################################
#  Local Continual Reassessment Methods for Dose Finding in Drug-Combination Trials  #
#         Jingyi Zhang, Fangrong Yan, Nolan A. Wages∗ and Ruitao Lin*                #
#              Nolan.Wages@vcuhealth.org and rlin@mdanderson.org                     #
######################################################################################
# -----------------------------------------------------------------------------------#
#        Main input variables                                                        #
#        p.true.tox    --> True toxicity rate of drug combination                    #
#        target.tox    --> target toxicity rate                                      #
#        cutoff.tox    --> cutoff probability for overly-toxic combination           #
#        ntrial        --> number of replications                                    #
#        namx          --> maximum sample size                                       #
#        cohortsize    --> number of patients in one cohort                          #
# -----------------------------------------------------------------------------------#

locrm <- function(p.true.tox,target.tox,cutoff.tox,ntrial,nmax,cohortsize){
  library(pocrm)
  library(parallel)
  library(rjags)
  ncohort <- nmax/cohortsize
  ndrugA <- dim(p.true.tox)[2]; ndrugB <- dim(p.true.tox)[1]
  r <- as.vector(t(p.true.tox))
  #Specify the skeleton valuesses
  skeleton <- getprior(0.05,0.3,3,5)
  skeleton3 <- getprior(0.05,0.3,2,3)
  skeleton4 <- getprior(0.05,0.3,3,4)
  
  comb.select <- matrix(rep(0,length(p.true.tox)), nrow=ndrugA)
  TOX <- PTS <- array(0,c(ndrugA,ncol=ndrugB,ntrial))
  
  modelstring <- "
        model {
          
          for(i in 1:ndose){ 
            yy[i] ~ dbin(pt[i],nn[i])
            pt[i] <- skeleton[i]^exp(a)
          }
          a ~ dnorm(0,1/var.a)
          
        }
        "
  
  locationC1 <- NULL
  for(id in 1:ndrugB){locationC1 <- c(locationC1,rep(id,ndrugA))}
  locationC <- cbind(locationC1, rep(1:ndrugA,ndrugB))
  locationM <- matrix(1:length(p.true.tox),nrow=ndrugB, byrow=T)
  
  for(trial in 1:ntrial){
    y=n=matrix(0,nrow=ndrugA,ncol=ndrugB)
    comb.curr = 1;  # current dose level	 
    for(icohort in 1:ncohort){
      # monitoring dose (1,1)
      if(pbeta(target.tox,1+y[1],1+n[1]-y[1],lower.tail = FALSE)>cutoff.tox){break}
      # data generation
      y[comb.curr] = y[comb.curr] + sum(rbinom(cohortsize,1,r[comb.curr]))
      n[comb.curr] = n[comb.curr] + cohortsize
      
      orders<-matrix(nrow=4,ncol=5)
      d.curr <- locationC[comb.curr,]
      adj <- c(0,0,comb.curr,0,0)
      if(d.curr[1]-1>0){adj[1] <- locationM[d.curr[1]-1,d.curr[2]]}
      if(d.curr[2]-1>0){adj[2] <- locationM[d.curr[1],d.curr[2]-1]}
      if(d.curr[1]+1<=ndrugB){adj[5] <- locationM[d.curr[1]+1,d.curr[2]]}
      if(d.curr[2]+1<=ndrugA){adj[4] <- locationM[d.curr[1],d.curr[2]+1]}
      orders[1,]<-c(adj[1],adj[2],adj[3],adj[4],adj[5])
      orders[2,]<-c(adj[2],adj[1],adj[3],adj[4],adj[5])
      orders[3,]<-c(adj[1],adj[2],adj[3],adj[5],adj[4])
      orders[4,]<-c(adj[2],adj[1],adj[3],adj[5],adj[4])
      if((d.curr[1]==1 & d.curr[2]==ndrugA) | (d.curr[1]==ndrugB & d.curr[2]==1)){
        orders <- orders[1,]
      }else{
        if(d.curr[1]==1 | d.curr[2]==1 ){ orders <- orders[c(1,3),] }
        if(d.curr[1]==ndrugB | d.curr[2]==ndrugA ){ orders <- orders[c(1,2),] }
      }
      if(is.vector(orders)){
        # if a single ordering is inputed as a vector, convert it to a matrix
        orders <- t(as.matrix(orders))
      }
      ske.use <- NULL
      if(sum(orders[1,]>0)==3){ ske.use <- skeleton3 }
      if(sum(orders[1,]>0)==4){ ske.use <- skeleton4 }
      if(sum(orders[1,]>0)==5){ ske.use <- skeleton }
      orders.re <- matrix(orders[which(orders>0)],nrow=nrow(orders))
      ske <- getwm(orders.re,ske.use)
      if(is.vector(ske)){ske <- t(as.matrix(ske))}
      orders.nz <- ptox.hat <- lik <- NULL
      # model average
      for(ip in 1:nrow(orders)){
        jags.data <- list(ndose=length(orders.re[1,]), # number of combos
                          yy=y[orders.re[1,]],# number of toxicities
                          nn=n[orders.re[1,]], # number of patients
                          skeleton=ske[ip,],
                          var.a=2)
        jags <- jags.model(textConnection(modelstring),data =jags.data,n.chains=1,n.adapt=5000,quiet=TRUE)
        t.sample <- coda.samples(jags,c('pt'),n.iter=2000,progress.bar="none")
        ptox.hat <- rbind(ptox.hat, colMeans(as.matrix(t.sample)))
        
        ndose_bma=length(orders.re[1,])
        yy_bma=y[orders.re[1,]]
        nn_bma=n[orders.re[1,]]
        skeleton_bma=ske[ip,]
        a1 <- rnorm(10000,0,sqrt(2))
        loglik_dose <- matrix(0,nrow=10000,ncol=ndose_bma)
        for(i in 1:ndose_bma){
          loglik_dose[,i] <- yy_bma[i]*exp(a1)*log(skeleton_bma[i])+(nn_bma[i]-yy_bma[i])*log((1-skeleton_bma[i]^exp(a1))) }
        lik <- c(lik,mean(rowSums(exp(loglik_dose))))
      }
      mprior.tox = rep(1/nrow(orders),nrow(orders));  # prior for each toxicity ordering
      postprob.tox = (lik*mprior.tox)/sum(lik*mprior.tox);
      pp.tox <- matrix(rep(postprob.tox,dim(ptox.hat)[2]), nrow=dim(ptox.hat)[1])
      adj.nz <- adj[which(adj>0)]
      loss=abs(colSums(pp.tox*ptox.hat)-target.tox)
      # dose assignment
      comb.curr <- sort(orders.re[1,])[which(loss==min(loss))]
      if(length(which(loss==min(loss)))>1){
        comb.curr <- orders[sample(which(loss==min(loss)),size=1)] }
    }
    if(sum(n) == nmax){
      # Isotonic regression
      phat = (y + 0.05)/(n + 0.1)
      phat = Iso::biviso(phat, n + 0.1, warn = TRUE)[,]
      phat = phat + (1e-05) * (matrix(rep(1:dim(n)[1], each = dim(n)[2], len = length(n)), dim(n)[1], byrow = T) + 
                                 matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)), dim(n)[1]))
      phat[n == 0] = 10
      comb.curr2 = which(abs(phat - target.tox) == min(abs(phat - target.tox)))
      comb.select[comb.curr2]=comb.select[comb.curr2]+1;
    }
    TOX[,,trial]=y
    PTS[,,trial]=n
  }
  sel <- t(comb.select)*100/ntrial
  tox <- t(round(apply(TOX,c(1,2),mean),2));ntox=sum(tox)
  pts <- t(round(apply(PTS,c(1,2),mean),2));npts=sum(pts)
  
  dlt.rate <- ntox*100/npts
  p.stop <- 100-sum(sel)
  
  result.list <- list(sel = sel, tox = tox, pts = pts, ntox = ntox, npts = npts, p.stop = round(p.stop,2),
                      dlt.rate = dlt.rate );result.list
  return(result.list)
}

## example ##
# p.true.tox <- rbind(c(0.15,0.30,0.45,0.50,0.60),
#                     c(0.30,0.45,0.50,0.60,0.75),
#                     c(0.45,0.55,0.60,0.70,0.80))
# locrm(p.true.tox,target.tox=0.3,cutoff.tox=0.95,ntrial=100,nmax=51,cohortsize=3)

