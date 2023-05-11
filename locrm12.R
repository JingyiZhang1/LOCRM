###########################################################################################
#  Local Continual Reassessment Methods for Dose Optimization in Drug-Combination Trials  #
#           Jingyi Zhang, Fangrong Yan, Nolan A. Wages∗ and Ruitao Lin*                   #
#                Nolan.Wages@vcuhealth.org and rlin@mdanderson.org                        #
###########################################################################################
# ----------------------------------------------------------------------------------------#
#        Main input variables                                                             #
#        p.true.tox    --> True toxicity rate of drug combination                         #
#        p.true.eff    -->True efficacy rate of drug combination                          #   
#        target.tox    --> upper limit of toxicity rate                                   #
#        cutoff.tox    --> cutoff probability for overly-toxic combination                # 
#        target.eff    --> lower limit of efficacy rate                                   #
#        cutoff.eff    --> cutoff probability for futile combination                      #
#        ntrial        --> number of replications                                         #
#        namx          --> maximum sample size                                            #
#        cohortsize    --> number of patients in one cohort                               #
# ----------------------------------------------------------------------------------------#

locrm12 <- function(p.true.tox,p.true.eff,target.tox,cutoff.tox,target.eff,cutoff.eff,ntrial,nmax,cohortsize){
  library(pocrm)
  library(parallel)
  library(rjags)
  ncohort <- nmax/cohortsize
  r <- as.vector(t(p.true.tox))
  ndrugA <- dim(p.true.tox)[2]; ndrugB <- dim(p.true.tox)[1]
  #Specify the skeleton values.
  skeleton <- getprior(0.05,target.tox,3,5);skeleton
  skeleton3 <- getprior(0.05,target.tox,2,3);skeleton3
  skeleton4 <- getprior(0.05,target.tox,3,4);skeleton4
  
  SEL <- matrix(rep(0,length(p.true.tox)), nrow=ndrugA)
  TOX <- EFF <- PTS <- ELIMI <- array(0,c(ndrugA,ncol=ndrugB,ntrial))
  comb.select.e <- matrix(0,nrow=ndrugB,ncol=ndrugA)
  PPE <- PHAT <- array(0,c(ndrugB,ncol=ndrugA,ntrial))
  
  modelstring.T <- "
        model {
          for(i in 1:ndose){ 
            yy[i] ~ dbin(pt[i],nn[i])
            pt[i] <- skeleton[i]^exp(a)
          }
          
          a ~ dnorm(mean.a,1/var.a)
          
        }
        "
  modelstring.E <- "
        model {
          for(i in 1:ndrugB){
            for(j in 1:ndrugA){
              pe[i,j] <- pt(alpha + beta1 * doseA[j] + beta2 * doseB[i] + gamma1 * doseA[j] * doseA[j] + gamma2 * doseB[i] * doseB[i],mean.t,var.t,dof)
              yy[i,j] ~ dbin(pe[i,j],nn[i,j])
            }
          }
          alpha ~ dnorm(mean.alpha,1/var.alpha)
          beta1 ~ dnorm(mean.beta1,1/var.beta1)
          beta2 ~ dnorm(mean.beta2,1/var.beta2)
          gamma1 ~ dnorm(mean.gamma1,1/var.gamma1)
          gamma2 ~ dnorm(mean.gamma2,1/var.gamma2)
          dof ~ dunif(dof1,dof2)
        }
        "
  
  locationC1 <- NULL
  for(id in 1:ndrugB){locationC1 <- c(locationC1,rep(id,ndrugA))}
  locationC <- cbind(locationC1,
                     rep(1:ndrugA,ndrugB))
  locationM <- matrix(1:length(p.true.tox),nrow=ndrugB, byrow=T)
  row.doseA <- c(0.08,0.16,0.24,0.32,0.40)
  doseA <- ((row.doseA[1:ndrugA]-mean(row.doseA[1:ndrugA]))/sqrt(var(row.doseA[1:ndrugA])))
  doseB <- ((row.doseA[1:ndrugB]-mean(row.doseA[1:ndrugB]))/sqrt(var(row.doseA[1:ndrugB])))
  
  for(trial in 1:ntrial){
    y=yE=n=elimi=ptox.bma.matrix=peff.est.matrix=peff.monitor.matrix=matrix(0,nrow=ndrugA,ncol=ndrugB)
    comb.curr = 1;  # current dose level	 
    d.curr <- locationC[comb.curr,]
    
    # startup stage
    traceRunin <- d <- c(1,1)
    for(i in 1:ncohort){
      if(d[1]<ndrugB & d[2]<ndrugA){
        rand2 <- runif(1)
        if(rand2 < 0.5) {d <- c(d[1],d[2]+1)}else{d <- c(d[1]+1,d[2])}
        traceRunin <- rbind(traceRunin,d) }
      if(d[1]==ndrugB & d[2]<ndrugA){
        d <- c(d[1],d[2]+1)
        traceRunin <- rbind(traceRunin,d) }
      if(d[1]<ndrugB & d[2]==ndrugA){
        d <- c(d[1]+1,d[2])
        traceRunin <- rbind(traceRunin,d) }
      if(d[1]==ndrugB & d[2]==ndrugA){break}
    }
    
    for(ione in 1:dim(traceRunin)[1]){
      y[traceRunin[ione,2],traceRunin[ione,1]] = y[traceRunin[ione,2],traceRunin[ione,1]] + sum(rbinom(cohortsize,1,(t(p.true.tox))[traceRunin[ione,2],traceRunin[ione,1]]))
      yE[traceRunin[ione,2],traceRunin[ione,1]] = yE[traceRunin[ione,2],traceRunin[ione,1]] + sum(rbinom(cohortsize,1,(t(p.true.eff))[traceRunin[ione,2],traceRunin[ione,1]]))
      n[traceRunin[ione,2],traceRunin[ione,1]] = n[traceRunin[ione,2],traceRunin[ione,1]] + cohortsize
      comb.curr <- which(n>0)[length(which(n>0))]
      if(sum(y)>0){break}
    }
    
    # main stage
    n1 <- sum(n)
    for(itwo in 0:((nmax-n1)/cohortsize)){
      d.curr <- locationC[comb.curr,]
      if(itwo>0){
        y[comb.curr] = y[comb.curr] + sum(rbinom(cohortsize,1,r[comb.curr]))
        yE[comb.curr] = yE[comb.curr] + sum(rbinom(cohortsize,1,(t(p.true.eff))[comb.curr]))
        n[comb.curr] = n[comb.curr] + cohortsize
      }
      
      # toxicity
      orders<-matrix(nrow=4,ncol=5)
      adj <- c(0,0,comb.curr,0,0)
      if(d.curr[1]-1>0){adj[1] <- locationM[d.curr[1]-1,d.curr[2]]}
      if(d.curr[2]-1>0){adj[2] <- locationM[d.curr[1],d.curr[2]-1]}
      if(d.curr[1]+1<=ndrugB){adj[5] <- locationM[d.curr[1]+1,d.curr[2]]}
      if(d.curr[2]+1<=ndrugA){adj[4] <- locationM[d.curr[1],d.curr[2]+1]}
      
      # elimination rule
      if(pbeta(target.tox,1+y[1],1+n[1]-y[1],lower.tail = FALSE)>cutoff.tox){break}
      if(pbeta(target.eff,1+yE[comb.curr],1+n[comb.curr]-yE[comb.curr],lower.tail = TRUE)>cutoff.eff){ elimi[comb.curr] <- 1}
      if(sum(elimi[adj])==length(elimi[adj])){break}
      
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
      #Initial guesses of toxicity probabilities for each ordering.
      if(is.vector(orders)){# if a single ordering is inputed as a vector, convert it to a matrix
        orders <- t(as.matrix(orders))}
      ske.use <- NULL
      if(sum(orders[1,]>0)==3){ske.use <- skeleton3}
      if(sum(orders[1,]>0)==4){ske.use <- skeleton4}
      if(sum(orders[1,]>0)==5){ske.use <- skeleton}
      orders.re <- matrix(orders[which(orders>0)],nrow=nrow(orders))
      ske <- getwm(orders.re,ske.use)
      if(is.vector(ske)){ske <- t(as.matrix(ske))}
      
      orders.nz <- ptox.hat <- lik <- NULL
      for(ip in 1:nrow(ske)){
        jags.data <- list(ndose=length(orders.re[1,]),
                          yy=y[orders.re[1,]],
                          nn=n[orders.re[1,]], 
                          skeleton=ske[ip,],
                          mean.a=0,var.a =2) 
        jags <- jags.model(textConnection(modelstring.T),data =jags.data,n.chains=1,n.adapt=5000,quiet=TRUE)
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
      mprior.tox = rep(1/nrow(ske),nrow(ske));  # prior for each toxicity ordering
      postprob.tox = (lik*mprior.tox)/sum(lik*mprior.tox);
      # BMA
      pp.tox <- matrix(rep(postprob.tox,dim(ptox.hat)[2]), nrow=dim(ptox.hat)[1])
      ptox.bma <- colSums(pp.tox*ptox.hat)
      adj.nz <- adj[which(adj>0)]
      ptox.bma.matrix[adj.nz] <- ptox.bma
      # select MTD among candidate set
      temp.mtd <- which(abs(ptox.bma.matrix[adj.nz]-target.tox)==min(abs(ptox.bma.matrix[adj.nz]-target.tox)))
      cand.set <- adj.nz[which(ptox.bma<=ptox.bma[temp.mtd])]
      cand.set <- cand.set[which(elimi[cand.set]==0)]
      
      # efficacy
      cand.set9 <- NULL
      for(idose in 1:length(adj.nz)){ cand.set9 <- rbind(cand.set9,locationC[adj.nz[idose],]) }
      nn <- (t(n))[min(cand.set9[,1]):max(cand.set9[,1]),min(cand.set9[,2]):max(cand.set9[,2])]
      yy <- (t(yE))[min(cand.set9[,1]):max(cand.set9[,1]),min(cand.set9[,2]):max(cand.set9[,2])]
      jags.data <- list(yy=yy, # number of efficities
                        nn=nn, # number of patients
                        ndrugA=length(min(cand.set9[,2]):max(cand.set9[,2])), ndrugB=length(min(cand.set9[,1]):max(cand.set9[,1])),
                        doseA=doseA[min(cand.set9[,2]):max(cand.set9[,2])],
                        doseB=doseB[min(cand.set9[,1]):max(cand.set9[,1])],
                        mean.alpha=0,   var.alpha =1.3,
                        mean.beta1=0.8, var.beta1 =1.3,
                        mean.beta2=0.8, var.beta2 =1.3,
                        mean.gamma1=0,  var.gamma1 =1.3,
                        mean.gamma2=0,  var.gamma2 =1.3,
                        mean.t=0,       var.t =2,
                        dof1=2,         dof2=10) 
      jags <- jags.model(textConnection(modelstring.E),data =jags.data,n.chains=1,n.adapt=5000,quiet=TRUE)
      e.sample <- coda.samples(jags,c('pe'),n.iter=2000,progress.bar="none")
      peff.hat <- colMeans(as.matrix(e.sample))
      peff.bma <- matrix(peff.hat,nrow=length(min(cand.set9[,2]):max(cand.set9[,2])),byrow=T)
      peff.est.matrix[min(cand.set9[,2]):max(cand.set9[,2]),min(cand.set9[,1]):max(cand.set9[,1])] <- peff.bma
      
      # dose assignment
      if(length(cand.set)==0){break}
      candi.eff <- rbind(cand.set,n[cand.set],peff.est.matrix[cand.set],((nmax-sum(n)-n1)/(nmax-n1))^2,(peff.est.matrix[cand.set]>((nmax-sum(n))/nmax)^2))
      if(ncol(candi.eff)>1){candi.eff <- candi.eff[,order(candi.eff[3,],decreasing=T)]}
      for(ican in 1:length(cand.set)){
        if(candi.eff[2,ican]==0|all(candi.eff[2,]>0)==TRUE){
          comb.curr <- candi.eff[1,ican];break
        }else{
          if(candi.eff[5,ican]==1){
            comb.curr <- candi.eff[1,ican];break
          }else{candi.eff[1,ican] <- 0}
        }
      }
    }
    
    if(sum(n) == nmax){
      # Isotonic regression
      phat = (y + 0.05)/(n + 0.1)
      phat = Iso::biviso(phat, n + 0.1, warn = TRUE)[,]
      PHAT[,,trial]=t(phat)
      phat[elimi == 1] = 1.1
      phat = phat + (1e-05) * (matrix(rep(1:dim(n)[1], each = dim(n)[2], len = length(n)), dim(n)[1], byrow = T) + 
                                 matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)), dim(n)[1]))
      phat[n == 0] = 10
      comb.curr2 = which(abs(phat - target.tox) == min(abs(phat - target.tox)))
      if(length(comb.curr2)>1){
        comb.curr2 <- comb.curr2[sample(1:length(comb.curr2),1)] }
      mtd <- locationC[comb.curr2,]
      
      # Efficacy analysis using all data
      jags.data <- list(yy=t(yE), nn=t(n), ndrugA=ndrugA, ndrugB=ndrugB,doseA=doseA,doseB=doseB,
                        mean.alpha=0,   var.alpha =1.3,
                        mean.beta1=0.8, var.beta1 =1.3,
                        mean.beta2=0.8, var.beta2 =1.3,
                        mean.gamma1=0,  var.gamma1 =1.3,
                        mean.gamma2=0,  var.gamma2 =1.3,
                        mean.t=0,       var.t =2,
                        dof1=2,         dof2=10) 
      jags <- jags.model(textConnection(modelstring.E),data =jags.data,n.chains=1,n.adapt=5000,quiet=TRUE)
      ppe.sample <- coda.samples(jags,c('pe'),n.iter=2000,progress.bar="none")
      ppe <- matrix(colMeans(as.matrix(ppe.sample)),nrow=ndrugB)
      PPE[,,trial]=ppe
      for(i in 1:ndrugB){
        for(j in 1:ndrugA){
          if(i>mtd[1]){ppe[i,j] <- -100}
          if(j>mtd[2]){ppe[i,j] <- -100}
          if(t(elimi)[i,j]>cutoff.eff){ppe[i,j] <- -100}
        }
      }
      comb.curr3 <- which(ppe==max(ppe))
      comb.select.e[comb.curr3]=comb.select.e[comb.curr3]+1;
    }
    TOX[,,trial]=y
    EFF[,,trial]=yE
    PTS[,,trial]=n
    ELIMI[,,trial]=elimi
  }
  
  sel = comb.select.e*100/ntrial
  tox <- t(round(apply(TOX,c(1,2),mean),2));ntox=sum(tox)
  eff <- t(round(apply(EFF,c(1,2),mean),2));neff=sum(eff)
  pts <- t(round(apply(PTS,c(1,2),mean),2));npts=sum(pts)
  
  true.obd0 <- which(p.true.tox<=target.tox)
  true.obd <- which(p.true.eff[true.obd0]==max(p.true.eff[true.obd0]))
  sel.obd <- sum(sel[true.obd])
  pts.obd <- sum(pts[true.obd])
  
  true.tar <- which(p.true.eff[true.obd0]>=0.45)
  sel.tar <- sum(sel[true.tar])
  pts.tar <- sum(pts[true.tar])
  
  overtox <- which(p.true.tox>target.tox)
  sel.overtox <- sum(sel[overtox])
  pts.overtox <- sum(pts[overtox])
  
  result.list <- list(sel = sel, tox=tox, eff = eff, pts = pts, 
                      ntox = sum(y), neff = sum(yE), npts = sum(n), p.stop = 100-sum(sel) );result.list
  return(result.list)
}

# ## example ##
# p.true.tox <- matrix(c(0.05,0.15,0.30,0.45,0.55,  0.15,0.30,0.45,0.55,0.65,  0.30,0.45,0.55,0.65,0.75),nrow = 3, byrow = TRUE)
# p.true.eff <- matrix(c(0.05,0.25,0.50,0.55,0.60,  0.25,0.50,0.55,0.60,0.65,  0.50,0.55,0.60,0.65,0.70),nrow = 3, byrow = TRUE)
# locrm12(p.true.tox,p.true.eff,target.tox=0.35,cutoff.tox=0.85,target.eff=0.2,cutoff.eff=0.9,nmax=51,cohortsize=3,ntrial=100)
