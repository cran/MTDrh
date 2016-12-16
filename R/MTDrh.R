
MTDrh <- function(scenarios, observation, prob=NULL,debias=FALSE,transformation=FALSE)
{
  
  n.instance <- dim(observation)[2]
  n.scen <- dim(scenarios)[2]
  n.length <- dim(observation)[1]
  if (is.null(prob)){
    prob <- array(rep(1/n.scen,times=n.scen*n.instance),dim=c(n.scen,n.instance))}
  if(all(apply(prob,2,sum)<1.0001&apply(prob,2,sum)>0.9999)){

  if (dim(observation)[1]==dim(scenarios)[1]&&dim(observation)[2]==dim(scenarios)[3]&&dim(prob)[1]==dim(scenarios)[2]&&dim(prob)[2]==dim(scenarios)[3]) {
  if (debias){
    z<-array(,dim=c(n.instance,n.length))
    for (i in 1:n.instance)
    {
      a<-rbind(observation[,i],t(scenarios[,,i]))
      for (h in 1:n.length){
        z[i,h]<-mean(as.matrix(a)[-1,h])-as.matrix(a)[1,h]
      }
    }
  }

#calculate MTD ranks and build histogram

  combine <- as.matrix(matrix(nrow=(n.scen+1),ncol=n.length))
  distances<-array(,dim=c(n.scen+1,n.scen+1,n.instance))
  moving.cost<-matrix(,n.scen+1,n.instance)
  rank.actual<-vector(,n.instance)

  for (i in 1:n.instance){
    
    if(n.length==1){
    a<- t(cbind(observation[,i],t(scenarios[,,i])))
    combine<-a
    }else{
    a<-rbind(observation[,i],t(scenarios[,,i]))
    combine<-a
    }

    if(debias){
      combine <- a
      for (s in 2:(n.scen+1)){
        for (h in 1:n.length){
          combine[s,h]<-as.matrix(a)[s,h]-mean(z[,h])
        }
      }

    }

    if(transformation){
      ctr   <- colMeans(combine)
      S     <- cov(combine)
      Seig  <- eigen(S)

      if ( sum(Seig$values<0) == 0 ) {

        sqrtD <- sqrt(Seig$values)
        
        if (n.length==1){
          SsqrtInv <- Seig$vectors %*% (1/sqrtD) %*% t(Seig$vectors)
        }else{
          SsqrtInv <- Seig$vectors %*% diag(1/sqrtD) %*% t(Seig$vectors)
        }
        ctr_matrix<-matrix(nrow=n.scen+1,ncol=n.length)
        for (s in 1:(n.scen+1)){
          ctr_matrix[s,]<-ctr}
        Xdot<-as.matrix(combine)-as.matrix(ctr_matrix)
        combine   <- t(SsqrtInv %*% t(Xdot)) }

    }


    distances[,,i]<-as.matrix(dist(combine,method="euclidean",diag=FALSE,upper=FALSE,p=2))
    for (s in 1:(n.scen+1)) {
      if (s==1) {
        moving.cost[1,i]<-sum(distances[,s,i]*c(0,prob[,i]))
        next}
      moving.cost[s,i]<-sum(distances[,s,i]*c(prob[s-1,i],prob[,i]))
    }
    rank.actual[i]<-(n.scen+1)-rank(moving.cost[,i])[1]+1
  }

  hist(rank.actual, breaks=c(0:(n.scen+1)),xlab="bin",ylab="frequency",col="gray",main="MTD rh",ylim=c(0,n.instance/2))
  return(rank.actual)
  }else{
    print("The number of scenario and observation instances in and the dimension of scenarios and observation should be equal.")
  }
  }else{
    print("The sum of the probabilities of scenarios should equal 1 for each instance.")
  }
}

