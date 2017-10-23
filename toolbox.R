library(nnet)

#function that loads the plant-pollinator network
#input: number = label of a network
#output: web = mutualsitic network
load_data <- function(){
    name <- paste('network.csv',sep='')
    d <- read.csv(file=name,header=FALSE)
    web <- as.matrix(d)
    web[web > 0] = 1
    return(web)
}

#computes the raw NODF
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% as.matrix(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)

  N <- as.matrix(web) %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}

#finds the maximizum raw value of NODF with given connectance and community size
#input: web = mutualistic network
#output: the maximun row value of NODF
max_nest <- function(web){
  #binarize the interaction matrix
  web_binary <- web
  web_binary[web_binary > 0] = 1
  #compute the number of pollinators, plants and interactions
  SA <- nrow(web_binary)
  SP <- ncol(web_binary)
  SI <- floor(sum(web_binary))
  #initialize the interaction matrix with minimum requirements
  web_opt <- matrix(0, nrow=SA, ncol=SP)  
  web_opt[1,] <- 1
  web_opt[,1] <- 1
  web_opt[2,2] <- 1
  #counting the number of
  SI_left <- SI-SP-SA
  if(SI_left>0){
    #search the best possible location
    for(j in 1:SI_left){
      #compare all possible locations and the maximum one
      position_potential <- websearch_NODF(web_opt)
      nest_poten <- c()
      for(i in 1:nrow(position_potential)) {
        web_poten <- web_opt
        web_poten[position_potential[i,1],position_potential[i,2]] <- 1
        nest_poten[i] <- nestedness_NODF(web_poten)
      } 
      position_the <- which.is.max(nest_poten)
      web_opt[position_potential[position_the,1],position_potential[position_the,2]] <- 1
    }
    return(nestedness_NODF(web_opt))
  }
  #this is to prevent the trivial case
  else{
	  return(-1)
  }
}

#finds all possible positions to add an interaction
#input: web = mutualistic network
#output: all achievable positions of adding an interaction to the network
websearch_NODF <- function(web){
  SA <- nrow(web)
  SP <- ncol(web)
  domain <- web
  position <- which(domain == 1, arr.ind=T)
  position <- subset(position, position[,2] != 1)
  position <- subset(position, position[,1] != 1)
  boundary <- matrix(0, nrow=2*nrow(position),ncol=2)
  j=1
  #choose boundary points
  for(i in 1:nrow(position)){
    if(position[i,1]<nrow(domain)&&position[i,2]<ncol(domain))
    if(domain[position[i,1]+1,position[i,2]]+
    domain[position[i,1]-1,position[i,2]]+
    domain[position[i,1],position[i,2]+1]+
    domain[position[i,1],position[i,2]-1]<=3){
      boundary[j,1] <- position[i,1]+1
      boundary[j,2] <- position[i,2]
      boundary[j+1,1] <- position[i,1]
      boundary[j+1,2] <- position[i,2]+1
      j <- j+2
      } 
  }
  #delete those with zero entries which entered as auxiliary in the first place
  keep <- c()
  for(i in 1:nrow(boundary)){
    if(boundary[i,1]+boundary[i,2]>0) keep <- append(keep,i)
  }
  boundary <- boundary[keep,]
  #choose true boundary points
  stay <- c()
  for(i in 1:nrow(boundary)){
    if(boundary[i,1]<SA&&boundary[i,2]<SP){
      if(domain[boundary[i,1]+1,boundary[i,2]]+
    domain[boundary[i,1]-1,boundary[i,2]]+
    domain[boundary[i,1],boundary[i,2]+1]+
    domain[boundary[i,1],boundary[i,2]-1]==2
    && domain[boundary[i,1],boundary[i,2]]==0){
      stay <- append(stay,i)
          }
    }
  }
  boundary <- boundary[stay,]
  return(boundary)
}

#calculates the combined NODF statistic 
#inputs: web = mutualistic network, raw NODF, maximum raw NODF
#output: the combined NODF statistic
comb_nest <- function(web,NODF,max_NODF){
  C <- sum(web)/(ncol(web)*nrow(web))
  S <- sqrt(ncol(web) * nrow(web) )
  out <- NODF / (max_NODF * C * log10(S))
  return(out)
}


