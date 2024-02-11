################################################################################

cpdag_list <- function(list.inc,E){    # E: end of burnIn phase
  L <- list()
  G <- list()
  
 
  nodes <- dim(list.inc[[1]][[1]])[1]
  
  
  mat.sum <- matrix(numeric(nodes*nodes),nrow=nodes)
  
 
  for (i in E:length(list.inc[[1]])){
   
    k <- cpdag(list.inc[[1]][[i]])
    dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
    
    if(length(nrow(k))!=0){
      dummy[k[,1]] <- k[,5]
      L[[i]] <- dummy
    }
    if(length(nrow(k))==0 && length(k)>0){
      dummy[k[1]] <- k[5]
      L[[i]] <- dummy
    }
    mat.com <-matrix(numeric(nodes*nodes),nrow=nodes)
    mat.re <- matrix(numeric(nodes*nodes),nrow=nodes)
    com <- which(L[[i]]>0)
    re  <- which(L[[i]]<0)
    
    mat.com[com] <- 1
    mat.re[re]   <- 1
    
    mat <- mat.com + mat.re + t(mat.re)
    
    G[[i]] <- mat
    
    mat.sum <- mat.sum + mat
  }
  return(list(L,G, (mat.sum/(length(list.inc[[1]])- E+1))))
}


################################################################################

extract_cpdag_from_dag <- function(true_incidence){    
  
  L <- list()
  
  nodes <- dim(true_incidence)[1]
  
  k <- cpdag(true_incidence)
  
  
  dummy <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  if(length(nrow(k))!=0){
    dummy[k[,1]] <- k[,5]
    L <- dummy
  }
  
  if(length(nrow(k))==0 && length(k)>0){
    dummy[k[1]] <- k[5]
    L <- dummy
  }
  
  mat.com <- matrix(numeric(nodes*nodes),nrow=nodes)
  mat.re  <- matrix(numeric(nodes*nodes),nrow=nodes)
  
  com <- which(L>0)
  re  <- which(L<0)
  
  mat.com[com] <- 1
  mat.re[re] <- 1
  mat <- mat.com + mat.re + t(mat.re)
  
  return(mat)
}


##########################################################################################