fICA <- function(X, g="tanh", dg=NULL, G=NULL, init=NULL, n.init=1, method="sym2", inR=TRUE, eps=1e-06, maxiter=1000)
{
   n <- nrow(X)
   p <- ncol(X)
   eps <- p*eps 
   name <- c("pow3","tanh","gaus") 

   method <- match.arg(method,c("sym2","sym","def"))
    
   S0 <- cov(X)
   EVD <- eigen(S0, symmetric = TRUE)
   S0.5 <- EVD$vectors %*% tcrossprod(diag(EVD$values^(0.5)), EVD$vectors)
   S0.5inv <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
   Z <- tcrossprod(sweep(X,2,colMeans(X)), S0.5inv)

    
   if(is.null(init)){
     VN <- diag(p)
   }else{ 
     VN <- crossprod(t(init),S0.5)
     VN <- crossprod(solve(mat.sqrt(tcrossprod(VN,VN))),VN)
   }

   V <- switch(method,
        "sym2"={ 
               if(inR){
                if(!(is.function(g)&&is.function(dg)&&is.function(G))){                  
                  gi <- which(name==g[1])
                  g1 <- gf[[gi]]
                  dg1 <- dgf[[gi]]    
                  G1 <- Gf[[gi]] 
                }else{
                  g1 <- g
                  dg1 <- dg
                  G1 <- G
                }
                
                V <- fICA.sym2(Z, VN, g=g1, dg=dg1, G=G1, n.init=n.init, eps=p*eps, maxiter=maxiter)
              }else{
               if(!(is.function(g)&&is.function(dg)&&is.function(G))){                  
                  gi <- which(name==g[1])
                  g1 <- gf[[gi]]
                  dg1 <- dgf[[gi]]    
                  G1 <- Gf[[gi]] 
                }else{
                  g1 <- g
                  dg1 <- dg
                  G1 <- G
                }
  
                warning("c++ code not available for method sym2 yet. Computations are done in R.")
                V <- fICA.sym2(Z, VN, g=g1, dg=dg1, G=G1, n.init=n.init, eps=p*eps, maxiter=maxiter)
               }
              } 
        ,
        "sym"={ 
               if(inR){
                if(!(is.function(g)&&is.function(dg))){                  
                  gi <- which(name==g[1])
                  g1 <- gf[[gi]]
                  dg1 <- dgf[[gi]] 
                  G1 <- Gf[[gi]]     
                }else{
                  g1 <- g
                  dg1 <- dg
                  G1 <- G 
                }
                
                V <- fICA.sym(Z, VN, g=g1, dg=dg1, G=G1, n.init=n.init, eps=p*eps, maxiter=maxiter)
              }else{
                gi <- which(name==g[1])
                V <- .Call("ficasym",Z,gi,VN,p*eps,maxiter,PACKAGE="fICA")
               }
              } 
        ,
        "def"={
                if(inR){
                 if(!(is.function(g)&&is.function(dg))){
                   gi <- which(name==g[1])
                   g1 <- gf[[gi]]
                   dg1 <- dgf[[gi]]    
                 }else{
                   g1 <- g
                   dg1 <- dg
                 }
                     
                 V <- fICA.def(Z, VN, g=g1, dg=dg1, eps=eps, maxiter=maxiter)
                }else{
                 gi <- which(name==g[1])
                 V <- .Call("ficadef",Z,gi,VN,eps,maxiter,PACKAGE="fICA")
                }
              }
         )
     
   if(sum(abs(V$W))>0){    
    W <- crossprod(V$W,S0.5inv)
    W <- crossprod(diag(sign(rowMeans(W))),W)
    S <- tcrossprod(sweep(X,2,colMeans(X)),W)
   }else stop("no convergence")
     
   res <- list(W=W, g=g,method=method, S=S)
   class(res) <- "bss"
   res
}

fICA.sym2<-function(Z, VN, g, dg, G, n.init, eps, maxiter)
{
   n <- nrow(Z)
   p <- ncol(Z)           
   ofm <- 0
   Vm <- VN
   itr <- FALSE
   if(n.init < 1.5){
    V0 <- VN
    V <- VN
    W <- matrix(0,p,p)

    iter <- 0
    while (TRUE){
      iter <- iter+1
  
      ZV <- tcrossprod(Z,V)
      W <- colMeans(G(ZV))*(crossprod(g(ZV),Z)/n-crossprod(diag(colMeans(dg(ZV))),V))
          
      V <- crossprod(solve(mat.sqrt(tcrossprod(W,W))),W)
      if(mat.norm(abs(V)-abs(V0))<eps) break
      if(iter==maxiter) stop("maxiter reached without convergence")
      V0 <- V
    } 
   }else{
    for(i in 1:n.init){
     V0 <- rort(p)
     V <- V0
     W <- matrix(0,p,p)

     iter <- 0
     while (TRUE){
      iter <- iter+1
  
      ZV <- tcrossprod(Z,V)
      ofv <- sum((apply(G(ZV),2,mean))^2)
      if(ofv > ofm){
        ofm <- ofv
        Vm <- V 
      }
      W <- colMeans(G(ZV))*(crossprod(g(ZV),Z)/n-crossprod(diag(colMeans(dg(ZV))),V))
           
      V <- crossprod(solve(mat.sqrt(tcrossprod(W,W))),W)
      if(mat.norm(abs(V)-abs(V0))<eps){
        itr <- TRUE
        ofm <- ofv 
        break
      }
      if(iter==maxiter){
        break
      }
      V0 <- V
     }
     if(itr) break 
    }
    if(!itr){
      warning("no convergence")
      V <- Vm
    } 
   }

   list(W=t(V)) 
} 


fICA.sym <- function(Z, VN, g, dg, G, n.init, eps, maxiter)
{
   n <- nrow(Z)
   p <- ncol(Z)  
   ofm <- 0
   Vm <- diag(p)
   itr <- FALSE
   if(n.init < 1.5){
    V0 <- VN
    V <- VN
    W <- matrix(0,p,p)

    iter <- 0
    while (TRUE){
      iter <- iter+1
  
      ZV <- tcrossprod(Z,V)
      W <- crossprod(g(ZV),Z)/n-crossprod(diag(colMeans(dg(ZV))),V)
           
      V <- crossprod(solve(mat.sqrt(tcrossprod(W,W))),W)
      if(mat.norm(abs(V)-abs(V0))<eps) break
      if(iter==maxiter) stop("maxiter reached without convergence")
      V0 <- V
    } 
   }else{
    for(i in 1:n.init){
     V0 <- rort(p)
     V <- V0
     W <- matrix(0,p,p)

     iter <- 0
     while (TRUE){
      iter <- iter+1
  
      ZV <- tcrossprod(Z,V)
      ofv <- sum(abs(apply(G(ZV),2,mean)))
      if(ofv>ofm){
        ofm <- ofv
        Vm <- V 
      }
      W <- crossprod(g(ZV),Z)/n-crossprod(diag(colMeans(dg(ZV))),V)
           
      V <- crossprod(solve(mat.sqrt(tcrossprod(W,W))),W)
      if(mat.norm(abs(V)-abs(V0))<eps){
        itr <- TRUE
        ofm <- ofv 
        break
      }
      if(iter==maxiter){
        break
      }
      V0 <- V
     }
     if(itr) break 
    }
    if(!itr){
      warning("no convergence")
      V <- Vm
    } 
   }

   list(W=t(V)) 
} 

fICA.def <- function(Z, VN, g, dg, eps, maxiter)
{
   p <- ncol(Z)
   V <- matrix(0,p,p)
  
   for(i in 1:p){
     vn <- VN[,i]
     a <- 0
     iter <- 0
     while(TRUE){
       iter <- iter+1 
       v <- vn
       Zv <- crossprod(t(Z),v)
       vn <- colMeans(sweep(Z,1,g(Zv),"*"))-mean(dg(Zv))*v
       vn <- vn/sqrt(sum(vn^2))
       if(sqrt(sum((vn-v)^2))>1) vn <- -vn
       if((a>0)&&(floor(iter/a)==iter/a)){
        vn <- (1-1/5)*vn+v/5
       }else vn <- (1-1/(iter+20))*vn+v/(iter+20) 
       vn <- vn-crossprod(tcrossprod(V,V),vn)  
       vn <- vn/sqrt(sum(vn^2))
       if(sqrt(sum((v-vn)^2))<eps || sqrt(sum((v+vn)^2))<eps) break
       if(iter==maxiter){
        a <- a+1
        iter <- 0
        if(a>10) stop("maxiter reached without convergence")  
       } 
      }
     V[,i] <- t(vn)
   }
   list(W=V) 
}


