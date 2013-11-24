adapt_fICA<-function(X, gs=gf, dgs=dgf, name=gnames, kj=2, eps=1e-06, maxiter=100)
{
   n<-dim(X)[1]
   p<-dim(X)[2]
   ng<-length(gs) 
   usedg<-NULL 
   init_est<-"k-JADE"   

   if(!(kj %in% 1:p)){
     W0<-FOBI(X)$W  
     init_est<-"FOBI"  
   }else{
     W0<-k_JADE(X,k=kj,eps=eps,maxiter=maxiter)$W 
     init_est<-paste(kj,"-JADE",sep="")
   }   
   Z <- tcrossprod(X,W0)
   Z <- sweep(Z,2,colMeans(Z))

   alphas <- matrix(0,ng,p) 
   for(i in 1:p){
    alphas[,i]<-compute_alphas(Z[,i],gs,dgs)
   }
   ca <- alphas 

   ord<-NULL
   VN<-diag(p)
   V<-matrix(0,p,p)
 
   for(i in 1:(p-1)){
     mina<-which(ca==min(ca))
     comp<-ceiling(mina/ng)
     gc<-mina-(comp-1)*ng
     vn<-VN[,comp]
     iter<-0
     it0<-1
     if(min(ca)==Inf){ stop("no convergence")
     }else{
     while(TRUE){
       iter<-iter+1
       v<-vn
       vn<-colMeans(sweep(Z,1,gs[[gc]](crossprod(t(Z),v)),"*"))-mean(dgs[[gc]](crossprod(t(Z),v)))*v
       vn<-vn/sqrt(sum(vn^2))
       if(sqrt(sum((vn-v)^2))>1) vn<--vn
    
       vn<-(1-1/(iter+20))*vn+v/(iter+20)
       vn<-vn-crossprod(tcrossprod(V,V),vn) 
       vn<-vn/sqrt(sum(vn^2))
    
       if(sqrt(sum((v-vn)^2))<eps || sqrt(sum((v+vn)^2))<eps){
         usedg[i] <- name[mina-(comp-1)*ng]
         ord[i] <- comp
         ca[,comp] <- Inf 
         break
       }  
       if(iter==maxiter){
         ca[mina]<-Inf
         mina<-which(ca==min(ca))[1]
         comp<-ceiling(mina/ng)
         iter<-0   
         vn<-VN[,comp]
         if(min(ca)==Inf) stop("no convergence")
       }
     }
     V[,i]<-t(vn)
     }
   }
   vn<-VN[,1]   
   vn<-vn-crossprod(tcrossprod(V,V),vn)   
   vn<-vn/sqrt(sum(vn^2))
   V[,p]<-t(vn)

   cnam<-NULL
   for(i in 1:p){
     cnam[i]<-paste("comp",i)
   }

   if(length(ord)==(p-1)){
     ord[p] <- sum(1:p)-sum(ord)
   }else ord <- 1:p

   W<-crossprod(V,W0)
   W<-crossprod(diag(sign(rowMeans(W))),W)
   S<-tcrossprod(sweep(X,2,colMeans(X)),W)

   alph<-matrix(alphas[,ord],ng,p)
   dimnames(alph)<-list(name, cnam)
   res<-list(W=W, gs=name, used_gs=usedg, alphas=alph, init_est=init_est, S=S)
   class(res)<-"bss"
   res
}

compute_alphas<-function(x,gs,dgs)
{
   alphas<-NULL
   for(i in 1:length(gs)){
   Eg<-mean(gs[[i]](x))
   Eg2<-mean(gs[[i]](x)^2)
   Egx<-mean(gs[[i]](x)*x)
   Edg<-mean(dgs[[i]](x))
   if((Egx-Edg)==0){
    alphas[i]<-Inf
   }else alphas[i]<-(Eg2-Eg^2-Egx^2)/(Egx-Edg)^2
    
   if(alphas[i]<0) alphas[i]<-Inf
  }
  alphas
}

