fICA<-function(X,g=NULL,dg=NULL,init=NULL,method="sym",eps=1e-06,maxiter=100)
{
   n<-dim(X)[1]
   p<-dim(X)[2]
 
   method<-match.arg(method,c("sym","def"))
    
   invS0<-solve(mat.sqrt(cov(X)))
   Z<-tcrossprod(sweep(X,2,colMeans(X)),invS0) 
  
   if(is.null(init)){
    VN<-diag(p)
   }else VN<-init
 
   if(is.null(g)){
    g<-g_tanh
    dg<-dg_tanh
   }
   
   V<-switch(method,
        "sym"={
               fICA.sym(Z, VN, g=g, dg=dg, eps=eps, maxiter=maxiter)
               }
        ,
        "def"={
               fICA.def(Z, VN, g=g, dg=dg, eps=eps, maxiter=maxiter)
               }
        )
         
    W<-crossprod(V,invS0)
    W<-crossprod(diag(sign(rowMeans(W))),W)
    S<-tcrossprod(sweep(X,2,colMeans(X)),W)
    
    res <- list(W=W, g=g,method=method, S=S)
    class(res) <- "bss"
    res
}

fICA.sym<-function(Z,VN,g,dg,eps,maxiter)
{
   p<-dim(Z)[2]
  
   V0<-VN
   V<-VN
   W<-matrix(0,p,p)

   iter<-0
   while (TRUE){
     iter <- iter+1
     for (i in 1:p){
       w<-V[i,]
       wi<-colMeans(sweep(Z,1,g(crossprod(t(Z),w)),"*"))-mean(dg(crossprod(t(Z),w)))*w 
       W[i,]<-wi
     }  
         
     V<-crossprod(solve(mat.sqrt(tcrossprod(W,W))),W)
     if(mat.norm(abs(V)-abs(V0))<eps) break
     if(iter==maxiter) stop("maxiter reached without convergence")
     V0<-V
   }
   t(V) 
} 

fICA.def<-function(Z,VN,g,dg,eps,maxiter)
{
   p<-dim(Z)[2]
   V<-matrix(0,p,p)
  
   for(i in 1:p){
     vn<-VN[,i]
     for(it in 1:maxiter){
       v<-vn
       vn<-colMeans(sweep(Z,1,g(crossprod(t(Z),v)),"*"))-mean(dg(crossprod(t(Z),v)))*v
       vn<-vn/sqrt(sum(vn^2))
       if(sqrt(sum((vn-v)^2))>1) vn<--vn
       vn<-(1-1/(it+20))*vn+v/(it+20)
       vn<-vn-crossprod(tcrossprod(V,V),vn)  
       vn<-vn/sqrt(sum(vn^2))
       
       if(sqrt(sum((v-vn)^2))<eps || sqrt(sum((v+vn)^2))<eps) break
       if(it==maxiter) stop("maxiter reached without convergence") 
     }
    
     V[,i]<-t(vn)
   }
   V
}




