#include "fica.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;
  
   vec g1(vec x){
    return x%x%x;    
   }

   vec dg1(vec x){
    return 3*x%x;    
   }

   vec g2(vec x){
    return tanh(x);    
   }

   vec dg2(vec x){
    return 1-tanh(x)%tanh(x);    
   }

   vec g3(vec x){
    return x%exp(-x%x/2);    
   }

   vec dg3(vec x){
    return exp(-x%x/2)-x%x%exp(-x%x/2);    
   }

   mat g1m(mat x){
    return x%x%x;    
   }

   mat dg1m(mat x){
    return 3*x%x;    
   }

   mat g2m(mat x){
    return tanh(x);    
   }

   mat dg2m(mat x){
    return 1-tanh(x)%tanh(x);    
   }

   mat g3m(mat x){
    return x%exp(-x%x/2);    
   }

   mat dg3m(mat x){
    return exp(-x%x/2)-x%x%exp(-x%x/2);    
   }

   mat grm(mat X, double a){
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)>a)
       gx(i,j) = (X(i,j)-a)*(X(i,j)-a); 
      }
     }    
    return gx;
   }
  
   mat dgrm(mat X, double a){
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)>a)
       gx(i,j) = 2*(X(i,j)-a); 
      }
     }    
     return gx;
   }
   
   mat glm(mat X, double a){
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)<(-a))
       gx(i,j) = (X(i,j)+a)*(X(i,j)+a); 
      }
     }    
     return gx;
   }
  
   mat dglm(mat X, double a){
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)<(-a))
       gx(i,j) = 2*(X(i,j)+a); 
      }
     }    
     return gx;
   }

    mat gbm(mat X, double a){
     return grm(X,a)-glm(X,a);
   }
  
   mat dgbm(mat X, double a){
     return dgrm(X,a)-dglm(X,a);
   }

   vec gl(vec x, double a){
    unsigned int i;
    vec res(x.size());
    res.fill(0);
    for(i=0; i<x.size(); i++){
     if(x(i)<(-a))
      res(i) = (x(i)+a)*(x(i)+a); 
    }    
    return res;
   }
  
   vec dgl(vec x, double a){
    unsigned int i;
    vec res(x.size());
    res.fill(0);
    for(i=0; i<x.size(); i++){
     if(x(i)<(-a))
      res(i) = 2*(x(i)+a); 
    }    
    return res;
   }

   vec gr(vec x, double a){
    unsigned int i;
    vec res(x.size());
    res.fill(0);
    for(i=0; i<x.size(); i++){
     if(x(i)>a)
      res(i) = (x(i)-a)*(x(i)-a); 
    }    
    return res;
   }
  
   vec dgr(vec x, double a){
    unsigned int i;
    vec res(x.size());
    res.fill(0);
    for(i=0; i<x.size(); i++){
     if(x(i)>a)
      res(i) = 2*(x(i)-a); 
    }    
    return res;
   }

   vec gb(vec x, double a){
    return gr(x, a)-gl(x, a);
   }  

   vec dgb(vec x, double a){
    return dgr(x, a)-dgl(x, a);
   }
  
   mat comp_alph(mat X){
    int p=X.n_cols;
    int j;
    int k;
 
    vec Egs = zeros(14*4*p);
      
    mat alphas = zeros(14,p);    

    for(j=0; j<p; j++){ 
     Egs(j*56+0) = mean(g1(X.col(j)));
     Egs(j*56+1) = mean(g1(X.col(j))%g1(X.col(j)));
     Egs(j*56+2) = mean(g1(X.col(j))%X.col(j));
     Egs(j*56+3) = mean(dg1(X.col(j))); 

     Egs(j*56+4) = mean(g2(X.col(j)));
     Egs(j*56+5) = mean(g2(X.col(j))%g2(X.col(j)));
     Egs(j*56+6) = mean(g2(X.col(j))%X.col(j));
     Egs(j*56+7) = mean(dg2(X.col(j))); 
 
     Egs(j*56+8) = mean(g3(X.col(j)));
     Egs(j*56+9) = mean(g3(X.col(j))%g3(X.col(j)));
     Egs(j*56+10) = mean(g3(X.col(j))%X.col(j));
     Egs(j*56+11) = mean(dg3(X.col(j))); 

     Egs(j*56+12) = mean(gl(X.col(j),0.6));
     Egs(j*56+13) = mean(gl(X.col(j),0.6)%gl(X.col(j),0.6));
     Egs(j*56+14) = mean(gl(X.col(j),0.6)%X.col(j));
     Egs(j*56+15) = mean(dgl(X.col(j),0.6)); 
    
     Egs(j*56+16) = mean(gr(X.col(j),0.6));
     Egs(j*56+17) = mean(gr(X.col(j),0.6)%gr(X.col(j),0.6));
     Egs(j*56+18) = mean(gr(X.col(j),0.6)%X.col(j));
     Egs(j*56+19) = mean(dgr(X.col(j),0.6)); 

     Egs(j*56+20) = mean(gb(X.col(j),0));
     Egs(j*56+21) = mean(gb(X.col(j),0)%gb(X.col(j),0));
     Egs(j*56+22) = mean(gb(X.col(j),0)%X.col(j));
     Egs(j*56+23) = mean(dgb(X.col(j),0)); 
    
     Egs(j*56+24) = mean(gb(X.col(j),0.2));
     Egs(j*56+25) = mean(gb(X.col(j),0.2)%gb(X.col(j),0.2));
     Egs(j*56+26) = mean(gb(X.col(j),0.2)%X.col(j));
     Egs(j*56+27) = mean(dgb(X.col(j),0.2)); 
   
     Egs(j*56+28) = mean(gb(X.col(j),0.4));
     Egs(j*56+29) = mean(gb(X.col(j),0.4)%gb(X.col(j),0.4));
     Egs(j*56+30) = mean(gb(X.col(j),0.4)%X.col(j));
     Egs(j*56+31) = mean(dgb(X.col(j),0.4)); 

     Egs(j*56+32) = mean(gb(X.col(j),0.6));
     Egs(j*56+33) = mean(gb(X.col(j),0.6)%gb(X.col(j),0.6));
     Egs(j*56+34) = mean(gb(X.col(j),0.6)%X.col(j));
     Egs(j*56+35) = mean(dgb(X.col(j),0.6)); 

     Egs(j*56+36) = mean(gb(X.col(j),0.8));
     Egs(j*56+37) = mean(gb(X.col(j),0.8)%gb(X.col(j),0.8));
     Egs(j*56+38) = mean(gb(X.col(j),0.8)%X.col(j));
     Egs(j*56+39) = mean(dgb(X.col(j),0.8)); 

     Egs(j*56+40) = mean(gb(X.col(j),1));
     Egs(j*56+41) = mean(gb(X.col(j),1)%gb(X.col(j),1));
     Egs(j*56+42) = mean(gb(X.col(j),1)%X.col(j));
     Egs(j*56+43) = mean(dgb(X.col(j),1)); 

     Egs(j*56+44) = mean(gb(X.col(j),1.2));
     Egs(j*56+45) = mean(gb(X.col(j),1.2)%gb(X.col(j),1.2));
     Egs(j*56+46) = mean(gb(X.col(j),1.2)%X.col(j));
     Egs(j*56+47) = mean(dgb(X.col(j),1.2)); 
   
     Egs(j*56+48) = mean(gb(X.col(j),1.4));
     Egs(j*56+49) = mean(gb(X.col(j),1.4)%gb(X.col(j),1.4));
     Egs(j*56+50) = mean(gb(X.col(j),1.4)%X.col(j));
     Egs(j*56+51) = mean(dgb(X.col(j),1.4)); 

     Egs(j*56+52) = mean(gb(X.col(j),1.6));
     Egs(j*56+53) = mean(gb(X.col(j),1.6)%gb(X.col(j),1.6));
     Egs(j*56+54) = mean(gb(X.col(j),1.6)%X.col(j));
     Egs(j*56+55) = mean(dgb(X.col(j),1.6)); 

     for(k=0; k<14; k++){
      if(Egs(j*56+k*4+2)!=Egs(j*56+k*4+3))
       alphas(k,j) = (Egs(j*56+k*4+1)-Egs(j*56+k*4)*Egs(j*56+k*4)-Egs(j*56+k*4+2)*
Egs(j*56+k*4+2))/((Egs(j*56+k*4+2)-Egs(j*56+k*4+3))*(Egs(j*56+k*4+2)-Egs(j*56+k*4+3)));
     }
    }
    return alphas;
}


   mat sweep(mat X, vec y)  
  {
    int i;
    int n = X.n_rows; 
    int p = X.n_cols;
    mat Xy(n,p);
    for(i=0; i<p; i++) 
     Xy.col(i) = X.col(i)%y;
    
    return Xy;
}

   vec get_comp(mat X, int g, double eps, int maxiter, vec wn0, mat W)
  {
    int n = X.n_rows;
    int p = X.n_cols;     
    int iter = 0;
    int a;
    vec wn = wn0;
    vec w(p);   
    vec xw(n);  
   
    if(g==0){
     a=0;
     while(1){ 
      iter++;
      w=wn;
      
      wn=mean(trans(sweep(X,g1(X*w))),1)-3*w;
       
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }     
     
    }else if(g==1){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,g2(xw))),1)-mean(dg2(xw))*w;
     
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }    
     
    }else if(g==2){
     a = 0;
     while(1){
      iter++;
      w = wn;
    
      xw = X*w;
      wn = mean(trans(sweep(X,g3(xw))),1)-mean(dg3(xw))*w;

      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }     
     
    }else if(g==3){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gl(xw,0.6))),1)-mean(dgl(xw,0.6))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }     
     
    }else if(g==4){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gr(xw,0.6))),1)-mean(dgr(xw,0.6))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }     
     
    }else if(g==5){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,0))),1)-mean(dgb(xw,0))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }
     
    }else if(g==6){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,0.2))),1)-mean(dgb(xw,0.2))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     } 
    
    }else if(g==7){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,0.4))),1)-mean(dgb(xw,0.4))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }
    
    }else if(g==8){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,0.6))),1)-mean(dgb(xw,0.6))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }
    
    }else if(g==9){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,0.8))),1)-mean(dgb(xw,0.8))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }
    
    }else if(g==10){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,1))),1)-mean(dgb(xw,1))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }
    
    }else if(g==11){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,1.2))),1)-mean(dgb(xw,1.2))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }
    
    }else if(g==12){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,1.4))),1)-mean(dgb(xw,1.4))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
      }
     }
    
    }else if(g==13){
     a = 0;
     while(1){
      iter++;
      w = wn;
      xw = X*w;
      wn = mean(trans(sweep(X,gb(xw,1.6))),1)-mean(dgb(xw,1.6))*w;
      
      wn = wn/sqrt(sum(square(wn)));   
      if(sqrt(sum(square(wn-w)))>1) wn = -wn;
      if((a>0)&&(floor(iter/(a+1))==iter/(a+1))){
        wn = (1-1/5)*wn+w/5;
      }else wn = (1-1/(iter+20))*wn+w/(iter+20); 

      wn = wn-W*W.t()*wn; 
      wn = wn/sqrt(sum(square(wn)));
    
      if(sqrt(sum(square(w-wn)))<eps || sqrt(sum(square(w+wn)))<eps) break;
       
      if(iter==maxiter){
        a += 1;
        iter = 0;
        if(a>10){ 
         wn.fill(0);
         break;
        }
       }
      }
    }  
    return wn;   
}


   vec argmin_mat(mat X)
  {
    int i;
    int j;
    int n = X.n_rows;
    int p = X.n_cols;
    double minX = X(0,0);     
    vec am(2);
    am(0)=0;  
    am(1)=0;

    for(i=0; i<n; i++){
     for(j=0; j<p; j++){
      if(X(i,j)<minX){
       am(0) = i;
       am(1) = j;
       minX = X(i,j); 
      } 
     }
    }
    return am;
  }

   SEXP adfica(SEXP x, SEXP epsilon, SEXP maxiter)
  {
    mat X = as<arma::mat>(x);
    int p = X.n_cols;
    int i;
    int j; 
    double eps = as<double>(epsilon);
    int maxit = as<int>(maxiter); 
    vec wj(p); 
    mat W(p,p);  
    W.fill(0);
    mat Id(p,p);  
    Id.eye();
    vec wn0(p);
    vec am(2);
    int g;
    int comp;
    vec usedg(p-1);
    vec ord(p-1);
 
    mat alphas = comp_alph(X);
    mat alph = alphas; 
    for(i=0; i<14; i++){
     for(j=0; j<p; j++){
      if(alph(i,j)<=0) 
       alph(i,j) = 1000000;
     }
    }
    
    for(j=0; j<(p-1); j++){
      if(min(min(alph,1))>999999){
        W.fill(0);
        break;
      } 
      am = argmin_mat(alph);
      g = am(0);
      comp = am(1);
      wn0 = Id.col(comp);   
      while(1){
       wj = get_comp(X, g, eps, maxit, wn0, W);
       if(sum(abs(wj))>0){
        W.col(j) = wj;
        usedg(j) = g;
        ord(j) = comp;
        alph.col(comp).fill(1000000); 
        break;
       }else{
        alph(g,comp) = 1000000; 
        am = argmin_mat(alph);
        g = am(0);
        comp = am(1);
        wn0 = Id.col(comp);  
        if(min(min(alph,1))>999999){
         W.fill(0);
         break;
        } 
       } 
      }      
    }  
     
    wj = Id.col(0);
    wj = wj-W*W.t()*wj; 
    W.col(p-1) = wj/sqrt(sum(square(wj)));
  
    return Rcpp::List::create(Rcpp::Named("W") = W,
                              Rcpp::Named("alphas") = alphas, 
                              Rcpp::Named("ord") = ord,
                              Rcpp::Named("usedg") = usedg);
  }

  vec comp_alph_rel(mat X, int g){
    int p=X.n_cols;
    int j;
    
    vec Egs = zeros(4*p);
    vec alphas = zeros(p);    

    if(g==0){
     for(j=0; j<p; j++){ 
      Egs(j*4+0) = mean(g1(X.col(j)));
      Egs(j*4+1) = mean(g1(X.col(j))%g1(X.col(j)));
      Egs(j*4+2) = mean(g1(X.col(j))%X.col(j));
      Egs(j*4+3) = mean(dg1(X.col(j))); 
     }  
    }else if(g==1){
     for(j=0; j<p; j++){ 
      Egs(j*4+0) = mean(g2(X.col(j)));
      Egs(j*4+1) = mean(g2(X.col(j))%g2(X.col(j)));
      Egs(j*4+2) = mean(g2(X.col(j))%X.col(j));
      Egs(j*4+3) = mean(dg2(X.col(j))); 
     } 
    }else{
     for(j=0; j<p; j++){ 
      Egs(j*4+0) = mean(g3(X.col(j)));
      Egs(j*4+1) = mean(g3(X.col(j))%g3(X.col(j)));
      Egs(j*4+2) = mean(g3(X.col(j))%X.col(j));
      Egs(j*4+3) = mean(dg3(X.col(j))); 
     }
    }
    for(j=0; j<p; j++){ 
     if(Egs(j*4+2)!=Egs(j*4+3))
      alphas(j) = (Egs(j*4+1)-Egs(j*4)*Egs(j*4)-Egs(j*4+2)*Egs(j*4+2))/
                  ((Egs(j*4+2)-Egs(j*4+3))*(Egs(j*4+2)-Egs(j*4+3)));
    }
 
    return alphas;
  }


  int argmin_vec(vec x)
  {
    int j;
    int p = x.size();
    double minx = x(0);     
    int am = 0;

    for(j=0; j<p; j++){
     if(x(j)<minx){
       am = j;
       minx = x(j); 
     } 
    }
    return am;
  }


  SEXP relfica(SEXP x, SEXP g, SEXP epsilon, SEXP maxiter)
  {
    mat X = as<arma::mat>(x);
    int p = X.n_cols;
    int j; 
    double eps = as<double>(epsilon);
    int maxit = as<int>(maxiter); 
    int gi = as<int>(g)-1;
    vec wj(p); 
    mat W(p,p);  
    W.fill(0);
    mat Id(p,p);  
    Id.eye();
    vec wn0(p);
    int comp;
    vec ord(p-1);
 
    vec alphas = comp_alph_rel(X, gi);
    vec alph = alphas; 
    for(j=0; j<p; j++){
     if(alph(j)<=0) 
      alph(j) = 1000000;
    }
    
    for(j=0; j<(p-1); j++){
      if(min(alph)>999999){
        W.fill(0);
        break;
      } 
      comp = argmin_vec(alph);
      wn0 = Id.col(comp);   
      while(1){
       wj = get_comp(X, gi, eps, maxit, wn0, W);
       if(sum(abs(wj))>0){
        W.col(j) = wj;
        ord(j) = comp;
        alph(comp) = 1000000; 
        break;
       }else{
        alph(comp) = 1000000; 
        comp = argmin_vec(alph);       
        wn0 = Id.col(comp);  
        if(min(alph)>999999){
         W.fill(0);
         break;
        } 
       } 
      }      
    }  
     
    wj = Id.col(0);
    wj = wj-W*W.t()*wj; 
    W.col(p-1) = wj/sqrt(sum(square(wj)));
  
    return Rcpp::List::create(Rcpp::Named("W") = W,
                              Rcpp::Named("alphas") = alphas, 
                              Rcpp::Named("ord") = ord);
  }

  SEXP ficadef(SEXP x, SEXP g, SEXP w0, SEXP epsilon, SEXP maxiter)
  {
    mat X = as<arma::mat>(x);
    mat W0 = as<arma::mat>(w0);
    int p = X.n_cols;
    int j; 
    double eps = as<double>(epsilon);
    int maxit = as<int>(maxiter); 
    int gi = as<int>(g)-1; 
    vec wj(p); 
    mat W(p,p);  
    W.fill(0);
    vec wn0(p);
          
    for(j=0; j<(p-1); j++){
      wn0 = W0.col(j);   
      while(1){
       wj = get_comp(X, gi, eps, maxit, wn0, W);
       if(sum(abs(wj))>0){
        W.col(j) = wj;
        break;
       }else{
         W.fill(0);
         break;
        } 
      }        
    }  
     
    wj = W0.col(p-1);
    wj = wj-W*W.t()*wj; 
    W.col(p-1) = wj/sqrt(sum(square(wj)));
  
    return Rcpp::List::create(Rcpp::Named("W") = W);
  }

  mat msqrt(mat X){
   int p = X.n_cols;
   vec eval(p);
   mat evec(p,p);
   eig_sym(eval,evec,X);
   mat Sqrt = evec*diagmat(sqrt(eval))*evec.t(); 
   return(Sqrt);
  }
  

  SEXP ficasym(SEXP x, SEXP g, SEXP w0, SEXP epsilon, SEXP maxiter)
  {
    mat X = as<arma::mat>(x);
    mat W0 = as<arma::mat>(w0);
    int n = X.n_rows;
    int p = X.n_cols;
    double eps = as<double>(epsilon);
    int maxit = as<int>(maxiter); 
    int gi = as<int>(g)-1; 
    int iter = 0;
    mat W = W0;
    mat Id(p,p);  
    Id.eye();
    vec wn(p);
    vec w(p);            
    mat xw(n,p); 
 
    if(gi==0){
    // a = 0;
     while(1){ 
      iter++;
      xw = X*W.t();
      W = g1m(xw).t()*X/n-diagmat(mean(dg1m(xw),0))*W; 

      W = solve(msqrt(W*W.t()),Id)*W;
      if(sqrt(sum(sum(square(abs(W)-abs(W0)))))<eps) break;
       
      if(iter==maxit){
         W.fill(0);
         break;
      }
      W0 = W;
     }     

    }else if(gi==1){     
     while(1){ 
      iter++;
      xw = X*W.t();
      W = g2m(xw).t()*X/n-diagmat(mean(dg2m(xw),0))*W; 

      W = solve(msqrt(W*W.t()),Id)*W;
      if(sqrt(sum(sum(square(abs(W)-abs(W0)))))<eps) break;
       
      if(iter==maxit){
         W.fill(0);
         break;
      }
      W0 = W;
     }     
    }else{
     while(1){ 
      iter++;
      xw = X*W.t();
      W = g3m(xw).t()*X/n-diagmat(mean(dg3m(xw),0))*W; 

      W = solve(msqrt(W*W.t()),Id)*W;
      if(sqrt(sum(sum(square(abs(W)-abs(W0)))))<eps) break;
       
      if(iter==maxit){
         W.fill(0);
         break;
      }
      W0 = W;
     }     
    }
    return Rcpp::List::create(Rcpp::Named("W") = W.t());
  }

   SEXP gpow3(SEXP x){
    mat X = as<arma::mat>(x);
    mat gx = X%X%X;
    return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }
  
   SEXP dgpow3(SEXP x){
    mat X = as<arma::mat>(x);
    mat gx = 3*X%X;
    return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }
   
   SEXP Gpow3(SEXP x){
    mat X = as<arma::mat>(x);
    mat gx = X%X%X%X/4;
    return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }

   SEXP grn(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)>a)
       gx(i,j) = (X(i,j)-a)*(X(i,j)-a); 
      }
     }    
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }

   SEXP dgrn(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)>a)
       gx(i,j) = 2*(X(i,j)-a); 
      }
     }    
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }

   SEXP gln(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)<(-a))
       gx(i,j) = (X(i,j)+a)*(X(i,j)+a); 
      }
     }    
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }

   SEXP dgln(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)<(-a))
       gx(i,j) = 2*(X(i,j)+a); 
      }
     }    
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }
 
   SEXP gbn(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     mat gx = grm(X, a)-glm(X, a); 
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }

   SEXP dgbn(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     mat gx = dgrm(X, a)-dglm(X, a); 
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   } 

   
   SEXP Grn(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)>a)
       gx(i,j) = (X(i,j)-a)*(X(i,j)-a)*(X(i,j)-a)/3; 
      }
     }    
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }

   SEXP Gln(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)<(-a))
       gx(i,j) = (X(i,j)+a)*(X(i,j)+a)*(X(i,j)+a)/3; 
      }
     }    
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }

   mat Grm(mat X, double a){
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)>a)
       gx(i,j) = (X(i,j)-a)*(X(i,j)-a)*(X(i,j)-a)/3; 
      }
     }    
    return gx;
   }
     
   mat Glm(mat X, double a){
     int n = X.n_rows;
     int p = X.n_cols;
     mat gx(n,p); 
     gx.zeros();
     int i;
     int j; 
     for(i=0; i<n; i++){
      for(j=0; j<p; j++){
       if(X(i,j)<(-a))
       gx(i,j) = (X(i,j)+a)*(X(i,j)+a)*(X(i,j)+a)/3; 
      }
     }    
     return gx;
   }
  
 
   SEXP Gbn(SEXP x, SEXP a1){
     mat X = as<arma::mat>(x);
     double a = as<double>(a1);
     mat gx = Grm(X, a)-Glm(X, a); 
     return Rcpp::List::create(Rcpp::Named("gx") = gx);
   }



