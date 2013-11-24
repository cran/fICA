mat.sqrt<-function(A)
{
 eig<-eigen(A, symmetric=TRUE)
 eig$vectors%*%(diag(eig$values^(1/2)))%*%t(eig$vectors)
}

mat.norm<-function(A)
{
 sqrt(sum(A^2))
}

g_pow3 <- function(x){x^3}
g_tanh <- function(x){tanh(x)}
g_gaus <- function(x){(x)*exp(-(x)^2/2)}

g_rt <- function(x){ifelse(x>0,(x-0)^2,0)}
g_lt <- function(x){ifelse(x<(-0),(x+0)^2,0)}
g_bt <- function(x){g_rt(x)-g_lt(x)}

g_rt0.2 <- function(x){ifelse(x>0.2,(x-0.2)^2,0)}
g_lt0.2 <- function(x){ifelse(x<(-0.2),(x+0.2)^2,0)}
g_bt0.2 <- function(x){g_rt0.2(x)-g_lt0.2(x)}

g_rt0.4 <- function(x){ifelse(x>0.4,(x-0.4)^2,0)}
g_lt0.4 <- function(x){ifelse(x<(-0.4),(x+0.4)^2,0)}
g_bt0.4 <- function(x){g_rt0.4(x)-g_lt0.4(x)}

g_rt0.6 <- function(x){ifelse(x>0.6,(x-0.6)^2,0)}
g_lt0.6 <- function(x){ifelse(x<(-0.6),(x+0.6)^2,0)}
g_bt0.6 <- function(x){g_rt0.6(x)-g_lt0.6(x)}

g_rt0.8 <- function(x){ifelse(x>0.8,(x-0.8)^2,0)}
g_lt0.8 <- function(x){ifelse(x<(-0.8),(x+0.8)^2,0)}
g_bt0.8 <- function(x){g_rt0.8(x)-g_lt0.8(x)}

g_rt1.0 <- function(x){ifelse(x>1.0,(x-1.0)^2,0)}
g_lt1.0 <- function(x){ifelse(x<(-1.0),(x+1.0)^2,0)}
g_bt1.0 <- function(x){g_rt1.0(x)-g_lt1.0(x)}

g_rt1.2 <- function(x){ifelse(x>1.2,(x-1.2)^2,0)}
g_lt1.2 <- function(x){ifelse(x<(-1.2),(x+1.2)^2,0)}
g_bt1.2 <- function(x){g_rt1.2(x)-g_lt1.2(x)}

g_rt1.4 <- function(x){ifelse(x>1.4,(x-1.4)^2,0)}
g_lt1.4 <- function(x){ifelse(x<(-1.4),(x+1.4)^2,0)}
g_bt1.4 <- function(x){g_rt1.4(x)-g_lt1.4(x)}

g_rt1.6 <- function(x){ifelse(x>1.6,(x-1.6)^2,0)}
g_lt1.6 <- function(x){ifelse(x<(-1.6),(x+1.6)^2,0)}
g_bt1.6 <- function(x){g_rt1.6(x)-g_lt1.6(x)}


gf<-c(g_pow3,g_tanh,g_gaus,g_lt0.6,g_rt0.6,g_bt,g_bt0.2,g_bt0.4,g_bt0.6,
g_bt0.8,g_bt1.0,g_bt1.2,g_bt1.4,g_bt1.6)

dg_pow3 <- function(x){3*x^2}
dg_tanh <- function(x){1-tanh(x)^2}
dg_gaus <- function(x){exp(-(x)^2/2)-(x)^2*exp(-(x)^2/2)}

dg_rt <- function(x){ifelse(x>0,2*(x-0),0)}
dg_lt <- function(x){ifelse(x<(-0),2*(x+0),0)}
dg_bt <- function(x){dg_rt(x)-dg_lt(x)}

dg_rt0.2 <- function(x){ifelse(x>0.2,2*(x-0.2),0)}
dg_lt0.2 <- function(x){ifelse(x<(-0.2),2*(x+0.2),0)}
dg_bt0.2 <- function(x){dg_rt0.2(x)-dg_lt0.2(x)}

dg_rt0.4 <- function(x){ifelse(x>0.4,2*(x-0.4),0)}
dg_lt0.4 <- function(x){ifelse(x<(-0.4),2*(x+0.4),0)}
dg_bt0.4 <- function(x){dg_rt0.4(x)-dg_lt0.4(x)}

dg_rt0.6 <- function(x){ifelse(x>0.6,2*(x-0.6),0)}
dg_lt0.6 <- function(x){ifelse(x<(-0.6),2*(x+0.6),0)}
dg_bt0.6 <- function(x){dg_rt0.6(x)-dg_lt0.6(x)}

dg_rt0.8 <- function(x){ifelse(x>0.8,2*(x-0.8),0)}
dg_lt0.8 <- function(x){ifelse(x<(-0.8),2*(x+0.8),0)}
dg_bt0.8 <- function(x){dg_rt0.8(x)-dg_lt0.8(x)}

dg_rt1.0 <- function(x){ifelse(x>1.0,2*(x-1.0),0)}
dg_lt1.0 <- function(x){ifelse(x<(-1.0),2*(x+1.0),0)}
dg_bt1.0 <- function(x){dg_rt1.0(x)-dg_lt1.0(x)}

dg_rt1.2 <- function(x){ifelse(x>1.2,2*(x-1.2),0)}
dg_lt1.2 <- function(x){ifelse(x<(-1.2),2*(x+1.2),0)}
dg_bt1.2 <- function(x){dg_rt1.2(x)-dg_lt1.2(x)}

dg_rt1.4 <- function(x){ifelse(x>1.4,2*(x-1.4),0)}
dg_lt1.4 <- function(x){ifelse(x<(-1.4),2*(x+1.4),0)}
dg_bt1.4 <- function(x){dg_rt1.4(x)-dg_lt1.4(x)}

dg_rt1.6 <- function(x){ifelse(x>1.6,2*(x-1.6),0)}
dg_lt1.6 <- function(x){ifelse(x<(-1.6),2*(x+1.6),0)}
dg_bt1.6 <- function(x){dg_rt1.6(x)-dg_lt1.6(x)}


dgf<-c(dg_pow3,dg_tanh,dg_gaus,dg_lt0.6,dg_rt0.6,dg_bt,dg_bt0.2,dg_bt0.4,
dg_bt0.6,dg_bt0.8,dg_bt1.0,dg_bt1.2,dg_bt1.4,dg_bt1.6)

gnames<-c("pow3","tanh","gaus","lt0.6","rt0.6","bt","bt0.2","bt0.4",
"bt0.6","bt0.8","bt1.0","bt1.2","bt1.4","bt1.6")










