implicit_FD<- function(T,dt,dx,gam,lB,uB,type,K,other_asset,x0){
  r = 0.05
  s = 0.4
  time_steps <- as.integer(T/dt)
  space_steps <- as.integer((uB-lB)/dx)+1
  vetS <- lB + dx*(0:(space_steps-1))
  
  V = matrix(0,space_steps,time_steps) # intitialize matrix
  
  vetI <- lB/dx + 0:(space_steps-1)
  a_i <- dt/2*(s^2*(vetI*dx)^(2*gam)/(dx^2) - r*vetI)
  b_i <- -dt*(s^2*(vetI*dx)^(2*gam)/(dx^2) + r)
  c_i <- dt/2*(s^2*(vetI*dx)^(2*gam)/(dx^2) + r*vetI)
  
  Amatrix <- diag(1-b_i,space_steps)
  Amatrix[(row(Amatrix) - col(Amatrix)) == 1] <- -a_i[2:(space_steps)]
  Amatrix[(row(Amatrix) - col(Amatrix)) == -1] <- -c_i[1:(space_steps-1)]
  if(type=="cor"){
    V[,ncol(V)] <- 1 - ((vetS>lB | vetS<uB)*1)
  }
  if(type=="call"){
    V[,ncol(V)] <- pmax(vetS-K,0)
    V[nrow(V),] <- uB - K*exp(-r*T)
  }
  if(type=="compound"){
    V[,ncol(V)] <- pmax(other_asset-K,0)
    V[nrow(V),] <- other_asset[length(other_asset)] - K*exp(-r*seq(dt,T,dt)) #hardcoded strike for other asset for convenience
  }
  
  V[1,] = 0
  
  Bmatrix <- diag(1,space_steps)
  for (k in 1:time_steps-1){
    t = time_steps - k
    g <- c(rep(0,space_steps-1),c_i[length(c_i)]*(V[nrow(V),ncol(V)-1]))
    g[1] = a_i[1]*V[1]
    V[,t-1]<-solve(Amatrix, Bmatrix %*% V[,t] - g)
    
  }
  return(list(V=V,price=V[which(vetS==x0)]))
}

lB = 15
uB = 25
x0 = 20
T = 1/2
dt = 0.01
dx = 0.1
gam = 0.8
temp <- implicit_FD(T,dt,dx,gam,lB,uB,0,"cor",0,x0)

plot(seq(lB,uB,dx),temp$V[,1],type="l")



gam<-0.8
dt <- 0.01
T = 1/2
L1 <- 0
K_end = 20
L2 <- 3*K_end
K_compound<-2


options(scipen=999)
prices = c()
dx_size = c(0.1,0.05,0.02)
for(i in 1:length(dx_size)){
  dx = dx_size[i]
  v_call_End <-implicit_FD(1/4,dt,dx,gam,L1,L2,type="call",K_end,0,20)$V
  v_total <- implicit_FD(1/4,dt,dx,gam,L1,L2,type="compound",K_compound,v_call_End[,1],20)$price
  prices<-c(prices,v_total)
}

prices