explicit_FD_put <- function(dt, ds, dv, K){
  
  # Hardcoded parameters
  r <- 0.05
  kappa <-1
  theta <- .2
  eta <- .5
  rho <- -0.4
  T <- 1
  S_max <- 1.5*K
  V_min <- 0.05
  V_max <- 0.6
  
  
  time_steps <- as.integer(T/dt)
  s_steps <- round(S_max/ds)+1
  v_steps <- round((V_max - V_min)/dv)+1
  S <- ds *(0:(s_steps-1))
  V <- V_min + dv *(0:(v_steps -1))
  
  
  #Terminal Condition
  term_vals <- pmax(K- S,0) 
  f <- matrix(rep(term_vals,v_steps),nrow = v_steps, ncol = s_steps,byrow = TRUE)
  rownames(f) <- V; colnames(f) <- S
  
  for(i in 1:time_steps){
    f_upd <- matrix(0,ncol = s_steps,nrow=v_steps)
    colnames(f_upd)<-S
    rownames(f_upd)<-V
    for(j_ind in 2:(v_steps-1)){ #Exclude first and last S value
      for( k_ind in 2:(s_steps-1)){  #Exclude first and last v value
        j = V[j_ind]; k = S[k_ind]
        fjk <- 1 - r*dt - j*k^2 *dt *dv - eta^2 * j *(dt/dv)
        fj_km1 <- (dt/2)*(j*k^2*dv - r*k)
        fj_kp1 <- (dt/2)*(r*k + j*k^2*dv)
        fjp1_k <- (dt/(2*dv))*(eta^2*j + kappa*(theta - j*dv))
        fjm1_k <- (dt/(2*dv))*(eta^2*j - kappa*(theta - j*dv))
        cross_same <- rho * eta * j *k *dt
        cross_different <- - cross_same
        f_upd[j_ind,k_ind] <- f[j_ind,k_ind]*fjk +f[j_ind,k_ind-1]*fj_km1 +
          f[j_ind,k_ind+1]*fj_kp1 + f[j_ind+1,k_ind]*fjp1_k + 
          f[j_ind-1,k_ind]*fjm1_k + 
          cross_same*(f[j_ind+1,k_ind+1]+ f[j_ind-1,k_ind-1])+
          cross_different*(f[j_ind-1,k_ind+1]+ f[j_ind+1,k_ind-1])
      }
    }
    
    # Enforce Boundary Conditions
    # vol upper
    f_upd[1,2:(s_steps-1)] <- 2*f_upd[2,2:(s_steps-1)]- f_upd[3,2:(s_steps-1)]
    # vol lower
    #f_upd[v_steps,2:(s_steps-1)] <- 2*f_upd[v_steps-1,2:(s_steps-1)]- 
          f_upd[v_steps-2,2:(s_steps-1)]
    
    f_upd>
    # S lower 
    f_upd[,1] <- K*exp(-r * i *dt)
    
    # S Upper - should be initialized to 0, but reset here to be safe
    f_upd[,s_steps]<-0
    #update F
    f<-f_upd[1,2:(s_steps-1)]<- -r*f[1,2:(s_steps-1)]
  }
  return(list(put_price = f,underlying_price = S, vol = V))
}

dt = 0.0001
dv = 0.01
ds = 5
s0 = 100
v0 = 0.25
out<-explicit_FD_put(dt,ds,dv,K)
v_out <- which(out$vol==v0)
s_out<- which(out$underlying_price==s0)
print(out$put_price[s_out,v_out])