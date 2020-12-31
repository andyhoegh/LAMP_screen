#### SEIR for waiting groups within true SEIR model


sir_simple_step <- function(waiting.group, sims,
                            I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                            theta, gamma_I1I2, gamma_I2R,
                            beta_vec.on, beta_vec.off, delta.t = 1) {
  
  dN_SE.on <- rbinom(n=sims,size=waiting.group[,1],prob=1-exp(-beta_vec.on*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) # add random introductions
  dN_EI1.on <- rbinom(n=sims,size=waiting.group[,2],prob=1-exp(-theta*delta.t))
  dN_I1I2.on <- rbinom(n=sims,size=waiting.group[,3],prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.on <- rbinom(n=sims,size=waiting.group[,4],prob=1-exp(-gamma_I2R*delta.t))
  dN_SE.off <- rbinom(n=sims,size=waiting.group[,6],prob=1-exp(-beta_vec.off*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) # add random introductions
  dN_EI1.off <- rbinom(n=sims,size=waiting.group[,7],prob=1-exp(-theta*delta.t))
  dN_I1I2.off <- rbinom(n=sims,size=waiting.group[,8],prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.off <- rbinom(n=sims,size=waiting.group[,9],prob=1-exp(-gamma_I2R*delta.t))
  
  # update classes
  waiting.group[,1] <- waiting.group[,1] - dN_SE.on 
  waiting.group[,2] <- waiting.group[,2] + dN_SE.on - dN_EI1.on 
  waiting.group[,3] <- waiting.group[,3] + dN_EI1.on - dN_I1I2.on
  waiting.group[,4] <- waiting.group[,4] + dN_I1I2.on - dN_I2R.on
  waiting.group[,5] <- waiting.group[,5] + dN_I2R.on
  
  waiting.group[,6] <- waiting.group[,6] - dN_SE.off 
  waiting.group[,7] <- waiting.group[,7] + dN_SE.off - dN_EI1.off 
  waiting.group[,8] <- waiting.group[,8] + dN_EI1.off - dN_I1I2.off
  waiting.group[,9] <- waiting.group[,9] + dN_I1I2.off - dN_I2R.off
  waiting.group[,10] <- waiting.group[,10] + dN_I2R.off
  return(waiting.group)
}