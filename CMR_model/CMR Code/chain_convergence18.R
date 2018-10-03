#------ call Library
library("jagsUI")
library("rjags")
library("coda")

# Call the function to simulate the data
source('1_Data_sim_gamma.R')

# This code simulates the data for a Jolly-Seber model and then aggregates the data for an N-mixture model

#--- Survey conditions
n.occasions <- 10   # Number of primary sampling occasions
n.surveys   <- 3     # Secondary surveys

#---- Population parameters
# Super population size
NN <- 300  # Number uninfected
NI <- 300  # Number infected

N <- NN + NI  # Total population size

# Number of true states
n.states <- 4

# Number of observed states
n.obs <- 3

#--- Define parameter values
# Uninfected survival probability
phiN <- 0.9

# Infected survival probability
phiI <- 0.8

# Entry probability     
# Must sum to one across all sampling ocassions
# With n.occasions, n.occasions entry probability must be defined
gammaN <- (1/n.occasions)/4
gammaI <- (1/n.occasions)/4

# Transition probability 
psi_NI <- 0.2 # Going from Uninfected to Infected
psi_IN <- 0.2 # Going from Infected to Uninfected

# Detection probability
pN <- 0.9  # Uninfected host
pI <- 0.9  # Infected host

# Number of sites
nsites <- 500

# Execute simulation function
sim <- simul.js(n.occasions = n.occasions, 
                n.surveys = n.surveys,
                gammaN = gammaN, 
                gammaI = gammaI, 
                N = N,
                pN = pN, 
                pI = pI,
                NN = NN, 
                NI = NI,
                nsites = nsites
)

####################################
######### Bundle the data ##########
####################################

# Analysis of the JS model as a multistate model
CH <- sim$CH.p

# Augment data
nz <- 500

CH.ms <- array(0, dim = c(dim(CH)[1]+nz, dim(CH)[2]+1, dim(CH)[3]))

for(i in 1:dim(CH)[1]){
  for(j in 1:dim(CH)[2]){
    for(t in 1:dim(CH)[3]){
      CH.ms[i,j+1,t] <- CH[i, j, t]
    }
  }
}

CH.du <- CH.ms[1:dim(CH)[1],,]

# Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 3        # Not seen = 3, seen = 1 or 2

# Bundle data
jags.data <- list(y = CH.ms, 
                  M = dim(CH.ms)[1],
                  K = dim(CH.ms)[2],
                  n.surv = dim(CH.ms)[3],
                  state = 3
)

# Initial values
ch <- CH.du

js.multistate.init <- function(ch, nz){
# 4 state system
    # 1 = Not entered
    # 2 = uninfected
    # 3 = infected
    # 4 = Dead 
  
# 3 observation states
    # seen uninfected
    # seen infected
    # not seen

# Put an NA when an individual was not seen ( = 0 or 3)
  ch[ch==0] <- NA
  ch[ch==3] <- NA
  
  state <- ch
  colnames(state) <- 1:dim(state)[2]
  
  # When the individual is known to be alive between first and last capture fill it in with 2s
  for(i in 1:dim(ch)[1]){    # For each individual
    for(j in 1:dim(ch)[3]){
    n1 <- min(which(ch[i,,j] < 3))
    n2 <- max(which(ch[i,,j] < 3))
    
    fill <- which(is.na(state[i,n1:n2,j]) == TRUE)
    state[i, names(fill),j] <- 2
    }
  }
  state <- state + 1
  
  f <- array(NA, dim = c(dim(ch)[1], dim(ch)[3]))
  
    for(i in 1:dim(ch)[1]){
      for(j in 1:dim(ch)[3]){
        f[i, j] <- min(which(!is.na(ch[i, , j])))
      }
    }

  l <- array(NA, dim = c(dim(ch)[1], dim(ch)[3]))
  
  for(i in 1:dim(ch)[1]){
    for(j in 1:dim(ch)[3]){
      l[i, j] <- max(which(!is.na(ch[i, , j])))
    }
  }

  for (i in 1:dim(ch)[1]){
    for(j in 1:dim(ch)[3]){
# Before initial observation- not observed
    state[i,1:(f[i,j]-1),j] <- 1
    
      # If the last time the animal was seen != the last survey date
      # Then add 1 to the survey date, and fill the rest to the last occasion with 4
      if(l[i,j]!= dim(ch)[2]){state[i, (l[i,j]+1):dim(ch)[2], j] <- 4}
    
      if(ch[i, f[i,j], j] == 1){state[i, f[i,j], j] <- 2}
      if(ch[i, f[i,j], j] == 2){state[i, f[i,j], j] <- 3} 
    }
  }   
 
  state2 <- array(NA, dim = c(dim(ch)[1]+nz, dim(ch)[2], dim(ch)[3]))
  state3 <- array(NA, dim = c(dim(ch)[1]+nz, dim(ch)[2]))

# Copy over the matrix you generated to the one with the right dimensions    
  for(i in 1:dim(ch)[1]){
    for(j in 1:(dim(ch)[2])){
      for(t in 1:dim(ch)[3]){
        state2[i,j,t] <- state[i, j, t]
      }
    }
  }    

  # For all the data augmented individuals- their state == 1
  for(i in (dim(state)[1]+1):dim(state2)[1]){ 
    for(j in 1:(dim(state)[2])){ # For all occasions
      for(t in 1:dim(state)[3]){ # For all surveys
        state2[i,j,t] <- 1
      }
    }
  }  
  
#  # Now colapse the 3D matrix to 2D for the initial values of z
#  for(i in 1:dim(state2)[1]){
#    for(j in 1:dim(state2)[3]){
#      state3[i, j] <- sample(state2[i, j, ], 1)
#    }
#  }
#  
#  # Make the first column == NA
#  state3[,1] <- NA
  
  return(state2)
}

n.occasions <- dim(CH.ms)[2]

zinit <- js.multistate.init(CH.du, nz)

dim(zinit)
zinit <- apply(zinit, c(1, 2), max)
zinit[,1] <- NA
dim(zinit)


zz <- rbind(sim$CH.sur2, matrix(1, ncol = ncol(sim$CH.sur2), nrow = nz))
zz <- cbind(rep(NA, times = nrow(zz)), zz)

dim(zz)
dim(CH.ms)
#--- Bundle the inits 

inits <- function(){list(phiN = phiN, 
                         phiI = phiI, 
                         
                         pN = pN, 
                         pI = pI, 
                         
                         psi_NI = psi_NI,
                         psi_IN = psi_IN,
                         
                         gammaN = rep(gammaN, times = 1), 
                         gammaI = rep(gammaI, times = 1),
                         
                         z = zz
                         )}    

#------- Parameters monitored

params <- c("phiN", 
            "phiI", 
            "pN", 
            "pI",
            "psi_NI",
            "psi_IN",
            "gammaN",
            "gammaI",
            "NN", "NI", 
            "Ntot", "Nsuper"
            )
#------- MCMC settings
nc <- 18
ni <- 3278
na <- 10000
nb <- 1000
nt <- 1

#------- Call JAGS from R
cat("Starting JAGS...\n")
start_time <- Sys.time()
js.ms <- jags(data = jags.data, 
              inits = inits, 
              parameters.to.save = params, 
              model.file = "JS.txt", 
              n.chains = nc, 
              n.thin = nt, 
              n.iter = ni, 
              n.burnin = nb, 
              n.adapt = na,
              parallel = TRUE)
cat("Completed JAGS...\n")
print(Sys.time() - start_time)

save(js.ms, file = paste("js_", nc, ".rda", sep = ""))



