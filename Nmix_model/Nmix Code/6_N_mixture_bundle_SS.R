#------ call Library
library("jagsUI")
library("rjags")
library("coda")

# Call the function to simulate the data
source('1_Data_sim.R')

# This code simulates the data for a Jolly-Seber model and then aggregates the data for an N-mixture model

#--- Survey conditions
n.occasions <- 10   # Number of primary sampling occasions
n.surveys   <- 3     # Secondary surveys

#---- Population parameters
# Super population size
NN <- 500  # Number uninfected
NI <- 500  # Number infected

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


###########################################
#--- Now bundle the data for the N-mixture model
###########################################

# Bundle data
win.data <- list(yN = sim$NNobs, 
                 yI = sim$NIobs,
                 
                 R = dim(sim$NNobs)[1], 
                 T = dim(sim$NNobs)[2],
                 K = dim(sim$NNobs)[3])

# Initial values
R = dim(sim$NNobs)[1]
T = dim(sim$NNobs)[2]
K = dim(sim$NNobs)[3]

NIst <- apply(sim$NIobs, c(1, 3), max, na.rm = TRUE)
NNst <- apply(sim$NNobs, c(1, 3), max, na.rm = TRUE)

SI <- GI <- array(NA, c(R, K-1))
SN <- GN <- array(NA, c(R, K-1))
TNn <- TIi <- array(NA, c(R, K-1))

SI[] <- 1
SN[] <- 1

TIi[] <- 0
TNn[] <- 0

SI[which(NIst[,-n.occasions] == 0)] <- 0
SN[which(NNst[,-n.occasions] == 0)] <- 0

GI <- NIst[,-1] - SI 
GI[GI == -1] <- 0
GN <- NNst[,-1] - SN 
GN[GN == -1] <- 0


NNst[ , -1] <- NA
NIst[ , -1] <- NA

inits <- function() {list( 
  NN = NNst, 
  NI = NIst, 
  
  GI = GI,  
  GN = GN, 
  
  TI = TIi, 
  TN = TNn, 
  
  SI = SI,
  SN = SN,
  
  alpha.lamN = runif(1, -1, 1),
  pN = qlogis(runif(1, 0.7, 1)),
  phiN = qlogis(runif(1, 0.7, 1)),
  gammaN = runif(1, -1, 1),
  psi_NI = qlogis(runif(1, 0.1, 0.4)),
  
  alpha.lamI = runif(1, -1, 1),
  pI = qlogis(runif(1, 0.7, 1)),
  phiI = qlogis(runif(1, 0.7, 1)),
  gammaI = runif(1, -1, 1),                        
  psi_IN = qlogis(runif(1, 0.1, 0.4)) 
  
)}

# Monitor Parameters
params <- c("alpha.lamN", "alpha.lamI",
            "pN", "pI",
            "phiN", "phiI",
            "gammaN", "gammaI",
            "psi_NI", "psi_IN"

)

# MCMC settings
ni <- 9833
nb <- 1000
na <- 10000
nt <- 1
nc <- 6


library("jagsUI")

cat("Starting JAGS...\n")
start_time <- Sys.time()
out <- jags(win.data, inits, params, "Nmix.txt", 
            n.chains = nc, 
            n.thin = nt, 
            n.iter = ni, 
            n.adapt = na,
            n.burnin = nb, 
            parallel = TRUE)
cat("Completed JAGS...\n")
print(Sys.time() - start_time)


rhat <- max(unlist(out$Rhat))

counter <- 0

while (rhat > 1.1){
  
  counter <- counter + 1
  
    if (counter > 10){break}
  
    out <- update(out, 
                  n.iter = 50000)
    
    rhat <- max(unlist(out$Rhat))
    
}

save(out, file = paste("Nmix_", nc, ".rda", sep = ""))

