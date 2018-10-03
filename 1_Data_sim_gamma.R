
#---- Define function to simulate capture-recapture data 
simul.js <- function(
  n.occasions,  # Number of primary seasons
  n.surveys,    # Number of secondary surveys within primary seasons
  gammaN,      # Uninfected entry probability 
  gammaI,      # Infected entry probability 
  N,            # Total abundance
  NN,          # Uninfected abundance
  NI,          # Infected abundance
  pN,          # Uninfected detection probability
  pI,          # Infected detection probability
  nsites       # Number of sites
){

########################################################   
#############   Table of Contents ######################       
########################################################  
  
  #---- This function has 5 parts
  # 1. Set up model parameters
  # 2. Assign true infection state
  # 3. Create observed capture history
  # 4. Cleaning up the data

########################################################   
######   1.  Set up parameters for the model ###########    
########################################################  
 
#---- 1. Define matrix for state process
    # 4-dimensional matrix
        # 1. Leaving state
        # 2. Arriving state
        # 3. Host abundance
        # 4. Number of transitions
PSI.state <- array(NA, dim = c(n.states, n.states, N, n.occasions-1))

for(i in 1:N){
  for(j in 1:n.occasions-1){
    PSI.state[,,i,j] <- matrix(c(1 - gammaN - gammaI, gammaN,  gammaI,  0,
                                 0,  phiN * (1- psi_NI),    phiN * psi_NI,   1-phiN,
                                 0,  phiI * psi_IN,    phiI * (1- psi_IN),   1-phiI,
                                 0,  0,        0,       1), 
                               nrow = 4, ncol = 4, byrow = T)
  }
}
  
# Set the number of sampling occasions  
  n.occasions <- dim(PSI.state)[4] + 1
  
# Generate no. of entering hosts between occasions
  B_U <- rmultinom(1, NN, rep(gammaN, times = n.occasions))
  B_I <- rmultinom(1, NI, rep(gammaI, times = n.occasions))

# Calculate total number of entering hosts each season
  B <- B_U + B_I

# Next step is to: scatter the number of individuals entering each season among the sites
# Empty vector to confirm that the number of individuals entering (B) is correct
entering_U <- numeric(length(B_U))  
# Matrix with the number of individuals per site that are entering
entering_sites <- matrix(NA, nrow = nsites, ncol = length(B_U)) 

# Fill in the matrix with the number of individuals per site that are entering each season
for(i in 1:n.occasions){          # For each occasion
  while (entering_U[i] != B_U[i]){  # While the number entering does not equal B_U
    entering_sites[,i] <- rpois(n = nsites, lambda = B_U[i]/nsites) # Generate a vector
    entering_U[i] <- sum(entering_sites[,i])
  }
}
# Create a vector from 1 to noccassions
n.occ <- 1:n.occasions  
# Create a vector from 1 to nsites
sit <- 1:nsites
# Create an empty vector to hold the occasion that each individual arrives in
ent.occ_U <- site.U <- numeric()    

for(i in 1:nsites){        # For each site
  for(j in 1:n.occasions){ # For each occassion
    # Occasion each individual is entering
    ent.occ_U <- c(ent.occ_U, rep(n.occ[j], each = entering_sites[i, j]))
    # Site ID that each individual is entering into
    site.U <- c(site.U, rep(sit[i], each = entering_sites[i, j]))
  }
}

# Copy everything from uninfected to infected individuals
ent.occ_I <- ent.occ_U
site.I <- site.U

# Combine together the uninfected and infected data
ent.occ <- c(ent.occ_U, ent.occ_I)
sites.all <- c(site.U, site.I)  

########################################################   
######   2.  Assign true infection state  ##############    
########################################################  
# Create empty matricies to store true and survey data 
CH.sur <- matrix(0, ncol = n.occasions, nrow = N)

#-- Each uninfected hosts enters the population as uninfected (= 2)
for (i in 1:length(ent.occ_U)){     # For each individual in uninfected
  CH.sur[i, ent.occ_U[i]] <- 2   # 2 = uninfected
}

#-- Each infected hosts enters the population as infected (= 3)
for(i in (length(ent.occ_U)+1):length(ent.occ)){    
  CH.sur[i, ent.occ[i]] <- 3   # 3 = infected
}

#-- Now determine what happens to that individual after it enters the population
  # Does it survive?
  # Does it transistion disease states?
for (i in 1:N){    # For each individual in the super population       
  # If the entry occasion = the last occasion, go next
  if (ent.occ[i] == n.occasions) next   
  # For each occasion from when an individuals enters all the way to the the last sampling occasion    
    for(t in (ent.occ[i]+1):n.occasions){  
      # Determine individual state history
      sur <- which(rmultinom(1, 1, PSI.state[CH.sur[i, t-1], , i, t-1]) == 1)
      CH.sur[i, t] <- sur
    } #t
} #i

# Calculate the total number of individuals in the matricies created
Nt <- numeric(n.occasions)
for(i in 1:n.occasions){
  Nt[i] <- length(which(CH.sur[,i] > 0))
}
 
########################################################   
######   3.  Create observed capture history  ##########   
########################################################  

# Calculate the probability of misdiagnosing individuals
# This occurs when the swab or qPCR sample results in 0 zoospores
# You would use ySp and Spore to create the misclassification probability

#---- Define matrix for observation process
  # 4-dimensional matrix
    # 1. true state
    # 2. Observation state
    # 3. Host abundance
    # 4. Number of occasions
PSI.obs <- array(NA, dim = c(n.states, n.obs, N, n.occasions, n.surveys))

for(i in 1:N){
  for(j in 1:(n.occasions)){
    for(t in 1:n.surveys){
    PSI.obs[,,i,j,t] <- matrix(c(0,       0,   1,
                              pN,      0,   1-pN,
                              0,      pI,   1-pI,
                              0,        0,   1), 
                             nrow = 4, ncol = 3, byrow = T)
    }
  }
}

#-- Create empty matrix to hold the data
# Dimensions of data
  # row = individual
  # column = primary period
  # sheets = secondary surveys
CH.p <- array(0, dim = c(N, n.occasions, n.surveys))


#-- 
for(i in 1:N){    # For each individual in the super population       
  # For each occasion from when an individuals enters to the last occasion    
    for(j in ent.occ[i]:n.occasions){  
        for(t in 1:n.surveys){  # For each survey
      # Determine if the individual is seen or not
        event <- which(rmultinom(1, 1, PSI.obs[CH.sur[i, j], , i, j, t]) == 1)
        CH.p[i, j, t] <- event
        } #ts
    } #js
} #is

# Note the numbers are 

# Copy the original data frame
CH.p2 <- CH.p



########################################################   
######   4.  Cleaning up the data  ##################### 
########################################################  

# Dimensions of data
  # row = individual
  # column = primary period
  # sheets = secondary surveys

##########
# Remove individuals never captured
##########

# Create an empty vector to save the individuals very captured in the population
never_cap <- numeric()

for(i in 1:dim(CH.p)[1]){   # For each indiviudal
  for(j in 1:dim(CH.p)[3]){  # For each survey
  
  # This only works if the entire capture history is == to 3 or ==0 (i.e., dead or not seen)
  n_seen <- grep("[03]", CH.p[i,,j])
  # If there were individuals not seen, in each 
  if(length(n_seen) == dim(CH.p)[2]){ # If the length of n_seen == the total number of sampling occasions
    never_cap <- c(never_cap, i) # then add the individuals ID # to the vector never_cap
  }
  }
}

never_cap <- unique(never_cap)  

# Remove individuals never captured from the observation matrix
if(length(never_cap) > 0){  # If the length of never_cap is greater than 0 then
  CH.p <- CH.p[-never_cap,,]   # Remove them from these matricies
  CH.sur2 <- CH.sur[-never_cap,]
}

# Individuals never captured == 1 in standard notation for Jolly-Seber model
CH.sur[CH.sur == 0] <- 1
CH.sur2[CH.sur2 == 0] <- 1


########################################################   
######   5.  Format data for N-mix mod ################# 
########################################################  

# Create a new matrix with the counts of uninfected & infected hosts at each site
  # each matrix = nsites by n.occasions
# TRUE Uninfected host matrix- with true abundance
NNT <- matrix(NA, ncol = n.occasions, nrow = nsites)
# TRUE Infected host matrix- with true abundance
NIT <- matrix(NA, ncol = n.occasions, nrow = nsites)
# Observed uninfected host abundance
NNobs <- array(NA, dim = c(nsites, n.surveys, n.occasions))
# Observed infected host abundance
NIobs <- array(NA, dim = c(nsites, n.surveys, n.occasions))
# The numbers in CH.sur mean:
# 4 true states
  # 1 = Not entered
  # 2 = uninfected
  # 3 = infected
  # 4 = Dead
# The numbers in CH.p2 mean:
# 3 observed states
  # 1 = seen uninfected
  # 2 = seen infected
  # 3 = not seen

for(i in 1:nsites){ # For each site
  RR <- which(sites.all == i)   # Subset all the rows with site == i
    for(j in 1:n.occasions){ # For each sampling occasion
    NNT[i, j] <- length(which(CH.sur[RR, j] == 2))  # Add up all the individuals that were uninfected
      # Sum up the number of rows 
      # Not infected in True matrix
    NIT[i, j] <- length(which(CH.sur[RR, j] == 3)) # Add up all the individuals that were infected
      # Infected in True matrix
    for(t in 1:n.surveys){  # For each survey
      NNobs[i, t, j] <- length(which(CH.p2[RR, j, t] == 1))
      #Add up all the individuals that were seen uninfected
      # 1 = seen uninfected
      NIobs[i, t, j] <- length(which(CH.p2[RR, j, t] == 2))
      # Add up all the individuals that were seen uninfected
      # 2 = seen infected
    }
  }
}

return(list(# Data for CMR model
            CH.p = CH.p, 
            CH.p2 = CH.p2, 
            CH.sur = CH.sur, 
            CH.sur2 = CH.sur2,
            B = B,
            Nt = Nt,
            
            # Data for N-mixture model
            NNobs = NNobs,
            NIobs = NIobs,
            NNT = NNT,
            NIT = NIT
            ))
}

