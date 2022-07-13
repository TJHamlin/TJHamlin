#Loads required packages
library(rstan)
library(bayesplot)
library(circular)

#Function that simulates a Dirichlet distribution that we need to simulate some initial values
rdirch = function(n,k,alpha=rep(1,k)){
  out = matrix(NA,n,k)
  for(i in 1:n){
    tmp = rgamma(k,alpha,1)
    out[i,] = tmp/sum(tmp)
  }
  return(out)
}

#Set where the data is located within your file directory
data_loc <- ""

#Set the individuals of interest. Alternatively, include all individuals.
ids <- c("N415A1", "N415A2", "N430A1", "N430A2", "N442A1", "N442A2", "N446A1", "N446A2")
#ids <- gsub(".csv", "", list.files(data_loc))
inds <- length(ids)

#Number of unknown states in the model
sts <- 3

#Number of iterations the model should run for
itrs <- 1000

#Function that draws in all of the data and formats it
data_fun <- function(){
  id <- NULL
  obs <- NULL
  dives <- NULL
  time <- NULL
  speed <- NULL
  t0 <- NULL
  turn <- NULL
  for(i in 1:inds){
    tdata <- read.csv(paste0(data_loc, "/", ids[i], ".csv"))
    
    tdata$speed[which(tdata$speed==0)] <- 0.0001
    tdata$time <- tdata$time/3600
    
    id <- c(id, rep(ids[i], nrow(tdata)-2))
    obs <- c(obs, nrow(tdata)-2)
    dives <- c(dives, tdata$dives[3:nrow(tdata)])
    time <- c(time, tdata$time[3:nrow(tdata)])
    speed <- c(speed, tdata$speed[3:nrow(tdata)])
    t0 <- c(t0, tdata$ta[2])
    turn <- c(turn, tdata$ta[3:nrow(tdata)])
  }
  accm <- NULL
  accm[1] <- obs[1]
  for(i in 2: length(obs)){
    accm[i] <- accm[i-1] + obs[i]
  }
   
  return(list("id"=id, "obs"=obs, "accm"=accm, "dives"=dives, "time"=time, "speed"=speed, "turn"=turn, "t0"=t0))
}
data <- data_fun()

#Function that sets all the initial values for the model
init_fun <- function() {
  #Lambda for the Poisson distributions on the number of dives
  i_dlam = matrix(rep(abs(rnorm(sts, c(0, 7, 20), 0.1)),inds), nrow=inds, ncol=sts, byrow=T, dimnames=NULL)
  #Mu for the log normal distributions on Lambda
  i_dlmu = abs(rnorm(sts, c(0.5, 2, 3), 0.1))
  #Sigma for the log normal distributions on Lambda
  i_dlsig = abs(rnorm(sts, c(0.2, 0.2, 0.2), 0.1))
  #
  #Mu for the log normal distributions on the step lengths
  i_mu = matrix(rep(rnorm(sts, c(0, 1, 2), 0.1),inds), nrow=inds, ncol=sts, byrow=T, dimnames=NULL)
  #Mu for the normal distributions on Mu
  i_nmu = rnorm(sts, 0, 1)
  #Sigma for the normal distributions on MU
  i_nsig = abs(rnorm(sts, 0, 1))
  #Sigma for the log normal distributions on the step lengths
  i_sig = matrix(rep(abs(rnorm(sts, c(0.5, 0.5, 0.5), 0.1)),inds), nrow=inds, ncol=sts, byrow=T, dimnames=NULL)
  #
  #Rho for the wrapped Cauchy distributions on the turning angles
  i_rho = matrix(rep(rbeta(sts, 2, 5),inds), nrow=inds, ncol=sts, byrow=T, dimnames=NULL)
  #Alpha for the beta distributions on Rho
  i_balp = abs(rnorm(sts, 0, 1))
  #Beta for the beta distributions on Rho
  i_bbet = abs(rnorm(sts, 0, 1))
  #
  #Initial probability distribution
  i_delta = rdirch(inds,sts)
  #Transition probability matrix
  i_G <- array(NA, c(inds,sts,sts))
  for(i in 1:inds){
    i_G[i,,] <- rdirch(sts,sts)
  }
  #
  list(d_lam = i_dlam, d_lmu = i_dlmu, d_lsig = i_dlsig, 
       s_mu = i_mu, s_sig = i_sig, s_mmu = i_nmu, s_msig = i_nsig,
       t_rho=i_rho, t_ralp = i_balp, t_rbet = i_bbet, 
       delta = i_delta, G = i_G
       )
} 


#Runs the HMM through stan, see ?stan for a description of the inputs
mod_stan = stan(file = "mod.stan",
                data = list(K = sts,
                            N = length(data$id),
                            I = inds,
                            y1 = data$dives,
                            y2 = data$speed,
                            y3 = data$turn,
                            y3_0 = data$t0,
                            y4 = data$time,
                            s = data$obs,
                            z = data$accm),
                init = init_fun,
                iter = itrs,
                warmup = itrs/5,
                cores = 4,
                pars = names(init_fun()),
                save_warmup = T
)

#Creates an output table for the HMM
print(mod_stan, pars=c("logalpha", "accumulator", "delta"), include=F, probs=c(0.025, 0.5, 0.975))
