#This file contains the model structure used for the HMM on the Adelie movement data.

#Stan requires the model to be stored in an external ".stan" file which we draw back into R when we want to run the model. The following code creates a temporal variable to store the model and then exports it to a ".stan" file wherever the current file is located

#Stan requires the model to be broken up into "blocks" separated by curly brackets "{}", which dictate different elements of the model. Each line of code is broken up with a semicolon ";".

#Note that in R comments are denoted by a preceeding hashtag "#" while in Stan double forward slashes "//" are used.

t_m = '

// This block allows us to define any functions that we may need in our model. As Stan does not have the wrapped Cauchy distribution native to it we need to put the log pdf in manually so that we can call it later.

functions {
  real wc_lpdf(real y, real mu, real rho){
    return log(1 / (2 * pi()) * (1 - pow(rho,2)) / (1 + pow(rho,2) - 2 * rho * cos(y - mu)));
  }
}

// This block allows us to identify what elements of the model should be considered data. It allows us to define the data type ("int" = integer, "real" = continuous), any constraints in angle brackets "<>", as well as the length of the data in square brackets "[]".

data {
  int<lower=1> K;                       // The number of underlying states
  
  int<lower=1> N;                       // The number of observations across individuals
  int<lower=1> I;                       // Number of individuals
  
  int<lower=0> y1[N];                   // The number of dives
  
  real<lower=0> y2[N];                  // The step lengths
  
  real<lower=0, upper=2*pi()> y3_0[I];  // The first turn angle for each individual
  real<lower=0, upper=2*pi()> y3[N];    // Each subsequent turn angle

  real<lower=0> y4[N];                  // The duration of each time step
  
  int s[I];                             // The number of observations per individual
  int z[I];                             // The accumulated number of observations. This allows us to keep track of which observation in the vectors storing the data does each individual start.
}

// This block allows us to identify what elements of the model should be considered parameters. Similar to the data section we can define the parameter type (simplex = all row elements add to 1), constraints in angle brackets "<>", and the number of parameters using square brackets "[]". 

// Where multiple sets of square brackets or multiple numbers/letters within the square are used, indicates that instead of a vector, we are instead creating a matrix to store the estimated parameter values. 


parameters {
  simplex[K] delta[I];                  // An I x K matrix that holds the initial distribution
  simplex[K] G[I,K];                    // An I x K x K matrix that holds the transition probabilities
  
  vector<lower=0>[K] d_lam[I];          // Lambda for Poisson distributions on the number of dives
  real d_lmu[K];                        // Mu for the log normal distributions on Lambda
  real<lower=0> d_lsig[K];              // Sigma for the log normal distributions on Lambda
  
  vector[K] s_mu[I];                    // Mu for log normal distributions on step length
  real s_mmu[K];                        // Mu for the normal distributions on Mu
  real<lower=0> s_msig[K];              // Sigma for the normal distributions on Mu
  
  vector<lower=0>[K] s_sig[I];          // Sigma for log normal distributions on step length
  
  vector<lower=0, upper=1>[K] t_rho[I]; // Rho for wrapped Cauchy distributions on turning angle
  real<lower=0> t_ralp[K];              // Alpha for the beta distributions on rho
  real<lower=0> t_rbet[K];              // Beta for the beta distributions on rho
}

// This block allows us to describe the model structure after establishing all the elements within it in previous blocks.

model {

  // Defines a position marker that allows us to keep track of which observation each individual should start at. The first starts at 1, and then after we run through the model for that first individual, "pos" is updated to be equal to the first observation for the second individual and so on.
  
  int pos;
  pos = 1;
  
  // Creates a for loop that will run the model for each individual iteratively.

  for(i in 1:I){
    
    // Because of how Stan functions we cannot draw a subset of the data out within our model structure and instead need to define new data elements that will temporarily store the subsetted data for us. Also because of Stan percularities we need to define each data element before we tell Stan what to fill it with.  
  
    int y1i[s[i]];  // The number of dives for individual "i"
    real y2i[s[i]]; // The step lengths for individual "i"
    real y3i[s[i]]; // The turning angles for individual "i"
    real y4i[s[i]]; // The duration of the time steps for individual "i"
      
    // We also need to define some parameters that are unconcern with their specific estimation but are required by the forward algorithm
      
    vector[K] logalpha[s[i]]; // Stores the log of the forward/alpha values in the forward algorithm. These are log-likelihood estimations that represent the likelihood over all possible state-to-state paths through the HMM
    real accumulator[K];      // Used as a temporary storage for each alpha value following the first.
      
    // Define what data is stored in each of the temporary data variable.
      
    y1i = y1[pos:z[i]];
    y2i = y2[pos:z[i]];
    y3i = y3[pos:z[i]];
    y4i = y4[pos:z[i]];
    
    // Where the forward algorithm is implemented. For a more detailed explanation see the accompanying documents.
      
    // For each state ("K") we need to calculate the first log-alpha value based on the first observed data points. Following this we can calculate the log-alpha value for each subsequent observation iteratively. The initial log-alpha values are the following summation of initial distribution and the log pdf/pmf of the state-dependent distributions on our data.
      
    for (h in 1:K){
      logalpha[1,h] = log(delta[i,h]) +                                 // Initial distribution
                      poisson_lpmf(y1i[1] | d_lam[i,h] * y4i[1]) +      // Poisson on number of dives
                      lognormal_lpdf(y2i[1] | s_mu[i,h], s_sig[i,h]) +  // Log normal on step lengths
                      wc_lpdf(y3i[1] | y3_0[i], t_rho[i,h]);            // Wrapped Cauchy on turning angle
    }
    
    // After the initial log-alpha values are calculated we can iteratively calculate the log-alpha values for each subsequent observations for that individual (from 2 to s[i]).
    
    for (t in 2:s[i]) {
    
    // At each observation we consider every state-to-state ("j" being the current state and "l" being the previous) transition possible. At each current state "j" we calculate a precursor log-alpha value (denoted as accumulator) given each possible previous state "l". This can be considered a probability of transitioning from state "i" to state "j" This value is based on the following summation of the previous log-alpha value for that state, the transition probability from "l" to "j", and the log pdf/pmf of the state-dependent distributions on our data.
    
      for (j in 1:K) {
        for (l in 1:K) {
          accumulator[l] =  logalpha[t-1, l] +                                // Previous log-alpha values
                            log(G[i, l, j]) +                                 // Transition probability from state "i" to "j"
                            poisson_lpmf(y1i[t] | d_lam[i,j] * y4i[t]) +      // Poisson on number of dives
                            lognormal_lpdf(y2i[t] | s_mu[i,j], s_sig[i,j]) +  // Log Normal on step length
                            wc_lpdf(y3i[t] | y3i[t-1], t_rho[i,j]);           // Wrapped Cauchy on turning angle
        }
      
        // To get the log-alpha values for state "j" at time "t" we take the logsumexp of "accumulator" or the probabilities of transitioning from any state "i" at time "t-1" to state "j" at time "t". This results in log-alpha storing log-likelhood values that represent the likelihood over all possible paths through the HMM. We use the logsumexp function to avoid numerical underflow as the probabilities we deal we can reduce to very small values.
      
        logalpha[t, j] = log_sum_exp(accumulator);
      
      }
    }
  
    // To calculate the overall log-likelihood we compute the logsumexp over the final log-likelihoods store in log-alpha. We once again use the logsumexp function to avoid numerical underflow.    
  
    target += log_sum_exp(logalpha[s[i]]);
  
    // Establish the priors for our state-dependent parameters and our population level parameters.
    
    for(k in 1:K){
      d_lam[i,k] ~ lognormal(d_lmu[k], d_lsig[k]);  // Lambda from the Poisson distributions on number of dives
      d_lmu[k] ~ normal(0, 10);                     // Mu from from the log normal distributions on Lambda
      d_lsig[k] ~ gamma(5, 1);                      // Sigma from log normal distributions on Lambda
    
      s_mu[i,k] ~ normal(s_mmu[k], s_msig[k]);      // Mu from the log normal distributions on the step lengths
      s_mmu[k] ~ normal(0, 5);                      // Mu from the normal distributions on Mu
      s_msig[k] ~ gamma(5, 1);                      // Sigma from the normal distributions on Mu
    
      t_rho[i,k] ~ beta(t_ralp[k], t_rbet[k]);      // Rho from the wrapped Cauchy distributions on turning angle
      t_ralp[k] ~ gamma(5, 1);                      // Alpha from the beta distributions on Rho
      t_rbet[k] ~ gamma(5, 1);                      // Beta from the beta distributions on Rho
    }
    s_sig[i,] ~ student_t(3,0,1);                   // Sigme from the log normal distributions on the step lengths
    
    // Update our position marker for the next individual.
    
    pos = pos + s[i];
  }
}
'
#Set working directory to the current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Sends the model to a ".stan" file
cat(t_m, file = "mod.stan")

#Removes the model from the R environment
remove(t_m)
