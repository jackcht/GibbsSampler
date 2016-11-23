# Machine Learning Assignment 7
# Problem 4

# Gibbs Sampler for the Hierachical Normal Model
#   with Gelman-Rubin diagnostic 

require(coda)

###################
## Gibbs Sampler ##
###################
#
#
# j groups (from j = 1 to J(=4), Diet A/B/C/D)
# i (from 1 to nj) measurements for each group (n1 = 4, n2 = 6, n3 = 6, n4 = 8)
# Diet A: 62, 60, 63, 59
# Diet B: 63, 67, 71, 64, 65, 66
# Diet C: 68, 66, 71, 67, 68, 68
# Diet D: 56, 62, 60, 61, 63, 64, 63, 59
y <- matrix(0, nrow = 8, ncol = 4)
y[,1] <- c(62, 60, 63, 59, 0, 0, 0, 0)
y[,2] <- c(63, 67, 71, 64, 65, 66, 0, 0)
y[,3] <- c(68, 66, 71, 67, 68, 68, 0, 0)
y[,4] <- c(56, 62, 60, 61, 63, 64, 63, 59)

n <- c(4,6,6,8)
N <- sum(colSums(y != 0))

# Each data: y_ij | theta_j ~ N(theta_j, sigma^2)       ~(mean: theta_j, common var: sigma^2)
# Group means: theta_j ~ N(mu, tau^2)
# Uniform prior for (mu, log(sigma), tau) 
# --> Posterior:    p(theta, mu, log(sigma), log(tau)|y) 
#                       -PROPORTIONAL_TO- 
#                   tau * ProdSum(N(theta_j|mu, tau^2) * ProdSum(N(y_ij|theta_j, sigma^2))

# starting points of theta_j & mu
theta.start <- colSums(y)/n         # here choose the mean as starting point
mu.start <- mean(theta.start)


# conditional posterior for theta_j
theta.update <- function (mu, sigma_sq, tau_sq, y, n){
    mean_y_j <- colSums(y)/n
    theta_estimate <- c(0,0,0,0)
    V <- c(0,0,0,0)
    for (j in 1: 4){
        theta_estimate[j] <- (mu/tau_sq + mean_y_j[j]*n[j]/sigma_sq) / (1/tau_sq + n[j]/sigma_sq)    
    }
    for (j in 1: 4){
        V[j] <- 1 / (1/tau_sq + n[j]/sigma_sq)    
    }
    value <- c()
    for (i in 1:4){
        value[i] <- rnorm(1, theta_estimate[i], sqrt(V[i]))
    }
    v_t <- value[1] < 0 || value[2] < 0 || value[3] < 0 || value[4] < 0
    while (v_t){
        for (i in 1:4){
            value[i] <- rnorm(1, theta_estimate[i], sqrt(V[i]))
        }
    }
    value
}

# conditional posterior for mu
mu.update <- function (theta, tau_sq){
    mu_estimate <- 1/4 * sum(theta)
    value <- rnorm(1, mu_estimate, sqrt(tau_sq/4))
    while (value < 0){
        value <- rnorm(1, mu_estimate, sqrt(tau_sq/4))
    }
    value
}

# conditional posterior for sigma^2
sigma_sq.update <- function(theta, y, n, N){
    residual <- matrix(0, nrow = 8, ncol = 4)
    # J = 4
    for (j in 1:4){
        for (i in 1:n[j]){
            residual[i,j] <- (y[i,j] - theta[j])^2
        }
    }
    SR <- sum(colSums(residual))
    sigma_sq_estimate <- 1/N * SR
    Y <- rchisq(1, N)
    value <- N * sigma_sq_estimate / Y
    value
}

# conditional posterior for tau^2
tau_sq.update <- function(theta, mu){
    # J-1 = 3
    tau_sq_estimate <- 1/3 * sum((theta - mu)^2)
    Y <- rchisq(1, 3)
    value <- 3 * tau_sq_estimate / Y
    value
}

# Gibbs Sampling function
gibbs <- function (num_sims, theta.start, mu.start, y, n, N){
    theta.draws <- matrix(0, nrow = num_sims, ncol = 4)
    mu.draws <- c()
    sigma.draws <- c()
    tau.draws <- c()
    
    theta.cur <- theta.start
    mu.cur <- mu.start
    
    for (k in 1:num_sims){
        tau_sq.cur <- tau_sq.update(theta.cur, mu.cur)
        sigma_sq.cur <- sigma_sq.update(theta.cur, y, n, N)
        mu.cur <- mu.update(theta.cur, tau_sq.cur)
        theta.cur <- theta.update(mu.cur, sigma_sq.cur, tau_sq.cur, y, n)
        theta.draws[k, ] <- theta.cur
        mu.draws[k] <- mu.cur
        sigma.draws[k] <- sqrt(sigma_sq.cur)
        tau.draws[k] <- sqrt(tau_sq.cur)
    }
    # remove first half simulations (as required for Gelman-Rubin diagnostic)
    theta.draws <- theta.draws[(floor(num_sims/2)+1) : num_sims, ]
    mu.draws <- mu.draws[(floor(num_sims/2)+1) : num_sims]
    sigma.draws <- sigma.draws[(floor(num_sims/2)+1) : num_sims]
    tau.draws <- tau.draws[(floor(num_sims/2)+1) : num_sims]
    
    results <- matrix(0, nrow <- floor(num_sims/2), 4+3)
    results[,1:4] <- theta.draws
    results[,5] <- mu.draws
    results[,6] <- sigma.draws
    results[,7] <- tau.draws
    results
}


# the following posterior is based on setting 'theta.start <- colSums(y)/n'
posterior <- gibbs(num_sims = 100, theta.start = theta.start, mu.start = mu.start, 
                   y = y, n = n, N = N)


#############################
## Gelman-Rubin diagnostic ##
#############################
#
# Note that DROP DRAWs are implemented in the Gibbs function,
# But it could also implemented using mcmc() function's parametr 'start' & 'end'
#

# set overdispersed starting points (as a function), by random sampling 
theta.start_new <- function(){
    theta_new <- c()
    for (i in 1:4){
        theta_new[i] <- sample(y[,i][which(y[,i] != 0)], 1)
    }
    theta_new
}

# if run 10 chains (as in the book)
posterior_1 <- gibbs(num_sims = 100, theta.start = theta.start_new(), 
                     mu.start = mu.start, y = y, n = n, N = N)
posterior_2 <- gibbs(num_sims = 100, theta.start = theta.start_new(), 
                     mu.start = mu.start, y = y, n = n, N = N)
posterior_3 <- gibbs(num_sims = 100, theta.start = theta.start_new(), 
                     mu.start = mu.start, y = y, n = n, N = N)
posterior_4 <- gibbs(num_sims = 100, theta.start = theta.start_new(), 
                     mu.start = mu.start, y = y, n = n, N = N)
# M = 4
chain1 <- mcmc(posterior_1)
chain2 <- mcmc(posterior_2)
chain3 <- mcmc(posterior_3)
chain4 <- mcmc(posterior_4)


mc_list <- mcmc.list(list(chain1, chain2, chain3, chain4))
summary(mc_list)
gelman.diag(mc_list)









