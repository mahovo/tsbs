# ==============================================================================
# 1. SETUP & DATA GENERATION (Simulating a "Real" Market History)
# ==============================================================================
set.seed(42)

# Parameters for two regimes
# Regime 1: Calm (Low vol, low persistence)
omega1 <- 0.00001; alpha1 <- 0.05; beta1 <- 0.90
# Regime 2: Crisis (High vol, high persistence)
omega2 <- 0.00005; alpha2 <- 0.10; beta2 <- 0.85

# Transition Matrix (Stay in state with 98% / 95% prob)
P <- matrix(c(0.98, 0.02, 0.05, 0.95), nrow=2, byrow=TRUE)

n_obs <- 3000
states <- numeric(n_obs)
returns <- numeric(n_obs)
sigma2 <- numeric(n_obs)

# Initialize
states[1] <- 1
sigma2[1] <- omega1 / (1 - alpha1 - beta1)
returns[1] <- rnorm(1, 0, sqrt(sigma2[1]))

# Generate History
for(t in 2:n_obs) {
  # 1. Markov Switch
  prev_state <- states[t-1]
  probs <- P[prev_state, ]
  states[t] <- sample(1:2, 1, prob = probs)
  
  # 2. GARCH Dynamics based on current state
  if(states[t] == 1) {
    sigma2[t] <- omega1 + alpha1 * returns[t-1]^2 + beta1 * sigma2[t-1]
  } else {
    sigma2[t] <- omega2 + alpha2 * returns[t-1]^2 + beta2 * sigma2[t-1]
  }
  
  # 3. Generate Return
  returns[t] <- rnorm(1, 0, sqrt(sigma2[t]))
}

# ==============================================================================
# 2. IMPLEMENTATION OF OPTION B: STATE-BLOCK BOOTSTRAP (The "Bad" Way)
# ==============================================================================
run_block_bootstrap <- function(hist_returns, hist_states, n_sim) {
  
  # A. Identify Blocks (Runs)
  runs <- rle(hist_states)
  end_idx <- cumsum(runs$lengths)
  start_idx <- c(1, head(end_idx, -1) + 1)
  
  blocks_s1 <- list()
  blocks_s2 <- list()
  
  for(i in 1:length(runs$values)) {
    blk <- hist_returns[start_idx[i]:end_idx[i]]
    if(runs$values[i] == 1) blocks_s1[[length(blocks_s1)+1]] <- blk
    else                    blocks_s2[[length(blocks_s2)+1]] <- blk
  }
  
  # B. Generate Target Schedule (Using same Transition Matrix P)
  # We simulate a state path first
  sim_states <- numeric(n_sim)
  sim_states[1] <- 1
  for(t in 2:n_sim) sim_states[t] <- sample(1:2, 1, prob = P[sim_states[t-1],])
  
  # C. Stitch Blocks
  sim_runs <- rle(sim_states)
  boot_series <- c()
  
  for(k in 1:length(sim_runs$values)) {
    target_state <- sim_runs$values[k]
    needed_len <- sim_runs$lengths[k]
    
    pool <- if(target_state == 1) blocks_s1 else blocks_s2
    
    # Fill the required length by stitching random blocks
    segment <- c()
    while(length(segment) < needed_len) {
      rand_blk <- pool[[sample(length(pool), 1)]]
      segment <- c(segment, rand_blk)
    }
    boot_series <- c(boot_series, segment[1:needed_len])
  }
  return(boot_series)
}

# ==============================================================================
# 3. IMPLEMENTATION OF OPTION A: RESIDUAL BOOTSTRAP (The "Good" Way)
# ==============================================================================
run_residual_bootstrap <- function(hist_returns, hist_sigma2, hist_states, n_sim) {
  
  # A. Whiten the residuals (z = r / sigma)
  z <- hist_returns / sqrt(hist_sigma2)
  z_s1 <- z[hist_states == 1]
  z_s2 <- z[hist_states == 2]
  
  # B. Simulate
  sim_ret <- numeric(n_sim)
  sim_sig2 <- numeric(n_sim)
  sim_states <- numeric(n_sim)
  
  # Init
  sim_states[1] <- 1
  sim_sig2[1] <- var(hist_returns)
  sim_ret[1] <- rnorm(1, 0, sqrt(sim_sig2[1]))
  
  for(t in 2:n_sim) {
    # 1. State Transition
    sim_states[t] <- sample(1:2, 1, prob = P[sim_states[t-1],])
    curr_s <- sim_states[t]
    
    # 2. Resample Innovation
    z_star <- if(curr_s == 1) sample(z_s1, 1) else sample(z_s2, 1)
    
    # 3. GARCH Recursion (The Magic Step: Smooth Volatility)
    # Note: We switch parameters based on the new state
    if(curr_s == 1) {
      sim_sig2[t] <- omega1 + alpha1 * sim_ret[t-1]^2 + beta1 * sim_sig2[t-1]
    } else {
      sim_sig2[t] <- omega2 + alpha2 * sim_ret[t-1]^2 + beta2 * sim_sig2[t-1]
    }
    
    sim_ret[t] <- sqrt(sim_sig2[t]) * z_star
  }
  return(sim_ret)
}

# ==============================================================================
# 4. RUN THE COMPARISON
# ==============================================================================

# Run both bootstraps
series_B <- run_block_bootstrap(returns, states, n_sim = 3000)
series_A <- run_residual_bootstrap(returns, sigma2, states, n_sim = 3000)

# ==============================================================================
# 5. VISUALIZE THE ACF (The Proof)
# ==============================================================================
par(mfrow=c(3,1), mar=c(4,4,2,1))

# 1. Original Data ACF
acf(returns^2, lag.max = 50, main = "1. Original Data (Squared Returns)\nNote the smooth decay", col="black")

# 2. Option B (Block Stitching)
acf(series_B^2, lag.max = 50, main = "2. Option B: Block Bootstrap\nNote the faster drop-off or irregularity", col="red")

# 3. Option A (Residual Bootstrap)
acf(series_A^2, lag.max = 50, main = "3. Option A: Residual Bootstrap\nPreserves the 'memory' structure perfectly", col="blue")