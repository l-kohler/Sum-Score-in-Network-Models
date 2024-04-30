# function to calculate expected node state
exp_state <- function(graph, thresholds, size, beta = 1, responses) {
  configP <- IsingLikelihood(graph, thresholds, beta, responses)
  e_state <- numeric(length = size)
  for (i in 1:size) {
    e_state[i] <- configP %>%
      dplyr::select(paste0("Var", i), Probability) %>%
      mutate(exp = .[[1]] * .[[2]]) %>%
      dplyr::select(exp) %>%
      sum()
  }
  return(e_state)
}

# compute error correlations
## register back end
cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)
clusterExport(cl = cl, c("exp_state", "data"))
clusterEvalQ(cl = cl, library("tidyverse"))
clusterEvalQ(cl = cl, library("IsingSampler"))

## simulation
error_corr_M <- 
  foreach (i = 1:nrow(data), .combine = "c") %dopar% {
    # read
    network <- read.table(paste0("~/BDS-Master/Master thesis/R/Networks/network_new/network_10_", i, ".txt")) %>%
    as.matrix() %>% matrix(nrow = 10, ncol = 10)
    config_o <- IsingSampler(n = 1000, network,
                             data$tholds[i,],
                             responses = c(data$encoding[i], 1))
  for (j in 1:10) {
    if (near(var(config_o[,j]), 0)) {
      config_o[,j] <- jitter(config_o[, j])
    }
  }
  config_e <- exp_state(network,
                        data$tholds[i,],
                        data$size[i],
                        responses = c(data$encoding[i], 1))
  error <- sweep(config_o, 2, config_e)
  error_corr <- cor(error)
  error_corr <- error_corr[upper.tri(error_corr)] %>% matrix(nrow = 1)
  error_corr_M[i] <- error_corr %>% mean()
}
# remove backend
parallel::stopCluster(cl)

## output
error_corr_M

#save object
#save(error_corr_M, file="~/BDS-Master/Master thesis/R/error_cor.Rdata")

# add mean error to data
data$error_corr_M <- NA
data$error_corr_M <- error_corr_M

# save new data frame
write.table(data, file = "~/BDS-Master/Master thesis/R/Sim_data/Simulation_10_error_cor.txt",
            sep = " ", row.names = FALSE)
