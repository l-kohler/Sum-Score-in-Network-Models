library("parSim")
library("tidyverse")

# function to calculate expected sum score
exp_sum <- function(graph, thresholds, beta = 1, responses) {
  sumP <- IsingSumLikelihood(graph, thresholds, beta, responses)
  e_sum <- sum(sumP$Sum * sumP$P)
  return(e_sum)
}

# Simulation:
parSim(
  size = 10, # change to simulate other network sizes
  encoding = c(0,-1),
  model = c("SW_0", "SW_25", "SW_50","SW_75", "SW_100", "CW"), # remove Sw_25/_75 if size = 5
  thresholds = c(-3, -2, -1, 0, 1),
  meanWeight = c(0, 0.25, 0.5, 0.75, 1),
  sdWeight = 0.1,
  sdTholds = 0.4,
  reps = 1000,
  write = TRUE,
  # change name for other network sizes
  name = "~/BDS-Master/Master thesis/R/Sim_data/Simulation_10_new",
  nCores = 12,
  export = c("exp_sum"),
  progressbar = TRUE,
  expression = {
    # Load libraries
    library("igraph")
    library("IsingSampler")
    library("dplyr")
    library("MASS")
    
    # Generate network structure
    if (model == "SW_0"){
      network <- matrix(0, nrow = size, ncol = size)
      
    } else if (model == "SW_25"){ # remove if size = 5
      # calculate connectivity
      nei <- round((0.25 * (size - 1))/2, 0) %>% max(1)
      network <- sample_smallworld(dim = 1, # lattice dimension
                                 size = size, # network size
                                 nei = nei, # the number of closest neighbors a node is connected to
                                 p = 0.04, # rewiring probability 
                                 loops = FALSE, # autocorrelations not allowed
                                 multiple = FALSE) %>% # multiple edges not allowed
        # convert structure to edge-weight matrix
        as_adjacency_matrix() %>% # convert to sparse matrix
        as.matrix() #%>% # convert to dense matrix
      
    } else if (model == "SW_50"){
      nei <- round((0.5 * (size - 1))/2, 0) %>% max(1)
      network <- sample_smallworld(dim = 1, # lattice dimension
                                 size = size, # network size
                                 nei = nei, # the number of closest neighbors a node is connected to
                                 p = 0.04, # rewiring probability 
                                 loops = FALSE, # autocorrelations not allowed
                                 multiple = FALSE) %>% # multiple edges not allowed
        # convert structure to edge-weight matrix
        as_adjacency_matrix() %>% # convert to sparse matrix
        as.matrix() #%>% # convert to dense matrix
      
    } else if (model == "SW_75"){ # remove if size = 5
      nei <- round((0.75 * (size - 1))/2, 0) %>% max(1)
      network <- sample_smallworld(dim = 1, # lattice dimension
                                 size = size, # network size
                                 nei = nei, # the number of closest neighbors a node is connected to
                                 p = 0.04, # rewiring probability 
                                 loops = FALSE, # autocorrelations not allowed
                                 multiple = FALSE) %>% # multiple edges not allowed
        # convert structure to edge-weight matrix
        as_adjacency_matrix() %>% # convert to sparse matrix
        as.matrix() #%>% # convert to dense matrix
      
    } else if (model == "SW_100"){
      nei <- round((1 * (size - 1))/2, 0) %>% max(1)
      network <- sample_smallworld(dim = 1, # lattice dimension
                                   size = size, # network size
                                   nei = nei, # the number of closest neighbors a node is connected to
                                   p = 0.04, # rewiring probability 
                                   loops = FALSE, # auto correlations not allowed
                                   multiple = FALSE) %>% # multiple edges not allowed
        # convert structure to edge-weight matrix
        as_adjacency_matrix() %>% # convert to sparse matrix
        as.matrix() #%>% # convert to dense matrix
      
    } else if (model == "CW"){
      network <- matrix(1, size, size)
      diag(network) <- 0
      network <- network * rnorm(1, meanWeight, sdWeight)
    }
    
    # Thresholds per node
    tholds <- rnorm(size, thresholds, sdTholds)
    
    # Add random weights to small-world networks
    if (model != "CW") {
      network[lower.tri(network)] <- network[lower.tri(network)] * rnorm(sum(lower.tri(network)),
                                                                         meanWeight,
                                                                         sdWeight)
      network[upper.tri(network)] <- t(network)[upper.tri(network)]
      diag(network) <- 0
    }
    
    # write network away
    ## change file name for other network sizes
    write.matrix(network, file = paste0("~/BDS-Master/Master thesis/R/Networks/network_new/network_10_", id, ".txt"))
    
    # calculate expected and observed sum scores
    sum_e <- exp_sum(network, tholds, responses = c(encoding, 1)) # expected sum score
    config_o <- IsingSampler(n = 1, network, tholds, responses = c(encoding, 1)) # observed configuration
    sum_o <- sum(config_o[1,]) # observed sum score
    
    # store results:
    data.frame(
      sum_e = sum_e,
      sum_o = sum_o,
      config_o = list(config_o),
      tholds = list(matrix(tholds, nrow = 1))
      #network = list(network)
    )
  })
