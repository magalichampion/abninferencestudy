#' Simulate Additive Bayesian Network (ABN) Data
#'
#' Generates synthetic datasets based on Directed Acyclic Graphs (DAGs) with specific
#' characteristics and node distributions using the \code{abn} package.
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Number of nodes in the graph.
#' @param nodedists Character. Distribution type for nodes. Must be one of 
#'   \code{"balanced"}, \code{"gaussian"}, \code{"poisson"}, \code{"binomial"}, or \code{"multinomial"}.
#' @param d Numeric. Density level (proportion of edges to keep in the graph, strictly between 0 and 1).
#' @param max_parents Numeric. Maximal number of parents allowed per node.
#' @param magnitude Numeric vector of length 2. Minimum and maximum range for uniform coefficients. 
#'   Default is \code{c(-1, 1)}.
#' @param MB List. Controls Markov blanket size restrictions. Contains \code{MB_size} (desired size) 
#'   and \code{node} (target node name). Defaults to \code{list(MB_size = NULL, node = NULL)}.
#' @param n_graphs Integer. Number of datasets to simulate and save.
#' @param path.outfile Character string. Root directory path where \code{data/} and \code{data-raw/} 
#'   folders will be created and saved.
#' @param rngSEED Integer. Random number generator seed for reproducibility.
#' @param counter Integer. Optional starting index for file naming counters. Starts by default at 1.
#'
#' @return Invisible. Validated simulation objects are saved directly to as `.Rdata` files 
#'   inside the specified \code{path.outfile/data/} directory.
#'
#' @importFrom igraph graph_from_adjacency_matrix as_adjacency_matrix is_connected
#' @import rjags
#' @export
simulateABNdata <- function(n, p, nodedists, d, max_parents, magnitude = NULL, MB = list(MB_size = NULL, node = NULL), n_graphs, path.outfile, rngSEED, counter = NULL) {
  if (!is.null(rngSEED)) {
    set.seed(rngSEED)
  }
  if (!is.numeric(d) || length(d) != 1 || d <= 0 || d >= 1) {
    stop("'d' must be a numeric value strictly between 0 and 1.")
  }
  
  if (!is.numeric(n) || n <= 0 || length(n) != 1 || floor(n) != n) {
    stop("'n' must be a positive integer.")
  }
  
  if (!is.numeric(p) || p < 0 || length(p) != 1 || floor(p) != p) {
    stop("'p' must be a positive integer.")
  }
  
  if (!is.numeric(max_parents) || length(max_parents) != 1 || max_parents < 0) {
    stop("'max_parents' must be a numeric value greater than 0.")
  }
  
  valid_types <- c("balanced","G","P","B","M")
  if (!nodedists %in% valid_types) {
    stop(paste("nodedists must be one of:", paste(valid_types, collapse = ", ")))
  }
  
  dists <- rep(1,4)
  names(dists) <- c("gaussian","poisson","binomial","multinomial")
  if (nodedists == "G"){
    dists["gaussian"] <- p-3
  } else if (nodedists == "P"){
    dists["poisson"] <- p-3
  } else if (nodedists == "B"){
    dists["binomial"] <- p-3
  } else if (nodedists == "M"){
    dists["multinomial"] <- p-3
  } else if (nodedists == "balanced"){
    dists <- rep(p/4,4)
    names(dists) <- c("gaussian","poisson","binomial","multinomial")
  }

  if (!is.numeric(n_graphs) || n_graphs <= 0 || length(n_graphs) != 1 || floor(n_graphs) != n_graphs) {
    stop("'n_graphs' must be a positive integer.")
  }
  
  if (!is.character(path.outfile) || length(path.outfile) != 1) {
    stop("'path.outfile' must be a single character string (path prefix).")
  }
  
  # Ensure target directories exist to prevent path errors
  dir.create(file.path(path.outfile, "data"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.outfile, "data-raw"), showWarnings = FALSE, recursive = TRUE)
  
  if (is.null(magnitude)){
    magnitude <- c(-1,1) 
  } 
  
  # 1st step: generate the nodes distribution
  p_distrib <- list(distrib = dists, levels_multi = sample(3:5,dists["multinomial"],replace = TRUE))
  dists <- as.list(rep(names(p_distrib$distrib),p_distrib$distrib))
  names(dists) <- paste0(substr(rep(names(p_distrib$distrib), p_distrib$distrib), 1, 1), # First letter of each distribution name
                         sequence(p_distrib$distrib)                     # Sequence of numbers for each distribution
  )
  
  # 2nd step: check if the combination is valid
  Check <- check_combinations(p = p, d, max_parents, nodedists, p_distrib)
  
  if (is.null(Check)){
    nbr_graphs <- 2000
    
    adj_matrix_total <- generate_RandomDAGs(p_distrib, d = d, max_parents = max_parents, nbr_graphs = nbr_graphs, MB = MB)
    
    dags <- lapply(adj_matrix_total, function(A){
      graph_from_adjacency_matrix(A)
    })
    
    if (is.null(counter)) {
      counter <- 1L
    }
    successful <- 0L
    
    for (g in seq_along(dags)){
      if (successful >= n_graphs) break
      
      DAG = dags[[g]]
      dag = t(igraph::as_adjacency_matrix(DAG,sparse=FALSE))
              
      if (!is.null(MB$MB_size)){
        outfile <- file.path(path.outfile, "/data-raw/", "n_",n,"_p_",p,"_nodedists_",nodedists,"_d_",gsub('[.]','',d),"_mp_",max_parents,"_MB_",MB$MB_size,"_node_",MB$node,"_graph",counter,".txt")
      } else {
        outfile <- file.path(path.outfile, "/data-raw/", "n_",n,"_p_",p,"_nodedists_",nodedists,"_d_",gsub('[.]','',d),"_mp_",max_parents,"_graph",counter,".txt")
      }
      
      # generation of the coefficients
      intercepts <- as.list(runif(p, min = -1, max = 1))
      names(intercepts) <- names(dists)
      if (length(grep("m",names(intercepts)))>0){
        for (i in (1:length(grep("m",names(intercepts))))){
          multi_node <- names(intercepts)[grep("m",names(intercepts))][i]
          intercepts[[multi_node]] <- runif(p_distrib$levels_multi[i]-1, min = -1, max = 1)
          names(intercepts[[multi_node]]) <- as.character(c(2:p_distrib$levels_multi[i]))
        }
      }
      if (length(grep("p",names(intercepts)))>0){
        poisson_nodes <- names(intercepts)[grep("p",names(intercepts))]
        intercepts[poisson_nodes] <- as.list(runif(length(poisson_nodes), min = -1, max = 0))
      }
      
      all_node_levels <- rep(1, length(dists))
      names(all_node_levels) <- names(dists)
      multi_indices <- grep("m", names(dists))
      if(length(multi_indices) > 0) {
        all_node_levels[multi_indices] <- p_distrib$levels_multi
      }
      
      parents_coefs <- list()
      
      for (i in seq_along(dists)) {
        nodename <- names(dists)[i]
        parent_names <- find_parents(graph = DAG, node = nodename)
        
        if (length(parent_names) == 0) {
          parents_coefs[[nodename]] <- numeric(0)
        } else {
          if (dists[[i]] == "multinomial") {
            L_child <- all_node_levels[nodename]
            
            col_info <- list()
            for (pa in parent_names) {
              L_p <- all_node_levels[pa]
              if (L_p > 1) {
                dummy_names <- paste0(pa, "", 2:L_p)
                for (dn in dummy_names) col_info[[dn]] <- dists[pa] 
              } else {
                col_info[[pa]] <- dists[pa]
              }
            }
            col_names <- names(col_info)
            
            mat_cols <- lapply(col_names, function(col_n) {
              if (length(grep("p", col_n)) > 0) {
                runif(L_child - 1, min = magnitude[1]/10, max = magnitude[2]/10)
              } else {
                runif(L_child - 1, min = magnitude[1], max = magnitude[2])
              }
            })
            mat <- do.call(cbind, mat_cols)
            
            colnames(mat) <- col_names
            parents_coefs[[nodename]] <- mat
          } else {
            coef_vector <- c()
            for (pa in parent_names) {
              L_p <- all_node_levels[pa]
              if (L_p > 1) {
                # Multinomial parent: create L-1 dummy coefficients
                dummies <- runif(L_p - 1, min = magnitude[1], max = magnitude[2])
                names(dummies) <- paste0(pa, "", 2:L_p)
                coef_vector <- c(coef_vector, dummies)
              } else if (dists[pa]=="poisson"){
                val <- runif(1, min = magnitude[1]/10, max = magnitude[2]/10)
                names(val) <- pa
                coef_vector <- c(coef_vector, val)
              } else {
                # Standard parent
                val <- runif(1, min = magnitude[1], max = magnitude[2])
                names(val) <- pa
                coef_vector <- c(coef_vector, val)
              }
            }
            parents_coefs[[nodename]] <- coef_vector
          }
        }
      }
      
      mymse <- rep(0,p)
      names(mymse) <- names(dists)
      mymse[which(dists=="gaussian")] <- 1
      
      mycoefs <- list()
      all_nodes <- names(intercepts)
      
      for (node in all_nodes) {
        current_mu <- intercepts[[node]]
        if (length(current_mu) > 1) {
          # Multinomial case: add suffixes like .2, .3, .4
          names(current_mu) <- paste0(node, "|intercept.", seq_along(current_mu) + 1)
        } else {
          # Standard case: Gaussian, Binomial, Poisson
          names(current_mu) <- paste0(node, "|intercept")
        }
        
        current_betas <- parents_coefs[[node]]
        
        if (is.matrix(current_betas)) {
          beta_names <- expand.grid(colnames(current_betas), seq_len(nrow(current_betas)) + 1)
          flat_beta_names <- paste0(beta_names$Var1, ".", beta_names$Var2)
          
          current_betas_vec <- as.vector(t(current_betas))
          names(current_betas_vec) <- flat_beta_names
        } else {
          # Already a vector or numeric(0)
          current_betas_vec <- current_betas
        }
        
        combined_vec <- c(current_mu, current_betas_vec)
        
        mycoefs[[node]] <- matrix(combined_vec,
                                  nrow = 1,
                                  byrow = TRUE,
                                  dimnames = list(NULL, names(combined_vec)))
      }
      
      AbnObject <- list("abnDag"=list("dag"=dag ,"data.dists"=dists), "mse" = mymse, "coef" = mycoefs)
      try.data <- try(simulateAbn(
        object = AbnObject, 
        n = n,
        seed = rngSEED
      ), TRUE)
      if (inherits(try.data, "try-error")) next
      
      data <- try.data
      simData <- list(data=data,DAG=DAG,adj.matrix = t(dag), dists = dists, mycoefs=mycoefs,mse=mymse,max.parents=max_parents)
      
      Check <- check_data(simData)
      
      if (is.null(Check)){
        
        # rename the object data
        if (!is.null(MB$MB_size)){
          newdataname <- paste0("n_",n,"_p_",p,"_nodedists_",nodedists,"_d_",gsub('[.]','',d),"_mp_",max_parents,"_MB_",MB$MB_size,"_node_",MB$node,"_graph",counter)
        } else {
          newdataname <- paste0("n_",n,"_p_",p,"_nodedists_",nodedists,"_d_",gsub('[.]','',d),"_mp_",max_parents,"_graph",counter)
        }
        assign(newdataname, simData)
        full_path <- file.path(path.outfile, paste0("data/",newdataname, ".RData"))
        save(list = newdataname, file = full_path)
        
        rm(newdataname)
        
        counter <- counter + 1L
        successful <- successful +1L
      }
    }
    
    if (successful < n_graphs) {
      warning(
        "Only ", successful, " valid datasets were generated, although ",
        n_graphs, " were requested. Available candidate DAGs: ",
        length(dags), "."
      )
    }
  }
}

#' Generate random Directed Acyclic Graphs (DAGs) with specified density, number of nodes, and max number of parents
#'
#' This function generates a list of random DAG adjacency matrices based on the
#' node distribution and desired number of nodes, density level, and max number of parents
#'
#' @param p_distrib List specifying the number of nodes per distribution
#'    and multinomial levels.
#' @param max_parents. Numeric. Maximal number of parents in the graph.
#' @param d Numeric. Density level (proportion of edges to keep in the graph, between 0 and 1).
#' @param nbr_graphs Integer. Number of random DAGs to generate.
#' @param MB List. Controls Markov blanket size restrictions. Contains \code{MB_size} (desired size) 
#'   and \code{node} (target node name). Defaults to \code{list(MB_size = NULL, node = NULL)}.
#' @return A list of adjacency matrices, each representing a DAG. Each matrix has
#'   dimension equal to the total number of nodes.
#'
#' @importFrom igraph graph_from_adjacency_matrix is_connected
#' @export
generate_RandomDAGs <- function(p_distrib, max_parents, d, nbr_graphs, MB = list(MB_size = NULL, node = NULL)){
  if (!is.list(p_distrib) || length(p_distrib) != 2 || length(p_distrib$levels_multi) != p_distrib$distrib["multinomial"] || length(p_distrib$distrib) <= 0 || !all(names(p_distrib$distrib) %in% c("gaussian","binomial","multinomial","poisson"))) {
    stop("'p_distrib' should be a list containing two elements: number of nodes per distribution (only Gaussian, binomial, multinomial and Poisson are supported) and multinomial levels for multinomial nodes.")
  }
  
  if (!is.numeric(d) || length(d) != 1 || d <= 0 || d >= 1) {
    stop("'d' must be a numeric value strictly between 0 and 1.")
  }
  
  if (!is.numeric(max_parents) || length(max_parents) != 1 || max_parents < 0 || max_parents > (sum(p_distrib$distrib)-1)) {
    stop(paste0("'max_parents' must be a numeric value between 0 and ",(sum(p_distrib$distrib)-1),"."))
  }
  
  if (!is.numeric(nbr_graphs) || nbr_graphs <= 0 || floor(nbr_graphs) != nbr_graphs) {
    stop("'nbr_graphs' must be a positive integer.")
  }
  
  p <- sum(p_distrib$distrib)
  node_names <- paste0(substr(rep(names(p_distrib$distrib), p_distrib$distrib), 1, 1), # First letter of each distribution name
                       sequence(p_distrib$distrib))     # Sequence of numbers for each distribution
  node_names_short <- substr(rep(names(p_distrib$distrib), p_distrib$distrib), 1, 1)
  
  adj_matrices <- list()
  count <- 1
  
  T <- matrix(0, ncol=p, nrow=p)
  T[lower.tri(T, diag = TRUE)] <- 1
  T <- T-diag(diag(T))
  
  target_edges <- round(p * (p - 1) / 2 * d)
  
  while (length(adj_matrices)<nbr_graphs && count < 1e5){
    T_tmp <- T
    
    for (j in 1:p) {
      potential_parents <- which(T_tmp[, j] == 1)
      
      if (length(potential_parents) > max_parents) {
        drop_rows <- sample(potential_parents, length(potential_parents) - max_parents)
        T_tmp[drop_rows, j] <- 0
      }
    }
    current_edges <- sum(T_tmp)
    
    if (current_edges > target_edges) {
      # --- CASE 1: Too many edges, trim excess ---
      active_edges <- which(T_tmp == 1, arr.ind = TRUE)
      num_to_drop <- current_edges - target_edges
      drop_indices <- sample(nrow(active_edges), num_to_drop)
      
      for (k in drop_indices) {
        T_tmp[active_edges[k, 1], active_edges[k, 2]] <- 0
      }
      
    } else if (current_edges < target_edges) {
      # --- CASE 2: Too few edges, add missing ones safely ---
      while (sum(T_tmp) < target_edges) {
        # Find all inactive slots in the lower triangle
        inactive_edges <- which(T_tmp == 0 & lower.tri(T_tmp), arr.ind = TRUE)
        if (nrow(inactive_edges) == 0) break # Safety check if grid is full
        
        # Filter to only keep slots whose target column hasn't hit max_parents yet
        col_sums <- colSums(T_tmp)
        valid_slots <- inactive_edges[col_sums[inactive_edges[, 2]] < max_parents, , drop = FALSE]
        
        if (nrow(valid_slots) == 0) break # No more slots can accept parents without breaking max_parents
        
        # Randomly pick one valid slot and turn it on
        add_idx <- sample(nrow(valid_slots), 1)
        T_tmp[valid_slots[add_idx, 1], valid_slots[add_idx, 2]] <- 1
      }
    }
    
    if (sum(T_tmp) != target_edges) {
      count <- count + 1
      next
    }
    
    #edges_to_remove <- sample(which(T==1),round(p*(p-1)/2 - p*(p-1)/2*d))
    #T_tmp[edges_to_remove] <- 0
    perm <- sample(1:p)
    perm <- node_names[perm]
    perm_matrix <- matrix(0, ncol=p, nrow=p)
    colnames(perm_matrix) <- node_names
    row_indices <- seq_along(perm)
    col_indices <- match(colnames(perm_matrix), perm)
    perm_matrix[cbind(row_indices, col_indices)] <- 1
    adj_matrix <- perm_matrix %*% T_tmp %*% t(perm_matrix)
    rownames(adj_matrix) = colnames(adj_matrix) <- node_names
    graph <- graph_from_adjacency_matrix(adj_matrix)
    
    is_valid <- (length(which(rowSums(adj_matrix) + colSums(adj_matrix) == 0)) == 0 && 
                   igraph::is_connected(graph, mode = "weak"))
    
    if (is_valid && !is.null(max_parents)) {
      if (max(colSums(adj_matrix)) > max_parents) {
        is_valid <- FALSE
      }
    }
    
    if (is_valid && !is.null(MB$MB_size)) {
      mb_nodes <- find_MB(graph, MB$node)
      if (length(mb_nodes) != MB$MB_size) {
        is_valid <- FALSE
      }
    }
    
    if (is_valid) {
      adj_matrices <- c(adj_matrices, list(adj_matrix))
    }
    count <- count +1
  }
  return(adj_matrices)
}

#' Check validity of simulation parameter combinations
#'
#' Performs several checks to validate if the combination of input parameters for
#' the simulation is valid according to predefined rules.
#'
#' @param p Integer. Number of nodes in the graph.
#' @param d Numeric. Density level (proportion of edges to keep in the graph, strictly between 0 and 1).
#' @param max_parents Numeric. Maximal number of parents allowed per node.
#' @param nodedists Character. Distribution type for nodes. Must be one of 
#'   \code{"balanced"}, \code{"gaussian"}, \code{"poisson"}, \code{"binomial"}, or \code{"multinomial"}.
#' @param p_distrib List. Output from \code{generate_distributions}, containing node counts per distribution and multinomial levels.
#'
#' @return Character vector of error messages if any checks fail, empty vector otherwise.
#'
#' @details
#' Checks include:
#' \itemize{
#'   \item Connectivity constraint: Ensures the density \code{d} yields enough edges to guarantee a connected graph ($\ge p - 1$ edges).
#'   \item Parent capacity constraint: Verifies that the required number of edges does not exceed the maximum capacity imposed by \code{max_parents}.
#'   \item Lower-triangular geometry limit: Checks that the edge count does not exceed the absolute maximum capacity of the lower-triangular adjacency matrix layout.
#'   \item Valid parameter range: Ensures that \code{max_parents} does not exceed \code{p - 1}.
#' }
#' @export
check_combinations <- function(p, d, max_parents, nodedists, p_distrib){
  # Input checks
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p != round(p)) {
    stop("`p` must be a single positive integer.")
  }
  if (!is.numeric(d) || length(d) != 1 || d < 0 || d > 1) {
    stop("`d` must be a numeric value between 0 and 1.")
  }
  
  Check <- c()
  
  # all nodes should be connected in the graph
  nbr_edges <- round(p*(p-1)/2*d,0)
  
  # 1. Check connectivity limit (too sparse)
  if (nbr_edges < (p - 1)) {
    Check <- c(Check,paste0("Density 'd = ", d, "' is too small for ", p, 
                            " nodes. It leaves only ", nbr_edges, 
                            " edges, but at least ", p - 1, " are required for connectivity."))
  }
  
  # 2. Check max_parents capacity limit (too dense for the parent cap)
  max_possible_edges_under_cap <- p * max_parents
  if (nbr_edges > max_possible_edges_under_cap) {
    Check <- c(Check,paste0("Combination impossible: density 'd = ", d, 
                            "' requires ", nbr_edges, " edges, but 'max_parents = ", 
                            max_parents, "' caps the total network at ", max_possible_edges_under_cap, " edges."))
  }
  
  # 3. Second check for max_parents capacity limit
  col_capacities <- p - (1:(p - 1)) 
  actual_max_edges <- sum(pmin(col_capacities, max_parents))
  
  if (nbr_edges > actual_max_edges) {
    Check <- c(Check, paste0("Combination impossible: density 'd = ", d, 
                             "' requires ", nbr_edges, " edges, but the lower-triangular geometry with 'max_parents = ", 
                             max_parents, "' caps the absolute maximum network capacity at ", actual_max_edges, " edges."))
  }
  
  # 4. Check that max_parents is not too large
  if (max_parents > (p-1)) {
    Check <- c(Check,paste0("Combination impossible: max_parents should be smaller than ",p-1))
  }
  return(Check)
}

#' Validate simulated ABN dataset structure and properties
#'
#' Performs checks on a generated simulation list object to ensure node levels, 
#' continuous ranges, missing values, and graph-to-data consistency meet structural requirements.
#'
#' @param data List. A simulation object containing at least elements \code{data} (a data frame 
#'   of observations) and \code{DAG} (the graph structure).
#'
#' @return A character vector of error messages if any checks fail, returns \code{NULL} 
#'   if all validations pass successfully.
#'
#' @export
check_data <-function(data){
  # input check
  if (!is.data.frame(data$data)) {
    stop("Input `data` must be a data.frame.")
  }
  
  bin_nodes <- grep("b",names(data$data))
  multi_nodes <- grep("m",names(data$data))
  gaussian_nodes <- grep("g",names(data$data))
  poisson_nodes <- grep("p",names(data$data))
  
  Check <- c()
  
  bin_nodes_levels <- sapply(bin_nodes, function(b){
    length(levels(data$data[[b]]))
  })
  if (length(which(bin_nodes_levels==1))>0){
    Check <- c(Check, paste("Checks fail since at least one bin node has only one level."))
  }
  
  bin_nodes_obs <- sapply(bin_nodes, function(b){
    t <- table(data$data[[b]])
    length(which(t==0))
  })
  if (any(bin_nodes_obs>0)){
    Check <- c(Check, paste("Checks fail since at least one bin node has one level with no observation."))
  }
  
  multi_nodes_levels <- sapply(multi_nodes, function(m){
    length(levels(data$data[[m]]))
  })
  if (length(which(multi_nodes_levels<3))>0){
    Check <- c(Check, paste("Checks fail since at least one multi node has only less than 3 levels."))
  }
  
  multi_nodes_obs <- sapply(multi_nodes, function(m){
    t <- table(data$data[[m]])
    length(which(t==0))
  })
  if (any(multi_nodes_obs>0)){
    Check <- c(Check, paste("Checks fail since at least one multi node has one level with no observation."))
  }
  
  cont_nodes <- c(gaussian_nodes,poisson_nodes)
  cont_obs <- sapply(cont_nodes, function(c){
    col_data <- data$data[[c]]
    max_abs <- max(abs(col_data), na.rm = TRUE)
    return(max_abs)
  })
  if (any(cont_obs>1e3)){
    Check <- c(Check, paste("Checks fail since at least one continuous node has exteme values."))
  }
  if (any(multi_nodes_obs>0)){
    Check <- c(Check, paste("Checks fail since at least one multi node has one level with no observation."))
  }
  
  if (any(is.na(data$data))){
    Check <- c(Check, paste("Checks fail since there is at least one missing value."))
  }
  
  if (ncol(data$data) != (length(data$DAG))){
    Check <- c(Check, paste("Checks fail since the number of nodes in the graph does not correspond to the number of variables in the data."))
  }
  
  comb <- sapply(colnames(data$data), function(i){
    parents <- find_parents(data$DAG, node = i)
    if (length(parents) == 0) {
      return(1) # Base case: No parents means 1 baseline state
    }
    multi_parents <- parents[grep("m", parents)]
    if (length(multi_parents)>0){
      levels <- sapply(multi_parents, function(m) {
        p_val <- data$data[[m]]
        length(levels(p_val))
      })
      prod(levels) * length(grep("m",parents))
    } else {
      length(grep("b",parents))
    }
  })
  
  if (max(comb) > 5000) {
    Check <- c(Check, paste("Checks fail since one node has too many parents."))
  }
  
  return(Check)
}