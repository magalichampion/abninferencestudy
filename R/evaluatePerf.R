#' Evaluate Node-Level ABN Inference Performance
#'
#' Computes prediction scores, bias, and (for discrete nodes) classification
#' accuracy for each node in an additive Bayesian network by comparing
#' observed values against sample-level or marginal predictions.
#'
#' @param observed A \code{data.frame} of observed values with \code{n} rows
#'   and columns corresponding to the node names in \code{dists}. 
#' @param predicted A named list indexed by node name containing either
#'   sample-level predictions (a list of length \code{n}, one entry per
#'   observation) or marginal predictions per node.
#' @param dists A named character vector or list mapping each node name to its
#'   distribution type (\code{"gaussian"}, \code{"binomial"}, \code{"poisson"},
#'   or \code{"multinomial"}).
#'
#' @return A \code{data.frame} (\code{PerfTable}) with row names set to node
#'   names and the following columns:
#' \describe{
#'   \item{node}{Character node name.}
#'   \item{node_type}{Distribution type for the node.}
#'   \item{bias}{Mean predicted value minus mean observed value.}
#'   \item{accuracy}{Classification accuracy for \code{"binomial"} and
#'     \code{"multinomial"} nodes, \code{NA} for \code{"gaussian"} and
#'     \code{"poisson"} nodes.}
#'   \item{score}{Node prediction score (MSE, Brier score, or relative
#'     deviance).}
#'   \item{scoreName}{Name of the score metric (\code{"MSE"},
#'     \code{"BrierScore"}, or \code{"Deviance"}).}
#' }
#'
#' @export
EvaluatePerf <- function(observed, predicted, dists){
  n <- nrow(observed)
  node_names <- names(dists)
  
  out_rows <- vector("list", length(node_names))
  
  for (j in seq_along(node_names)){
    node_name <- node_names[j]
    pred_node <- predicted[[node_name]]
    obs_node <- observed[[node_name]]
    dist_node <- dists[[node_name]]
    
    is_sample_level <- is.list(pred_node) && length(pred_node) == n 
    
    if (dist_node == "binomial") {
      obs_levels <- levels(obs_node)
      obs_numeric <- match(as.character(obs_node), levels(obs_node)) - 1
      
      if (is_sample_level) {
        pred_node <- lapply(pred_node, function(p_vec) {
          if (!is.null(names(p_vec))) p_vec[obs_levels] else p_vec
        })
      } else if (!is.null(names(pred_node))) {
        pred_node <- pred_node[obs_levels]
      }
    } else if (dist_node == "multinomial"){
      obs_levels  <- levels(obs_node)
      obs_numeric <- match(as.character(obs_node), obs_levels)
      
      if (is_sample_level) {
        pred_node <- lapply(pred_node, function(p_vec) {
          if (!is.null(names(p_vec))) {
            p_aligned <- p_vec[obs_levels]
            p_aligned[is.na(p_aligned)] <- 0
            return(p_aligned)
          } else {
            return(p_vec)
          }
        })
      } else {
        if (!is.null(names(pred_node))) {
          pred_node <- pred_node[obs_levels]
          pred_node[is.na(pred_node)] <- 0
        }
      }
    } else {
      obs_numeric <- as.numeric(obs_node)
    }
    
    avg_obs <- mean(obs_numeric, na.rm = TRUE) 
    
    if (dist_node == "gaussian"){
      mu <- if (is_sample_level) {
        vapply(pred_node, function(x) as.numeric(x)[1], numeric(1))
      } else {
        rep(as.numeric(pred_node)[1], n)
      }
      score_val <- mean((obs_numeric-mu)^2,na.rm=TRUE)
      bias <- mean(mu, na.rm = TRUE) - avg_obs
      score_name <- "MSE"
      accuracy <- NA
    } else if (dist_node == "binomial"){
      p <- if (is_sample_level) {
        vapply(pred_node, function(x) as.numeric(x)[2], numeric(1))
      } else {
        rep(as.numeric(pred_node)[2], n)
      }
      pred_class <- ifelse(p > 0.5, 1, 0)
      accuracy <- mean(pred_class == obs_numeric, na.rm = TRUE)
      
      score_val <- mean((obs_numeric - p)^2 + ((1 - obs_numeric) - (1 - p))^2, na.rm = TRUE)
      score_name <- "BrierScore"
      bias <- mean(p, na.rm = TRUE) - avg_obs
    } else if (dist_node == "poisson"){
      lambda <- if (is_sample_level) {
        vapply(pred_node, function(x) as.numeric(x)[1], numeric(1))
      } else {
        rep(as.numeric(pred_node)[1], n)
      }
      
      term1 <- obs_numeric * log(obs_numeric / (lambda + 1e-10))
      term1[obs_numeric == 0] <- 0
      dev_model <- 2 * sum(term1 - (obs_numeric - lambda))
      
      term2 <- obs_numeric * log(obs_numeric / (avg_obs + 1e-10))
      term2[obs_numeric == 0] <- 0
      dev_null <- 2 * sum(term2 - (obs_numeric - avg_obs))
      
      score_val <- dev_model/dev_null
      score_name <- "Deviance"
      bias <- mean(lambda, na.rm = TRUE) - avg_obs
      accuracy <- NA
    } else { # multi distribution
      if (is_sample_level) {
        p_actual <- sapply(1:n, function(i) pred_node[[i]][obs_numeric[i]])
        predictions_class <- sapply(pred_node, which.max)
        accuracy <- mean(predictions_class == obs_numeric, na.rm = TRUE)
        expected_category <- sapply(pred_node, function(p_vec) {
          sum(seq_along(p_vec) * p_vec)
        })
      } else {
        p_actual <- pred_node[obs_numeric]
        predictions_class <- which.max(pred_node)
        accuracy <- mean(predictions_class == obs_numeric, na.rm = TRUE)
        expected_category <- sum(seq_along(pred_node) * pred_node)
      }
      bias <- mean(expected_category, na.rm = TRUE) - avg_obs
      score_name <- "BrierScore"
      score_val  <- mean((1 - p_actual)^2)
    }
    
    out_rows[[j]] <- data.frame(
      node        = node_name,
      node_type   = dist_node,
      bias        = bias,
      accuracy  = accuracy,
      score    = score_val,
      scoreName    = score_name
    )
    
  }
  PerfTable <- do.call(rbind, out_rows)
  rownames(PerfTable) <- PerfTable$node
  return(PerfTable)
}
