#' Run ABN Inference Across Conditioning Evidence Types
#'
#' Computes predictions from a simulated additive Bayesian network (ABN) dataset
#' using \code{abn::predictAbn()}.
#' 
#' @param dataset A list containing the simulated ABN components: \code{data}
#'   (data frame of observations), \code{dists} (node distribution list),
#'   \code{DAG} (an \code{igraph} DAG object), \code{mycoefs} (coefficients),
#'   and \code{mse} (residual variances / MSEs).
#' @param evidence_type Character scalar specifying the conditioning set to use.
#'   Must be one of \code{"none"}, \code{"parents"}, \code{"children"},
#'   \code{"mb"}, or \code{"spouses"}.
#' @param rngSEED Optional integer seed passed to \code{set.seed()} for
#'   reproducibility. Defaults to \code{NULL}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{predictions}{When \code{evidence_type = "none"}, the per-node
#'     marginal predictions returned by \code{predictAbn(..., hypothesis = NULL)$predictions}.
#'     When conditioning on evidence, a named list of length equal to the number
#'     of nodes, where each element is a list of length \code{nrow(dataset$data)}
#'     containing per-sample hypothesis predictions.}
#'   \item{elapsed_sec}{Numeric scalar giving elapsed time in seconds.}
#' }
#'
#' @export
#' @importFrom future plan multisession
#' @importFrom furrr future_map furrr_options
#' @importFrom igraph as_adjacency_matrix neighbors
RunABN_inference <- function(dataset, evidence_type, rngSEED){
  t_start <- proc.time()
  
  if (!is.null(rngSEED)) {
    set.seed(rngSEED)
  }

  valid_types <- c("none", "parents", "children", "mb","spouses")
  if (!evidence_type %in% valid_types) {
    stop(paste("evidence_type must be one of:", paste(valid_types, collapse = ", ")))
  }
  
  df <- dataset$data
  dists <- dataset$dists
  dag <- t(igraph::as_adjacency_matrix(dataset$DAG,sparse=FALSE))
  coefs <- dataset$mycoefs
  vars <- dataset$mse

  AbnFit_object <- list(coef=coefs,mse=vars,abnDag=list(dag=dag,data.dists=dists))
  
  if (evidence_type == "none"){
    predictions <- predictAbn(object = AbnFit_object, data = df, hypothesis = NULL)
    predictions <- predictions$predictions
  } else {
    old_plan <- future::plan(future::multisession, workers = 2)
    on.exit(future::plan(old_plan), add = TRUE)
    
    predictions <- list()
    for (i in 1:ncol(dag)) {
      node_name <- colnames(dag)[i]
      
      if (evidence_type == "mb") {
        mb <- find_MB(graph = dataset$DAG, hypothesis = node_name)
        nodes_pool <- setdiff(mb, node_name)
      } else if (evidence_type == "parents"){
        nodes_pool <- igraph::neighbors(dataset$DAG, v = node_name, mode = "in")$name
      } else if (evidence_type == "children"){
        nodes_pool <- igraph::neighbors(dataset$DAG, v = node_name, mode = "out")$name
      } else if (evidence_type == "spouses"){
        children <- igraph::neighbors(dataset$DAG, v = node_name, mode = "out")$name
        spouses <- character(0)
        for (child in children) {
          child_parents <- names(igraph::neighbors(dataset$DAG, v=child, mode = "in"))
          spouses <- c(spouses, child_parents)
        }
        nodes_pool <- c(children, spouses)
        nodes_pool <- setdiff(nodes_pool, node_name)
      }
      
      predictions_per_sample <- furrr::future_map(1:nrow(df), function(j) {
        evidence <- as.list(df[j, nodes_pool, drop = FALSE])
        
        pred_out <- predictAbn(object = AbnFit_object,
                               data = df,
                               hypothesis = node_name,
                               evidence = evidence
        )
        
        return(pred_out$prediction_hypothesis)
      }, .options = furrr_options(seed = TRUE)) 
      predictions <- c(predictions,list(predictions_per_sample))
    }
    names(predictions) <- colnames(dag)
  }
  
  elapsed <- proc.time() - t_start
  
  return(list(predictions = predictions,timing = unname(elapsed["elapsed"])))
}
