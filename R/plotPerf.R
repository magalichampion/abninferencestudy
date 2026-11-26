#' Merge ABN Inference Performance Score Files
#'
#' Searches a results directory for `Scores_*.RData` files, optionally filters
#' by simulation design parameters (`p`, `d`, `mp`, `node_dists`), loads each
#' matching `Scores_<name_data>` object into an isolated environment, attaches
#' parsed metadata columns, and combines them into a single data frame.
#'
#' @param results_dir Character scalar giving the path to the directory
#'   containing `Scores_*.RData` files (typically `here::here("results")`).
#' @param p Optional numeric or integer vector of graph sizes (`p`) to keep.
#'   When `NULL` (default), all graph sizes are included.
#' @param d Optional numeric vector of target graph densities (`d`) to keep.
#'  When `NULL` (default), all densities are included.
#' @param mp Optional numeric or integer vector of maximum parent counts (`mp`)
#'   to keep. When `NULL` (default), all `mp` values are included.
#' @param node_dists Optional character vector of node distribution strings
#'   (e.g., `"balanced"`, `"G"`, `"P"`, `"B"`, `"M"`) to keep. When `NULL`
#'   (default), all node distribution settings are included.
#'
#' @return A `data.frame` combining all matching score tables across graphs and
#'   filter settings, with added columns `graph`, `p`, `d`, `mp`, and
#'   `node_dists`. 
#'
#' @importFrom dplyr bind_rows
#' @importFrom tools file_path_sans_ext
#' @export
merge_scores <- function(results_dir, p = NULL, d = NULL, mp = NULL, node_dists = NULL) {
  files <- list.files(
    path = results_dir,
    pattern = "Scores_.*\\.Rdata",
    full.names = TRUE
  )
  files <- files[!grepl("_MB_", basename(files), fixed = TRUE)] # remove MB files
  
  if (length(files) == 0) {
    warning("No score files found in: ", results_dir)
    return(data.frame())
  }
  
  if (!is.null(d)){
    d_formatted <- sprintf("%02d", as.integer(d * 10))
  }
  
  rows <- list()
  idx <- 1L
  
  for (f in files) {
    fname <- tools::file_path_sans_ext(basename(f))
    name_data  <- sub("^[Ss]cores_", "", fname)
    
    # Extract design factors from filename when present
    p_match <- regmatches(name_data, regexpr("(?<=^|_)p_[0-9]+(?=_|$)", name_data, perl = TRUE))
    p_val   <- as.numeric(sub("^p_", "", p_match)) 
    
    d_match <- regmatches(name_data, regexpr("(?<=^|_)d_[0-9]+(?=_|$)", name_data, perl = TRUE))
    d_code  <- sub("^d_", "", d_match) 
    d_val   <- as.numeric(d_code) / 10 
    
    mp_match <- regmatches(name_data, regexpr("(?<=^|_)mp_[0-9]+(?=_|$)", name_data, perl = TRUE))
    mp_val   <- as.numeric(sub("^mp_", "", mp_match)) 
    
    nd_match <- regmatches(name_data, regexpr("(?<=^|_)nodedists_.+?(?=_(?:d_|mp_|graph)[0-9_]|$)", name_data, perl = TRUE))
    nd_val   <- sub("^nodedists_", "", nd_match) 
    
    g_match  <- regmatches(name_data, regexpr("graph[0-9]+$", name_data))
    graph_id <- sub("^graph", "", g_match) 
    
    # Filter only when arguments are explicitly provided
    if (!is.null(p) && !(p_val == p)) next
    if (!is.null(d) && !(d_code == d_formatted)) next
    if (!is.null(mp) && !(mp_val == mp)) next
    if (!is.null(node_dists) && !(nd_val == node_dists)) next
    
    env <- new.env(parent = emptyenv())
    load(f, envir = env)
    scores_name <- grep(paste0("^Scores_", name_data, "$"), ls(envir = env), value = TRUE)
    if (length(scores_name) == 0) {
      stop("Object Scores_", name_data, " not found in file: ", f)
    }
    scores <- env[[scores_name[1]]]
    
    scores$graph <- graph_id
    scores$p <- p_val
    scores$d <- d_val
    scores$mp <- mp_val
    scores$node_dists <- nd_val
    
    rows[[idx]] <- scores
    idx <- idx + 1L
  }
  
  if (length(rows) == 0) {
    warning("No score files matched the specified filters in: ", results_dir)
    return(data.frame())
  }
  
  combinedScores <- dplyr::bind_rows(rows)
  return(combinedScores)
}

#' Merge ABN Inference Markov Blanket Performance Score Files
#'
#' Searches a results directory for `Scores_*_MB_*.RData` files, optionally filters
#' by simulation design parameters (`p`, `d`, `mp`, `node_dists`, `mb_size`, `target_node`),
#' loads each matching `Scores_<name_data>` object into an isolated environment, attaches
#' parsed metadata columns, and combines them into a single data frame.
#'
#' @param results_dir Character scalar giving the path to the directory
#'   containing `Scores_*.RData` files (typically `here::here("results")`).
#' @param p Optional numeric or integer vector of graph sizes (`p`) to keep.
#'   When `NULL` (default), all graph sizes are included.
#' @param d Optional numeric vector of target graph densities (`d`) to keep.
#'  When `NULL` (default), all densities are included.
#' @param mp Optional numeric or integer vector of maximum parent counts (`mp`)
#'   to keep. When `NULL` (default), all `mp` values are included.
#' @param node_dists Optional character vector of node distribution strings
#'   (e.g., `"balanced"`, `"G"`, `"P"`, `"B"`, `"M"`) to keep. When `NULL`
#'   (default), all node distribution settings are included.
#' @param target_node Character vector of target node names to include.
#'    When `NULL` (default), all target node names are included.
#' @param mb_size Integer vector of Markov blanket sizes to include.
#'    When `NULL` (default), all Markov blanket sizes are included.
#'
#' @return A `data.frame` combining all matching score tables across graphs and
#'   filter settings, with added columns `graph`, `p`, `d`, `mp`, `node_dists`,
#'   `mb_size`, and `target_node`
#'
#' @export
#' @importFrom dplyr bind_rows mutate
merge_mb_scores <- function(results_dir, p = NULL, d = NULL, mp = NULL, node_dists = NULL, target_node = NULL, mb_size = NULL) {
  files <- list.files(
    path = results_dir,
    pattern = "Scores_.*\\.Rdata",
    full.names = TRUE
  )
  files <- files[grepl("_MB_", basename(files), fixed = TRUE)] # keep only MB files
  
  if (length(files) == 0) {
    warning("No score files found in: ", results_dir)
    return(data.frame())
  }
  
  if (!is.null(d)){
    d_formatted <- sprintf("%02d", as.integer(d * 10))
  }
  
  rows <- list()
  idx <- 1L
  
  for (f in files) {
    fname <- tools::file_path_sans_ext(basename(f))
    name_data  <- sub("^[Ss]cores_", "", fname)
    
    # Extract design factors from filename when present
    p_match <- regmatches(name_data, regexpr("(?<=^|_)p_[0-9]+(?=_|$)", name_data, perl = TRUE))
    p_val   <- as.numeric(sub("^p_", "", p_match)) 
    
    d_match <- regmatches(name_data, regexpr("(?<=^|_)d_[0-9]+(?=_|$)", name_data, perl = TRUE))
    d_code  <- sub("^d_", "", d_match) 
    d_val   <- as.numeric(d_code) / 10 
    
    mp_match <- regmatches(name_data, regexpr("(?<=^|_)mp_[0-9]+(?=_|$)", name_data, perl = TRUE))
    mp_val   <- as.numeric(sub("^mp_", "", mp_match)) 
    
    nd_match <- regmatches(name_data, regexpr("(?<=^|_)nodedists_.+?(?=_(?:d_|mp_|graph)[0-9_]|$)", name_data, perl = TRUE))
    nd_val   <- sub("^nodedists_", "", nd_match) 
    
    mb_match <- regmatches(name_data, regexpr("(?<=^|_)MB_[0-9]+(?=_|$)", name_data, perl = TRUE))
    mb_val   <- as.numeric(sub("^MB_", "", mb_match))
    
    node_match  <- regmatches(name_data, regexpr("(?<=^|_)node_[^_]+(?=_|$)", name_data, perl = TRUE))
    node_val <- sub("^node", "", node_match) 
    
    g_match  <- regmatches(name_data, regexpr("graph[0-9]+$", name_data))
    graph_id <- sub("^graph", "", g_match) 
    
    if (!is.null(p) && !(p_val == p)) next
    if (!is.null(d) && !(d_code == d_formatted)) next
    if (!is.null(mp) && !(mp_val == mp)) next
    if (!is.null(node_dists) && !(nd_val == node_dists)) next
    if (!is.null(target_node) && !(node_val == target_node)) next
    if (!is.null(mb_size) && !(mb_val == mb_size)) next
    
    env <- new.env(parent = emptyenv())
    load(f, envir = env)
    scores_name <- grep(paste0("^Scores_", name_data, "$"), ls(envir = env), value = TRUE)
    if (length(scores_name) == 0) {
      stop("Object Scores_", name_data, " not found in file: ", f)
    }
    scores <- env[[scores_name[1]]]
    
    scores$graph <- graph_id
    scores$p <- p_val
    scores$d <- d_val
    scores$mp <- mp_val
    scores$node_dists <- nd_val
    scores$target_node <- node_val
    scores$mb_size <- mb_val
    
    rows[[idx]] <- scores
    idx <- idx + 1L
  }
  
  if (length(rows) == 0) {
    warning("No score files matched the specified filters in: ", results_dir)
    return(data.frame())
  }
  
  combinedScores <- dplyr::bind_rows(rows)
  return(combinedScores)
}

#' Plot ABN Inference Performance Across Distribution Settings
#'
#' Aggregates node-level evaluation metrics across simulation graphs and
#' generates a side-by-side comparison of prediction scores (panel A) and
#' classification accuracy for discrete nodes (panel B) across evidence sets and
#' node distribution configurations.
#'
#' @param scores A `data.frame` of combined evaluation scores
#'   (typically returned by [merge_scores()]), containing columns `p`, `d`,
#'   `mp`, `node_dists`, `node_type`, `evidence`, `graph`, `score`, and
#'   `accuracy`.
#' @param results_dir Character scalar giving the path to the results directory.
#'   Required when `save_figure = TRUE`, figures are written to
#'   `<results_dir>/figures`. Default is `NULL`.
#' @param save_figure Logical scalar. If `TRUE`, saves the combined figure as a
#'   PDF (`12 x 7.5` inches, `300` dpi) inside `<results_dir>/figures`. Default
#'   is `FALSE`.
#'
#'
#' @importFrom dplyr mutate group_by summarise filter
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point facet_wrap
#'   as_labeller expand_limits scale_x_discrete scale_y_continuous
#'   scale_color_manual scale_fill_manual theme_minimal theme element_text
#'   element_rect element_blank labs margin ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom stats sd
#' @export
#'
plot_compare_distributions <- function(scores, results_dir = NULL, save_figure= FALSE){
  p_val  <- unique(scores$p)[1]
  d_val  <- unique(scores$d)[1]
  mp_val <- unique(scores$mp)[1]
  
  dist_levels <- c("Balanced", "Gaussian", "Poisson", "Binomial", "Multinomial")
  evidence_labels <- c(
    "none"     = "No evidence",
    "children"  = "Children",
    "spouses" = "Spouses",
    "parents"  = "Parents",
    "mb"       = "Full Markov blanket"
  )
    
  scenario_colors <- c(
    "Gaussian"    = "#0072B2", # Bleu franc
    "Poisson"     = "#D55E00", # Vermillon
    "Binomial"    = "#009E73", # Vert bleuté
    "Multinomial" = "#CC79A7", # Violet rosé
    "Balanced"    = "#E69F00"  # Orange ambré
  )
    
  scores <- scores %>%
    mutate(
      evidence = factor(evidence, levels = c("none", "children", "spouses", "parents", "mb")),
      node_type = factor(node_type, levels = c("gaussian", "binomial", "poisson", "multinomial")),
      node_dists = dplyr::case_when(
        node_dists == "balanced" ~ "Balanced",
        node_dists == "G" ~ "Gaussian",
        node_dists == "P" ~ "Poisson",
        node_dists == "B" ~ "Binomial",
        node_dists == "M" ~ "Multinomial"
      ),
      node_dists = factor(node_dists, levels = dist_levels)
    )
  
  df_by_network <- scores %>%
    group_by(node_dists, node_type, evidence, graph) %>%
    summarise(network_mean_score = mean(score, na.rm = TRUE), .groups = 'drop')
    
  df_final_agg <- df_by_network %>%
    group_by(node_dists, node_type, evidence) %>%
    summarise(
      mean_score = mean(network_mean_score, na.rm = TRUE),
      se_score   = sd(network_mean_score, na.rm = TRUE) / sqrt(n()), 
      .groups    = 'drop'
    )
    
  plot_distrib_scores <- ggplot(df_final_agg, aes(x = evidence, y = mean_score, color = node_dists, group = node_dists)) +
    geom_ribbon(aes(ymin = mean_score - se_score, ymax = mean_score + se_score, fill = node_dists), 
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~ node_type, scales = "free_y", ncol = 2, labeller = as_labeller(c(
      gaussian = "Gaussian (MSE)",
      binomial = "Binomial (Brier score)",
      poisson = "Poisson (Deviance)",
      multinomial = "Multinomial (Brier score)"
    ))) +
    expand_limits(y = 0) +
    scale_x_discrete(
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    scale_color_manual(name = "Distribution setting",values = scenario_colors) +
    scale_fill_manual(name = "Distribution setting",values = scenario_colors) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )+
    labs(
      title = "A) Prediction scores",
      x = "",
      y = "Score (lower is better)"
    ) 
  
  acc_p_data <- scores %>%
    filter(grepl("binomial|multinomial", node_type))
  
  df_by_network_acc <- acc_p_data %>%
    group_by(node_dists, node_type, evidence, graph) %>%
    summarise(network_mean_score = mean(accuracy, na.rm = TRUE), .groups = 'drop')
  
  df_final_agg_acc <- df_by_network_acc %>%
    group_by(node_dists, node_type, evidence) %>%
    summarise(
      mean_score = mean(network_mean_score, na.rm = TRUE),
      # Calcul de l'erreur standard (sd/sqrt(n)) basée sur les 10 réseaux
      se_score   = sd(network_mean_score, na.rm = TRUE) / sqrt(n()), 
      .groups    = 'drop'
    )
    
  plot_distrib_accuracy <- ggplot(df_final_agg_acc, aes(x = evidence, y = mean_score, color = node_dists, group = node_dists)) +
    geom_ribbon(aes(ymin = mean_score - se_score, ymax = mean_score + se_score, fill = node_dists), 
                alpha = 0.15, color = NA) +
    
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~ node_type, scales = "free_y", ncol = 1,labeller = as_labeller(c(
      binomial = "Binomial",
      multinomial = "Multinomial"
    ))) +
      
    scale_x_discrete(
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(name = "Distribution setting",values = scenario_colors) +
    scale_fill_manual(name = "Distribution setting",values = scenario_colors) +
    
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1)
    ) +
    labs(
      title = "B) Classification accuracy",
      x = "",
      y = "Accuracy rate (higher is better)"
    ) 
  
  final_distrib_plot <- (plot_distrib_scores | plot_distrib_accuracy) +
    plot_layout(widths = c(1.8, 1), guides = "collect") +
    plot_annotation(
      title = paste0("Performance Evaluation - Effect of distribution settings (p = ",p_val, ", d = ", gsub('[.]','',d_val), ", mp = ",mp_val,")"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 4))
      )
    ) &
    theme(legend.position = "bottom")
  
  print(final_distrib_plot)
    
  if (save_figure){
    fig_dir <- file.path(results_dir, "figures")
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }
    ggsave(
      filename = paste0(fig_dir,"/scores_p_",p_val,"_d_",gsub('[.]','',d_val),"_mp_",mp_val,".pdf"),
      plot = final_distrib_plot,
      device = "pdf",
      width = 12,        # Width in inches
      height = 7.5,      # Height in inches
      units = "in",
      dpi = 300          # Sharp rendering resolution
    )
  }
}

#' Plot ABN Inference Performance Across Number of Nodes
#'
#' Aggregates node-level evaluation metrics across simulation graphs and
#' generates a side-by-side comparison of prediction scores (panel A) and
#' classification accuracy for discrete nodes (panel B) across evidence sets and
#' number of nodes.
#'
#' @param scores A `data.frame` of combined evaluation scores
#'   (typically returned by [merge_scores()]), containing columns `p`, `d`,
#'   `mp`, `node_dists`, `node_type`, `evidence`, `graph`, `score`, and
#'   `accuracy`.
#' @param results_dir Character scalar giving the path to the results directory.
#'   Required when `save_figure = TRUE`, figures are written to
#'   `<results_dir>/figures`. Default is `NULL`.
#' @param save_figure Logical scalar. If `TRUE`, saves the combined figure as a
#'   PDF (`12 x 7.5` inches, `300` dpi) inside `<results_dir>/figures`. Default
#'   is `FALSE`.
#'
#'
#' @importFrom dplyr mutate group_by summarise filter
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point facet_wrap
#'   as_labeller expand_limits scale_x_discrete scale_y_continuous
#'   scale_color_manual scale_fill_manual theme_minimal theme element_text
#'   element_rect element_blank labs margin ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom stats sd
#' @importFrom ggh4x facetted_pos_scales
#' @export
#'
plot_compare_p <- function(scores, results_dir = NULL, save_figure = FALSE){
  d_val  <- unique(scores$d)[1]
  mp_val <- unique(scores$mp)[1]
  
  dist_levels <- c("Balanced", "Gaussian", "Poisson", "Binomial", "Multinomial")
  evidence_labels <- c(
    "none"     = "No evidence",
    "children"  = "Children",
    "spouses" = "Spouses",
    "parents"  = "Parents",
    "mb"       = "Full Markov blanket"
  )
  evidence_colors <- c(
    "none"     = "#999999",
    "children"  = "#0072B2",
    "spouses" = "#E69F00",
    "parents"  = "#009E73",
    "mb"       = "#D55E00"
  )
  
  scores <- scores %>%
    mutate(
      p = factor(p, levels = sort(unique(as.numeric(as.character(p))))),
      evidence = factor(evidence, levels = c("none", "children", "spouses", "parents", "mb")),
      node_type = factor(node_type, levels = c("gaussian", "binomial", "poisson", "multinomial")),
      node_dists = dplyr::case_when(
        node_dists == "balanced" ~ "Balanced",
        node_dists == "G" ~ "Gaussian",
        node_dists == "P" ~ "Poisson",
        node_dists == "B" ~ "Binomial",
        node_dists == "M" ~ "Multinomial"
      ),
      node_dists = factor(node_dists, levels = dist_levels)
    )
  
  panel_caps <- scores %>%
    group_by(node_type, p, evidence) %>%
    summarise(
      upper_whisker = tryCatch(boxplot.stats(score)$stats[5], error = function(e) max(score, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    group_by(node_type) %>%
    summarise(
      # Take the maximum whisker across all scenarios and methods, plus 20% breathing room
      ymax = max(upper_whisker, na.rm = TRUE) * 1.2
    ) %>%
    arrange(node_type) # Ensures exact alignment with ggplot's alphabetical panel order
    
  panel_scales <- purrr::map2(panel_caps$node_type, panel_caps$ymax, function(nt, mx) {
    scale_y_continuous(limits = c(0, mx), oob = scales::squish, trans="pseudo_log")
  })
    
  # Plot A: Scores vs p
  plot_p_scores <- ggplot(scores, aes(x = p, y = score, fill = evidence)) +
    geom_boxplot(outlier.alpha = 0.3, linewidth = 0.3, position = position_dodge(0.8)) +
      facet_wrap(~ node_type, scales = "free_y", ncol = 2, labeller = as_labeller(c(
        gaussian = "Gaussian (MSE)",
        binomial = "Binomial (Brier score)",
        poisson = "Poisson (Deviance)",
        multinomial = "Multinomial (Brier score)"
      ))) +
      facetted_pos_scales(y = panel_scales) +
      scale_fill_manual(
        values = evidence_colors,
        breaks = names(evidence_labels),
        labels = evidence_labels
      ) +
      theme_minimal(base_size = 12) +
      theme(
        strip.text = element_text(face = "bold", size = 11),
        panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = "A) Prediction scores",
        x = "Number of nodes (p)",
        y = "Score (lower is better)"
      ) 
    
  # Plot B: Accuracy vs p for Discrete Nodes
  acc_p_data <- scores %>%
    filter(grepl("binomial|multinomial", node_type))
    
  plot_p_accuracy <- ggplot(acc_p_data, aes(x = p, y = accuracy, fill = evidence)) +
    geom_boxplot(outlier.alpha = 0.3, linewidth = 0.3, position = position_dodge(0.8)) +
    facet_wrap(~ node_type, ncol = 1,labeller = as_labeller(c(
      binomial = "Binomial",
      multinomial = "Multinomial"
    ))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(
      values = evidence_colors,
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "B) Classification accuracy",
      x = "Number of nodes (p)",
      y = "Accuracy rate (higher is better)"
    ) 
  
  # Combine using Patchwork
  final_p_plot <- (plot_p_scores | plot_p_accuracy) +
    plot_layout(widths = c(1.8, 1), guides = "collect") +
    plot_annotation(
      title = paste0("Performance Evaluation - Effect of Graph Size (p) (d = ", d_val, ", mp = ",mp_val,")"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 4))
      )
    ) &
    theme(legend.position = "bottom")
  
  print(final_p_plot)
    
  if (save_figure){
    fig_dir <- file.path(results_dir, "figures")
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }
    ggsave(
      filename = paste0(fig_dir,"/scores_p_d_",gsub('[.]','',d_val),"_mp_",mp_val,".pdf"),
      plot = final_p_plot,
      device = "pdf",
      width = 12,        # Width in inches
      height = 7.5,      # Height in inches
      units = "in",
      dpi = 300          # Sharp rendering resolution
    )
  }
}

#' Plot ABN Inference Performance Across Densities
#'
#' Aggregates node-level evaluation metrics across simulation graphs and
#' generates a side-by-side comparison of prediction scores (panel A) and
#' classification accuracy for discrete nodes (panel B) across evidence sets and
#' densities.
#'
#' @param scores A `data.frame` of combined evaluation scores
#'   (typically returned by [merge_scores()]), containing columns `p`, `d`,
#'   `mp`, `node_dists`, `node_type`, `evidence`, `graph`, `score`, and
#'   `accuracy`.
#' @param results_dir Character scalar giving the path to the results directory.
#'   Required when `save_figure = TRUE`, figures are written to
#'   `<results_dir>/figures`. Default is `NULL`.
#' @param save_figure Logical scalar. If `TRUE`, saves the combined figure as a
#'   PDF (`12 x 7.5` inches, `300` dpi) inside `<results_dir>/figures`. Default
#'   is `FALSE`.
#'
#'
#' @importFrom dplyr mutate group_by summarise filter
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point facet_wrap
#'   as_labeller expand_limits scale_x_discrete scale_y_continuous
#'   scale_color_manual scale_fill_manual theme_minimal theme element_text
#'   element_rect element_blank labs margin ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom stats sd
#' @importFrom ggh4x facetted_pos_scales
#' @export
#'
plot_compare_d <- function(scores, results_dir = NULL, save_figure = FALSE){
  p_val  <- unique(scores$p)[1]
  mp_val <- unique(scores$mp)[1]
  
  dist_levels <- c("Balanced", "Gaussian", "Poisson", "Binomial", "Multinomial")
  evidence_labels <- c(
    "none"     = "No evidence",
    "children"  = "Children",
    "spouses" = "Spouses",
    "parents"  = "Parents",
    "mb"       = "Full Markov blanket"
  )
  evidence_colors <- c(
    "none"     = "#999999",
    "children"  = "#0072B2",
    "spouses" = "#E69F00",
    "parents"  = "#009E73",
    "mb"       = "#D55E00"
  )
  
  scores <- scores %>%
    mutate(
      d = factor(d, levels = sort(unique(as.numeric(as.character(d))))),
      evidence = factor(evidence, levels = c("none", "children", "spouses", "parents", "mb")),
      node_type = factor(node_type, levels = c("gaussian", "binomial", "poisson", "multinomial")),
      node_dists = dplyr::case_when(
        node_dists == "balanced" ~ "Balanced",
        node_dists == "G" ~ "Gaussian",
        node_dists == "P" ~ "Poisson",
        node_dists == "B" ~ "Binomial",
        node_dists == "M" ~ "Multinomial"
      ),
      node_dists = factor(node_dists, levels = dist_levels)
    )
  
  panel_caps <- scores %>%
    group_by(node_type, d, evidence) %>%
    summarise(
      upper_whisker = tryCatch(boxplot.stats(score)$stats[5], error = function(e) max(score, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    group_by(node_type) %>%
    summarise(
      # Take the maximum whisker across all scenarios and methods, plus 20% breathing room
      ymax = max(upper_whisker, na.rm = TRUE) * 1.2
    ) %>%
    arrange(node_type) # Ensures exact alignment with ggplot's alphabetical panel order
  
  panel_scales <- purrr::map2(panel_caps$node_type, panel_caps$ymax, function(nt, mx) {
    scale_y_continuous(limits = c(0, mx), oob = scales::squish, trans="pseudo_log")
  })
  
  # Plot A: Scores vs d
  plot_d_scores <- ggplot(scores, aes(x = d, y = score, fill = evidence)) +
    geom_boxplot(outlier.alpha = 0.3, linewidth = 0.3, position = position_dodge(0.8)) +
    facet_wrap(~ node_type, scales = "free_y", ncol = 2, labeller = as_labeller(c(
      gaussian = "Gaussian (MSE)",
      binomial = "Binomial (Brier score)",
      poisson = "Poisson (Deviance)",
      multinomial = "Multinomial (Brier score)"
    ))) +
    facetted_pos_scales(y = panel_scales) +
    scale_fill_manual(
      values = evidence_colors,
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "A) Prediction scores",
      x = "Densities (d)",
      y = "Score (lower is better)"
    ) 
  
  # Plot B: Accuracy vs d for Discrete Nodes
  acc_d_data <- scores %>%
    filter(grepl("binomial|multinomial", node_type))
  
  plot_d_accuracy <- ggplot(acc_d_data, aes(x = d, y = accuracy, fill = evidence)) +
    geom_boxplot(outlier.alpha = 0.3, linewidth = 0.3, position = position_dodge(0.8)) +
    facet_wrap(~ node_type, ncol = 1,labeller = as_labeller(c(
      binomial = "Binomial",
      multinomial = "Multinomial"
    ))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(
      values = evidence_colors,
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "B) Classification accuracy",
      x = "Densities (d)",
      y = "Accuracy rate (higher is better)"
    ) 
  
  # Combine using Patchwork
  final_d_plot <- (plot_d_scores | plot_d_accuracy) +
    plot_layout(widths = c(1.8, 1), guides = "collect") +
    plot_annotation(
      title = paste0("Performance Evaluation - Effect of Densities (d) (p = ", p_val, ", mp = ",mp_val,")"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 4))
      )
    ) &
    theme(legend.position = "bottom")
  
  print(final_d_plot)
  
  if (save_figure){
    fig_dir <- file.path(results_dir, "figures")
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }
    ggsave(
      filename = paste0(fig_dir,"/scores_p_",p_val,"_d_mp_",mp_val,".pdf"),
      plot = final_d_plot,
      device = "pdf",
      width = 12,        # Width in inches
      height = 7.5,      # Height in inches
      units = "in",
      dpi = 300          # Sharp rendering resolution
    )
  }
}

#' Plot ABN Inference Performance Across Maximal Number of Parents
#'
#' Aggregates node-level evaluation metrics across simulation graphs and
#' generates a side-by-side comparison of prediction scores (panel A) and
#' classification accuracy for discrete nodes (panel B) across evidence sets and
#' maximal number of parents.
#'
#' @param scores A `data.frame` of combined evaluation scores
#'   (typically returned by [merge_scores()]), containing columns `p`, `d`,
#'   `mp`, `node_dists`, `node_type`, `evidence`, `graph`, `score`, and
#'   `accuracy`.
#' @param results_dir Character scalar giving the path to the results directory.
#'   Required when `save_figure = TRUE`, figures are written to
#'   `<results_dir>/figures`. Default is `NULL`.
#' @param save_figure Logical scalar. If `TRUE`, saves the combined figure as a
#'   PDF (`12 x 7.5` inches, `300` dpi) inside `<results_dir>/figures`. Default
#'   is `FALSE`.
#'
#'
#' @importFrom dplyr mutate group_by summarise filter
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point facet_wrap
#'   as_labeller expand_limits scale_x_discrete scale_y_continuous
#'   scale_color_manual scale_fill_manual theme_minimal theme element_text
#'   element_rect element_blank labs margin ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom stats sd
#' @importFrom ggh4x facetted_pos_scales
#' @export
#'
plot_compare_mp <- function(scores, results_dir = NULL, save_figure = FALSE){
  p_val  <- unique(scores$p)[1]
  d_val <- unique(scores$d)[1]
  
  evidence_labels <- c(
    "none"     = "No evidence",
    "children"  = "Children",
    "spouses" = "Spouses",
    "parents"  = "Parents",
    "mb"       = "Full Markov blanket"
  )
  evidence_colors <- c(
    "none"     = "#999999",
    "children"  = "#0072B2",
    "spouses" = "#E69F00",
    "parents"  = "#009E73",
    "mb"       = "#D55E00"
  )
  
  id_cols <- c("evidence", "node", "node_type", "graph", "p", "d", "node_dists")
  id_cols <- intersect(id_cols, names(scores))
  
  scores <- scores %>%
    mutate(
      mp = as.numeric(as.character(mp)),
      evidence = factor(evidence, levels = c("none", "children", "spouses", "parents", "mb")),
      node_type = factor(node_type, levels = c("gaussian", "binomial", "poisson", "multinomial")),
      node_dists = dplyr::case_when(
        node_dists == "balanced" ~ "Balanced",
        node_dists == "G" ~ "Gaussian",
        node_dists == "P" ~ "Poisson",
        node_dists == "B" ~ "Binomial",
        node_dists == "M" ~ "Multinomial"
      ),
      node_dists = factor(node_dists, levels = dist_levels)
    )
  
  ref_scores <- scores %>%
    filter(mp == p_val-1) %>%
    select(all_of(id_cols),ref_score = score, ref_accuracy = accuracy) %>%
    distinct()
  
  delta_scores <- scores %>%
    filter(mp != p_val-1) %>%
    inner_join(ref_scores, by = id_cols) %>%
    mutate(
      delta_score = score - ref_score,
      delta_accuracy = accuracy - ref_accuracy,
      mp = factor(mp, levels = sort(unique(mp)))
    )
  
  sum_scores <- delta_scores %>%
    group_by(node_type, mp, evidence) %>%
    summarise(
      med = median(delta_score, na.rm = TRUE),
      q25 = quantile(delta_score, 0.25, na.rm = TRUE),
      q75 = quantile(delta_score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  panel_caps <- sum_scores %>%
    group_by(node_type) %>%
    summarise(
      bnd = max(abs(c(q25, q75, med)), na.rm = TRUE) * 1.4,
      .groups = "drop"
    ) %>%
    arrange(node_type)
  
  panel_scales <- purrr::map2(panel_caps$node_type, panel_caps$bnd, function(nt, bnd) {
    bnd <- ifelse(is.finite(bnd) && bnd > 0, bnd, 0.05)
    scale_y_continuous(limits = c(-bnd, bnd), oob = scales::squish)
  })
  
  pd <- position_dodge(width = 0.65)
  
  # Plot A: Scores vs mp
  plot_mp_scores <- ggplot(sum_scores, aes(x = mp, y = med, color = evidence)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.6) +
    geom_linerange(aes(ymin = q25, ymax = q75), position = pd, linewidth = 0.8) +
    geom_point(position = pd, size = 2) +
    facet_wrap(~ node_type, scales = "free_y", ncol = 2, labeller = as_labeller(c(
      gaussian = "Gaussian (\u0394 MSE)",
      binomial = "Binomial (\u0394 Brier score)",
      poisson = "Poisson (\u0394 Deviance)",
      multinomial = "Multinomial (\u0394 Brier score)"
    ))) +
    facetted_pos_scales(y = panel_scales) +
    scale_color_manual(
      values = evidence_colors,
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "A) \u0394 Prediction scores",
      x = "Maximal number of parents (mp)",
      y = "Score (0 = same as not restricted)"
    ) 
  
  # Plot B: Accuracy vs d for Discrete Nodes
  sum_acc <- delta_scores %>%
    filter(node_type %in% c("binomial", "multinomial")) %>%
    group_by(node_type, mp, evidence) %>%
    summarise(
      med = median(delta_accuracy, na.rm = TRUE),
      q25 = quantile(delta_accuracy, 0.25, na.rm = TRUE),
      q75 = quantile(delta_accuracy, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  acc_bnd <- sum_acc %>%
    summarise(b = max(abs(c(q25, q75, med)), na.rm = TRUE) * 1.4) %>%
    pull(b)
  acc_bnd <- ifelse(length(acc_bnd) == 0 || !is.finite(acc_bnd) || acc_bnd == 0, 0.05, acc_bnd)
  
  plot_mp_accuracy <- ggplot(sum_acc, aes(x = mp, y = med, color = evidence)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.6) +
    geom_linerange(aes(ymin = q25, ymax = q75), position = pd, linewidth = 0.8) +
    geom_point(position = pd, size = 2) +
    facet_wrap(~ node_type, ncol = 1,labeller = as_labeller(c(
      binomial = "Binomial (\u0394 Accuracy)",
      multinomial = "Multinomial (\u0394 Accuracy)"
    ))) +
    scale_y_continuous(limits = c(-acc_bnd, acc_bnd), oob = scales::squish) +
    scale_color_manual(
      values = evidence_colors,
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "B) \u0394 Classification accuracy",
      x = "Maximal number of parents (mp)",
      y = "Accuracy rate (0 = same as not restricted)"
    ) 
  
  # Combine using Patchwork
  final_mp_plot <- (plot_mp_scores | plot_mp_accuracy) +
    plot_layout(widths = c(1.8, 1), guides = "collect") +
    plot_annotation(
      title = paste0("Performance Evaluation - Effect of the Maximal Number of Parents (p = ", p_val, ", d = ",d_val,")"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 4))
      )
    ) &
    theme(legend.position = "bottom")
  
  print(final_mp_plot)
  
  if (save_figure){
    fig_dir <- file.path(results_dir, "figures")
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }
    has_cairo <- isTRUE(capabilities("cairo"))
    ggsave(
      filename = paste0(fig_dir,"/scores_p_",p_val,"_d_",gsub('[.]','',d_val),"_mp.pdf"),
      plot = final_mp_plot,
      device = if (has_cairo) grDevices::cairo_pdf else "pdf",
      width = 12,        # Width in inches
      height = 7.5,      # Height in inches
      units = "in",
      dpi = 300          # Sharp rendering resolution
    )
  }
}

#' Plot ABN Inference Computational Time Across Maximal Number of Parents
#'
#' Aggregates node-level evaluation metrics across simulation graphs and
#' generates a comparison of computational time across evidence sets and
#' maximal number of parents.
#'
#' @param scores A `data.frame` of combined evaluation scores
#'   (typically returned by [merge_scores()]), containing columns `p`, `d`,
#'   `mp`, `node_dists`, `node_type`, `evidence`, `graph`, `score`, and
#'   `accuracy`.
#' @param results_dir Character scalar giving the path to the results directory.
#'   Required when `save_figure = TRUE`, figures are written to
#'   `<results_dir>/figures`. Default is `NULL`.
#' @param save_figure Logical scalar. If `TRUE`, saves the combined figure as a
#'   PDF (`12 x 7.5` inches, `300` dpi) inside `<results_dir>/figures`. Default
#'   is `FALSE`.
#' @param log_y Logical scalar. If `TRUE`, y-scale is log-transformed.
#'
#' @importFrom dplyr mutate group_by summarise filter
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point facet_wrap
#'   as_labeller expand_limits scale_x_discrete scale_y_continuous
#'   scale_color_manual scale_fill_manual theme_minimal theme element_text
#'   element_rect element_blank labs margin ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom stats sd
#' @importFrom ggh4x facetted_pos_scales
#' @export
#'
plot_compare_mp_time <- function(scores, results_dir = NULL, save_figure = FALSE, log_y = FALSE){
  p_val  <- unique(scores$p)[1]
  d_val <- unique(scores$d)[1]
  
  evidence_labels <- c(
    "none"     = "No evidence",
    "children"  = "Children",
    "spouses" = "Spouses",
    "parents"  = "Parents",
    "mb"       = "Full Markov blanket"
  )
  evidence_colors <- c(
    "none"     = "#999999",
    "children"  = "#0072B2",
    "spouses" = "#E69F00",
    "parents"  = "#009E73",
    "mb"       = "#D55E00"
  )
  
  scores <- scores %>%
    mutate(
      mp = as.numeric(as.character(mp)),
      evidence = factor(evidence, levels = c("none", "children", "spouses", "parents", "mb")),
      node_type = factor(node_type, levels = c("gaussian", "binomial", "poisson", "multinomial")),
      node_dists = dplyr::case_when(
        node_dists == "balanced" ~ "Balanced",
        node_dists == "G" ~ "Gaussian",
        node_dists == "P" ~ "Poisson",
        node_dists == "B" ~ "Binomial",
        node_dists == "M" ~ "Multinomial"
      ),
      node_dists = factor(node_dists, levels = dist_levels)
    )
  
  time_df <- scores %>% select("evidence", "comp_time", "graph", "mp", "node_dists") %>%
    distinct() %>%
    group_by(node_dists, mp, evidence) %>%
    summarise(
      med = median(comp_time, na.rm = TRUE),
      q25 = quantile(comp_time, 0.25, na.rm = TRUE),
      q75 = quantile(comp_time, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(mp = factor(mp, levels = sort(unique(mp))))
  
  p_time <- ggplot(time_df, aes(x = mp, y = med, color = evidence, group = evidence)) +
    geom_line(position = position_identity(), linewidth = 0.75, alpha = 0.9) +
    geom_linerange(aes(ymin = q25, ymax = q75), position = position_identity(), linewidth = 0.8) +
    geom_point(position = position_identity(), size = 2.1) +
    facet_wrap(~ node_dists, nrow = 1) +
    scale_color_manual(
      values = evidence_colors,
      breaks = names(evidence_labels),
      labels = evidence_labels
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0("Computational Time - Effect of the Maximal Number of Parents (p = ", p_val, ", d = ", d_val, ")"),
      x = "Maximal number of parents (mp)",
      y = if (log_y) "Computational time (s, log scale)" else "Computational time (s)"
    )
  
  if (log_y) {
    p_time <- p_time + scale_y_log10()
  }
  
  print(p_time)
  
  if (save_figure){
    fig_dir <- file.path(results_dir, "figures")
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }
    ggsave(
      filename = paste0(fig_dir,"/runtime_p_",p_val,"_d_",gsub('[.]','',d_val),"_mp.pdf"),
      plot = p_time,
      device = "pdf",
      width = 12,        # Width in inches
      height = 7.5,      # Height in inches
      units = "in",
      dpi = 300          # Sharp rendering resolution
    )
  }
}

#' Plot ABN Inference Performance Across Maximal Number of Parents
#'
#' Aggregates node-level evaluation metrics across simulation graphs and
#' generates a side-by-side comparison of prediction scores (panel A) and
#' classification accuracy for discrete nodes (panel B) across evidence sets and
#' maximal number of parents.
#'
#' @param scores A `data.frame` of combined evaluation scores
#'   (typically returned by [merge_scores()]), containing columns `p`, `d`,
#'   `mp`, `node_dists`, `node_type`, `evidence`, `graph`, `score`, and
#'   `accuracy`.
#' @param results_dir Character scalar giving the path to the results directory.
#'   Required when `save_figure = TRUE`, figures are written to
#'   `<results_dir>/figures`. Default is `NULL`.
#' @param save_figure Logical scalar. If `TRUE`, saves the combined figure as a
#'   PDF (`12 x 7.5` inches, `300` dpi) inside `<results_dir>/figures`. Default
#'   is `FALSE`.
#'
#'
#' @importFrom dplyr mutate group_by summarise filter
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point facet_wrap
#'   as_labeller expand_limits scale_x_discrete scale_y_continuous
#'   scale_color_manual scale_fill_manual theme_minimal theme element_text
#'   element_rect element_blank labs margin ggsave
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom stats sd
#' @importFrom ggh4x facetted_pos_scales
#' @export
#'
plot_compare_mb <- function(scores, results_dir = NULL, save_figure = FALSE){
  p_val  <- unique(scores$p)[1]
  d_val <- unique(scores$d)[1]
  mp_val <- unique(scores$mp)[1]
  
  dist_levels <- c("Balanced", "Gaussian", "Poisson", "Binomial", "Multinomial")
  evidence_labels <- c(
    "none"     = "No evidence",
    "children"  = "Children",
    "spouses" = "Spouses",
    "parents"  = "Parents",
    "mb"       = "Full Markov blanket"
  )
  
  scores <- scores %>%
    mutate(
      mb_size = factor(mb_size, levels = sort(unique(as.numeric(as.character(mb_size))))),
      target_node = factor(target_node, levels = sort(unique(as.character(target_node)))),
      evidence = factor(evidence, levels = c("none", "children", "spouses", "parents", "mb")),
      node_type = factor(node_type, levels = c("gaussian", "binomial", "poisson", "multinomial")),
      node_dists = dplyr::case_when(
        node_dists == "balanced" ~ "Balanced",
        node_dists == "G" ~ "Gaussian",
        node_dists == "P" ~ "Poisson",
        node_dists == "B" ~ "Binomial",
        node_dists == "M" ~ "Multinomial"
      ),
      node_dists = factor(node_dists, levels = dist_levels)
    )
  
  scores_filtered <- scores %>% 
    filter(as.character(node) == as.character(target_node))
  
  delta_data <- scores_filtered %>%
    filter(evidence %in% c("none", "mb")) %>%
    group_by(graph, mb_size, node_type, target_node, evidence) %>%
    summarise(Scores = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = evidence, values_from = Scores) %>%
    mutate(
      # Score metrics (MSE, Deviance, Brier): higher score is worse, so positive delta means 'mb' is better
      Score_Diff = none - mb
    )
  
  delta_summary <- delta_data %>%
    group_by(node_type, mb_size) %>%
    summarise(
      mean_diff = mean(Score_Diff, na.rm = TRUE),
      median_diff = median(Score_Diff, na.rm = TRUE),
      .groups = "drop"
    )
  
  plot_delta <- ggplot(delta_data, aes(x = mb_size, y = Score_Diff)) +
    # Boxplots showing variance across graphs for each MB size
    geom_boxplot(fill = "#8DA0CB", color = "#333333", alpha = 0.5, outlier.alpha = 0.3) +
    
    # Trend line connecting medians to highlight the trajectory across MB sizes
    geom_line(data = delta_summary, aes(x = mb_size, y = median_diff, group = 1), 
              color = "#2b5c8f", linewidth = 1.1) +
    geom_point(data = delta_summary, aes(x = mb_size, y = median_diff), 
               color = "#2b5c8f", size = 2) +
    
    facet_wrap(~ node_type, scales = "free_y", ncol = 2, labeller = as_labeller(c(
      gaussian = "Gaussian (MSE)",
      binomial = "Binomial (Brier score)",
      poisson = "Poisson (Deviance)",
      multinomial = "Multinomial (Brier score)"
    ))) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8)
    ) +
    labs(
      title = "A) Prediction scores",
      x = "Markov blanket size",
      y = "Delta score (higher means MB advantage)"
    ) 
  
  delta_acc_data <- scores_filtered %>%
    filter(node_type %in% c("binomial", "multinomial")) %>%
    filter(evidence %in% c("none", "mb")) %>%
    group_by(graph, mb_size, node_type, target_node, evidence) %>%
    summarise(Accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = evidence, values_from = Accuracy) %>%
    mutate(
      # Accuracy: higher is better, so mb - none means positive values favor mb
      Acc_Diff = mb - none
    )
  
  delta_acc_summary <- delta_acc_data %>%
    group_by(node_type, mb_size) %>%
    summarise(
      mean_diff = mean(Acc_Diff, na.rm = TRUE),
      median_diff = median(Acc_Diff, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3. Plot the performance gap (Delta) in Accuracy across MB sizes
  plot_delta_acc <- ggplot(delta_acc_data, aes(x = mb_size, y = Acc_Diff)) +
    # Boxplots showing variance across graphs for each MB size
    geom_boxplot(fill = "#FC8D62", color = "#333333", alpha = 0.5, outlier.alpha = 0.3) +
    
    # Trend line connecting medians to highlight the trajectory
    geom_line(data = delta_acc_summary, aes(x = mb_size, y = median_diff, group = 1), 
              color = "#d95f02", linewidth = 1.1) +
    geom_point(data = delta_acc_summary, aes(x = mb_size, y = median_diff), 
               color = "#d95f02", size = 2) +
    
    facet_wrap(~ node_type, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8)
    ) +
    labs(
      title = "B) Classification accuracy",
      x = "Markov blanket size",
      y = "Delta score (higher means MB advantage)"
    )
  
  # Combine using Patchwork
  final_mb_plot <- (plot_delta | plot_delta_acc) +
    plot_layout(widths = c(1.8, 1), guides = "collect") +
    plot_annotation(
      title = paste0("Performance Advantage of Markov Blanket (Score_none - Score_mb) (p = ", p_val,", d = ",d_val,", mp = ",mp_val,")"),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 4))
      )
    ) &
    theme(legend.position = "bottom")
  
  print(final_mb_plot)
  
  if (save_figure){
    fig_dir <- file.path(results_dir, "figures")
    if (!dir.exists(fig_dir)) {
      dir.create(fig_dir, recursive = TRUE)
    }
    ggsave(
      filename = paste0(fig_dir,"/scores_p_",p_val,"_d_",gsub('[.]','',d_val),"_mp_",mp_val,"_MB.pdf"),
      plot = final_mb_plot,
      device = "pdf",
      width = 12,        # Width in inches
      height = 7.5,      # Height in inches
      units = "in",
      dpi = 300          # Sharp rendering resolution
    )
  }
}