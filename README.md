# Additive Bayesian Networks Inference Simulation Study

The `abninferencestudy` package is designed to facilitate the simulation of data and evaluate the predictive performance of Additive Bayesian Networks (ABNs) inference in a controlled environment.

The primary goal is to evaluate network inference and conditional predictions against a known ground truth without requiring structural learning steps of ABN.

Simulated datasets are stored in the `data` directory, while results and generated plots from the simulation study are stored in the `results` directory.

## Installation

You can install the development version of `abninferencestudy` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("magalichampion/abninferencestudy")
```

## Experimental Design

The package implements a comprehensive simulation study framework to evaluate the performance of `predictAbn()` function under various data characteristics (number of nodes, densities, distributions) and evidence conditioning scenarios. 

## Workflow and Usage

The simulation study is organized into three main stages:

1. Data Generation (`generateSimData.Rmd`)

Simulates synthetic data from predefined ABN structures and regression parameters, saving each experimental condition into one `.Rdata` file.

  - File Naming Convention: `n_<sample_size>_p_<number_of_nodes>_nodedists_<node_distribution_type>_d_<density_without_decimal>_mp_<maximal_number_of_parents>_graph<number_of_DAG>.Rdata`

  - Parameter Explanation:

    - `<sample_size>`: Total number of observations.
    
    - `<number_of_nodes>`: Number of variables/nodes in the network ($p$).
    
    - `<node_distribution_type>`: Node distribution strategy:
    
        - `Balanced`: all distributions appear equally often.
        - `G`, `P`, `B`, `M`: One node of each distribution type with the remainder set to Gaussian, Poisson, Binomial, or Multinomial, respectively.
  
    - `<density_without_decimal>`: Structural graph density (without decimal point, e.g, `01` for 0.1).
    
    - `<number_of_DAG>`: Index of the replicate DAG (1 to 10).

  - Loaded dataset objects:
  
```r
# Load a dataset with 100 observations, 8 nodes, balanced distributions, 60% density
load("data/n_100_p_8_nodedists_Balanced_d_06_graph2.rda")
```
Each dataset object contains:

  - `data`: data.frame of observed node variables.

  - `DAG`: igraph object of the true network structure.

  - `adj.matrix`: Adjacency matrix representation.

  - `dists`: Named list of node distribution types.

  - `mycoefs`: List of true regression coefficients for each node's linear predictor.

  - `max.parents`: Maximum parent constraints per node.

  - `mse`: Conditional residual variances for the Gaussian nodes.

2. Inference and performance evaluation function (`RunABN_inference_benchmark.Rmd`)

Because raw sample-level predictions are memory-intensive, inference (`predictAbn()`) and evaluation (`EvaluatePerf()`) are executed together in a single parallelized pipeline without storing intermediate predictions.

- Core functions:

  - `RunABN_inference_benchmark(dataset, evidence_type)`: Iterates across every node in the network, isolates the requested conditional evidence set, runs `predictAbn()`, and immediately evaluates predictions against ground truth.

  - `EvaluatePerf(observed, predicted, distribution)`: Computes node-level predictive performance metrics tailored to each distribution family.

- Supported Evidence Scenarios (`evidence_type`):

  - `none`: Unconditional prior marginal predictions (no evidence).

  - `parents`: Condition on the target node's direct parents.

  - `children`: Condition on the target node's direct children.

  - `spouses`: Condition on the target node's co-parents (spouses).

  - `mb`: Condition on the target node's complete Markov blanket.
  
- Metrics comparing predictions against known ground-truth values across distribution families:

  - Gaussian: Mean Squared Error (MSE).

  - Poisson: Deviance ratio (dev_model / dev_null).

  - Binomial / Multinomial: Brier Score and Accuracy.

  - General metrics: Prediction bias (PredMean - ObsMean).

- Output storage: Evaluation summaries are aggregated across all 10 replicate DAGs and stored in the `results/` directory as: `Scores_n_<sample_size>_p_<number_of_nodes>_nodedists_<node_distribution_type>_d_<density_without_decimal>_mp_<maximal_number_of_parents>.Rdata`.

3. Results Analysis and Visualization (`abnplotinferenceresults.Rmd`)

Aggregates evaluation metrics across all experimental conditions and produces publication-ready figures stored in `results/`.

- Main Functions:
  
  - `merge_simulation_results()`: Combines all `Scores_*.Rdata` files from `results/` into a unified master data frame.
  
  - `plot_compare_distributions()`: Evaluates performance differences across node distribution strategies (`Balanced`, `G`, `P`, `B`, `M`).
  
  - `plot_compare_p()`, `plot_compare_d()`, and `plot_compare_mp()`: Visualizes scaling behavior as graph density, node count ($p$), and maximal number of parents ($mp$) increase.

  - `plot_diff_MB()`: Compares predictive scores with no evidence against conditioning on the full Markov blanket as the Markov blanket increases.
  
## Reproducing the study


To run the complete analysis pipeline, execute the vignettes in this order:

  - `generateSimData.Rmd`: Generate all synthetic DAGs, parameters, and simulated datasets in `data/`.

  - `runInferenceBenchmark.Rmd`: Run conditional inference (`predictAbn`) and compute performance metrics (`EvaluatePerf`) saved to `results/`.

  - `abnplotinferenceresults.Rmd`: Merge results and generate all final comparison plots.

⚠️ Computational Requirements:

  - Time: Full study requires several days of computation
    
  - Memory: Recommend 16+ GB RAM for larger datasets
    
  - Cores: Parallelization recommended (8+ cores optimal)
  
  - Quick Start: Set DEBUG <- TRUE in the vignettes to test the full pipeline on a reduced parameter subset.
    