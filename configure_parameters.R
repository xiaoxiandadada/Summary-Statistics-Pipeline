#!/usr/bin/env Rscript

# Minimal parameter configuration to satisfy lavaknock_pipeline.R sourcing.

get_default_params <- function() {
  list(
    data_dir = "./",
    output_dir = "./results",
    sample_size = 20000,
    window_size = 100000,
    window_step = 50000,
    knockoffs = 5,
    fdr = 0.1,
    ld_threshold = 0.75,
    prune_threshold = 99,
    threads = 4,
    run_mode = "full",
    use_example_data = FALSE,
    verbose = FALSE
  )
}

