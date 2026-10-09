# 02_run_optimize_SLIDE.R

# Safer SLIDE optimization retry for pILC pleomorphic-component (%) outcome

USER_LIB <- "/ix1/alee/LO_LAB/Personal/Rahul/R_library"
.libPaths(c(USER_LIB, .libPaths()))

SLIDE_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/SLIDE"
INPUT_YAML <- file.path(SLIDE_DIR, "Input_pleomorphic", "slide_pleomorphic.yaml")
RESULTS_DIR <- file.path(SLIDE_DIR, "Results_pleomorphic_stable_retry")
LOG_DIR <- file.path(SLIDE_DIR, "Logs_pleomorphic")
RETRY_YAML <- file.path(SLIDE_DIR, "Input_pleomorphic", "slide_pleomorphic_stable_retry.yaml")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

params <- yaml::read_yaml(INPUT_YAML)

# Write results separately; do not mix them with the failed 5,000-feature run.
params$out_path <- RESULTS_DIR

# Conservative tuning for the 24-sample pILC cohort.
# delta: higher value limits the number of candidate latent factors.
# lambda: higher value makes individual factors sparser.
# spec: higher value retains fewer outcome-associated latent factors.
params$delta <- c(0.01)
params$lambda <- c(0.5, 1.0)
params$spec <- 0.20

# Keep all other established settings explicit for reproducibility.
params$fdr <- 0.10
params$alpha_level <- 0.05
params$thresh_fdr <- 0.20
params$eval_type <- "corr"
params$SLIDE_iter <- 1000
params$SLIDE_top_feats <- 15
params$do_interacts <- FALSE
params$sampleCV_iter <- 500
params$sampleCV_K <- 4
params$y_factor <- FALSE
params$y_levels <- NULL

yaml::write_yaml(params, RETRY_YAML)

set.seed(42)
SLIDE::checkDataParams(params)

message("Starting stable SLIDE optimization at: ", Sys.time())
message("Output directory: ", RESULTS_DIR)
message("Tuning: delta = ", paste(params$delta, collapse = ", "),
        "; lambda = ", paste(params$lambda, collapse = ", "),
        "; spec = ", params$spec)

SLIDE::optimizeSLIDE(params, sink_file = TRUE)

message("SLIDE optimization finished at: ", Sys.time())
