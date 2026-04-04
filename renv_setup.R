# ============================================================
# renv_setup.R
#
# Run this script ONCE on a fresh clone of the repository to
# initialize the reproducible R environment (B1.2).
#
# Usage (from the project root):
#   Rscript renv_setup.R
#
# After this script completes, all subsequent `quarto render`
# calls will use the pinned package versions recorded in
# renv.lock. The renv.lock file must be committed to version
# control.
# ============================================================

# 1. Install renv if not present
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# 2. Initialize the project environment.
#    consent = TRUE suppresses the interactive prompt.
renv::init(
  bioconductor = FALSE,
  bare         = FALSE,
  restart      = FALSE
)

# 3. Install all packages used by the project.
#    Exact versions will be pinned in renv.lock after snapshot().
pkgs <- c(
  "knitr",      # table rendering
  "rmarkdown",  # document rendering
  "quarto",     # Quarto R interface
  "zoo"         # rolling statistics (optional but used in heteroscedasticity figure)
)

renv::install(pkgs)

# 4. Snapshot the installed versions to renv.lock.
#    This is the file that guarantees reproducibility.
renv::snapshot(prompt = FALSE)

message(
  "\n",
  "renv setup complete.\n",
  "Commit renv.lock and the renv/ directory to version control.\n",
  "To restore on a new machine: renv::restore()\n"
)
