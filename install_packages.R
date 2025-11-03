#!/usr/bin/env Rscript

# Install required R packages for the pipeline (CRAN + GitHub)
# - Core pipeline code uses summary_pipeline.R and run_pipeline.R
# - Installs CRAN deps first, then attempts GitHub install for LAVAKnock

message_line <- function(...) cat(paste0(..., "\n"))

failed_cran <- character()
failed_bioc <- character()
failed_python <- character()

ensure_repos <- function() {
  r <- getOption("repos")
  if (is.null(r) || length(r) == 0 || is.na(r[["CRAN"]]) || r[["CRAN"]] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
}

install_if_missing <- function(pkgs, upgrade = FALSE) {
  ensure_repos()
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message_line("Installing ", pkg, " from CRAN...")
      ok <- FALSE
      try({
        install.packages(pkg, dependencies = TRUE, quiet = TRUE)
        ok <- requireNamespace(pkg, quietly = TRUE)
      }, silent = TRUE)
      if (!ok) {
        # Fallback attempts: binary for mac/win, source for linux
        fallback_type <- ifelse(.Platform$OS.type == "windows" || Sys.info()["sysname"] == "Darwin", "binary", "source")
        message_line("  ... retry with type='", fallback_type, "'")
        try({
          install.packages(pkg, dependencies = TRUE, quiet = TRUE, type = fallback_type)
          ok <- requireNamespace(pkg, quietly = TRUE)
        }, silent = TRUE)
      }
      if (!ok) {
        message_line("  ✗ Failed to install ", pkg, " from CRAN")
        failed_cran <<- unique(c(failed_cran, pkg))
      } else {
        message_line("  ✓ Installed ", pkg)
      }
    } else if (isTRUE(upgrade)) {
      message_line("Upgrading ", pkg, " from CRAN...")
      tryCatch({
        install.packages(pkg, dependencies = TRUE, quiet = TRUE)
      }, error = function(e) {
        message_line("  ✗ Failed to upgrade ", pkg, ": ", e$message)
      })
    } else {
      message_line("✓ ", pkg, " already installed")
    }
  }
}

install_bioc_if_missing <- function(pkgs, bioc_version = NULL) {
  ensure_repos()
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message_line("Installing BiocManager from CRAN ...")
    install.packages("BiocManager", dependencies = TRUE, quiet = TRUE)
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message_line("  ✗ Failed to load BiocManager; cannot install Bioconductor packages")
    return(invisible())
  }
  if (!is.null(bioc_version)) {
    tryCatch({
      BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE)
    }, error = function(e) {
      message_line("  ⚠️  Unable to switch Bioconductor version: ", e$message)
    })
  }
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message_line("Installing ", pkg, " from Bioconductor ...")
      ok <- FALSE
      try({
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
        ok <- requireNamespace(pkg, quietly = TRUE)
      }, silent = TRUE)
      if (!ok) {
        message_line("  ✗ Failed to install ", pkg, " via Bioconductor")
        failed_bioc <<- unique(c(failed_bioc, pkg))
      } else {
        message_line("  ✓ Installed ", pkg, " (Bioconductor)")
      }
    } else {
      message_line("✓ ", pkg, " already installed (Bioconductor)")
    }
  }
}

install_lavaknock <- function() {
  if (!requireNamespace("LAVAKnock", quietly = TRUE)) {
    message_line("Installing LAVAKnock from GitHub: shiyangm/LAVA-Knock ...")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install_if_missing("devtools")
    }
    ok <- FALSE
    try({ devtools::install_github("shiyangm/LAVA-Knock", upgrade = "never", dependencies = TRUE, quiet = TRUE); ok <- TRUE }, silent = TRUE)
    if (!ok) {
      message_line("  ✗ Could not install LAVAKnock from GitHub. Please try manually:")
      message_line("    remotes::install_github('shiyangm/LAVA-Knock')")
    }
  } else {
    message_line("✓ LAVAKnock already installed")
  }
}

install_python_packages <- function(packages = c("numpy", "numba", "pandas"), envname = NULL, method = "auto") {
  if (length(packages) == 0) return(invisible())
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    message_line("⚠️  reticulate 未安装，跳过 Python 依赖自动安装")
    return(invisible())
  }
  message_line("Checking Python configuration via reticulate ...")
  cfg <- tryCatch(reticulate::py_config(), error = identity)
  if (inherits(cfg, "error")) {
    message_line("  ⚠️  无法检测现有 Python 解释器: ", cfg$message)
    message_line("     可手动设置 RETICULATE_PYTHON 或提供 envname 以创建虚拟环境")
  }
  message_line("Installing Python packages via reticulate::py_install ...")
  ok <- TRUE
  tryCatch({
    reticulate::py_install(packages, envname = envname, method = method, pip = TRUE)
  }, error = function(e) {
    ok <<- FALSE
    message_line("  ✗ Failed to install Python packages via reticulate: ", e$message)
  })
  if (ok) {
    message_line("✓ Python packages installed/verified: ", paste(packages, collapse = ", "))
  } else {
    failed_python <<- unique(c(failed_python, packages))
    message_line("    建议手动执行: pip install ", paste(packages, collapse = " "))
  }
}

main <- function() {
  message_line("============================================================")
  message_line("Installing packages for GhostKnockoff + LAVA-Knock pipeline")
  message_line("============================================================")

  args <- commandArgs(trailingOnly = TRUE)
  python_path <- NULL
  envname <- NULL
  for (arg in args) {
    if (grepl("^--python=", arg)) python_path <- sub("^--python=", "", arg)
    if (grepl("^--envname=", arg)) envname <- sub("^--envname=", "", arg)
  }

  if (!is.null(python_path) && nzchar(python_path)) {
    Sys.setenv(RETICULATE_PYTHON = python_path)
    message_line("Using RETICULATE_PYTHON = ", python_path)
  } else {
    message_line("RETICULATE_PYTHON not provided; reticulate will auto-detect.")
  }
  if (!is.null(envname) && nzchar(envname)) {
    message_line("reticulate envname: ", envname)
  }

  cran_core <- c(
    # base utilities
    "Matrix", "MASS", "corpcor", "parallel", "doParallel", "foreach",
    # data + tidy
    "data.table", "dplyr", "readr", "stringr",
    # CLI
    "optparse",
    # stats / methods used by pipeline
    "SKAT", "SPAtest", "CompQuadForm", "irlba", "matrixsampling",
    "qqman",
    # knockoff components
    "GhostKnockoff",
    # optional acceleration
    "reticulate",
    # helpers for GitHub installs
    "devtools"
  )

  message_line("\n==> Installing CRAN packages")
  install_if_missing(cran_core)

  bioc_core <- c("snpStats", "graph", "MatrixGenerics")
  message_line("\n==> Installing Bioconductor packages")
  install_bioc_if_missing(bioc_core, bioc_version = "3.17")

  # Attempt GitHub-only dependency
  message_line("\n==> Installing GitHub packages")
  install_lavaknock()

  # Optional Python acceleration dependencies
  message_line("\n==> Installing Python packages via reticulate")
  install_python_packages(packages = c("numpy", "numba", "pandas"), envname = envname)

  if (requireNamespace("reticulate", quietly = TRUE)) {
    message_line("\nPython configuration summary:")
    print(tryCatch(reticulate::py_config(), error = identity))
  }

  # Final check summary
  pkgs_to_check <- c(cran_core, bioc_core, "LAVAKnock")
  missing <- pkgs_to_check[!vapply(pkgs_to_check, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
  if (length(missing) == 0 && length(failed_python) == 0) {
    message_line("\nAll required packages are installed. ✅")
  } else {
    if (length(missing) > 0) {
      message_line("\nSome R packages are still missing:")
      for (m in missing) message_line("  - ", m)
    }
    if (length(failed_python) > 0) {
      message_line("\nPython packages pending: ", paste(unique(failed_python), collapse = ", "))
    }
    if ("Matrix" %in% missing) {
      message_line("\nHint: installing 'Matrix' often requires a compiler toolchain (gfortran/gcc). On macOS you can install Command Line Tools via 'xcode-select --install'; on Linux ensure build-essential and libblas/lapack headers are present.")
    }
    message_line("\nPlease install the above manually, then re-run the pipeline.")
  }
}

if (sys.nframe() == 0) main()
