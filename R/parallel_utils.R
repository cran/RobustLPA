#' Dispatch a Loop of Independent Work Units in Parallel, Reproducibly
#'
#' Shared internal parallel-backend helper used by \code{\link{robust_lpa}}
#' (EM restarts / MCMC chains), \code{\link{estimate_profiles_robust}}
#' (model grid), \code{\link{blrt_robust}} and \code{\link{bch_robust}}
#' (bootstrap replicates).
#'
#' Before dispatching, one integer seed per work unit is drawn from the
#' caller's random-number stream, and each unit starts by calling
#' \code{set.seed()} with its own seed and the caller's RNG kind. As a
#' result, results depend only on the caller's \code{set.seed()} and are
#' identical whether the units run sequentially, forked
#' (\code{parallel::mclapply()}) or on a PSOCK cluster (Windows), and for
#' any number of cores. The caller's RNG state is restored afterwards
#' (advanced only by the seed draw).
#'
#' @param cores Integer, requested number of CPU cores.
#' @param n_units Integer, the number of independent units of work.
#' @param FUN A function of one argument (the unit index, \code{1:n_units}).
#' @param export_vars Character vector of variable names that \code{FUN}
#'   depends on, exported to PSOCK workers.
#' @param export_env The environment in which \code{export_vars} live.
#' @return A list of length \code{n_units}, the result of \code{FUN} applied
#'   to each unit, in order.
#' @keywords internal
#' @noRd
.run_parallel <- function(cores, n_units, FUN, export_vars = character(0), export_env = parent.frame()) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
  unit_seeds <- sample.int(.Machine$integer.max, n_units)
  rng_kind <- RNGkind()
  saved_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", saved_seed, envir = globalenv()), add = TRUE)

  run_unit <- function(k) {
    set.seed(unit_seeds[k], kind = rng_kind[1], normal.kind = rng_kind[2], sample.kind = rng_kind[3])
    FUN(k)
  }

  if (cores > 1 && n_units > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("`cores` > 1 was requested but the 'parallel' package is unavailable; running sequentially.")
    } else {
      n_workers <- min(cores, n_units)
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(n_workers)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        .setup_psock_workers(cl)
        if (length(export_vars) > 0) {
          parallel::clusterExport(cl, varlist = export_vars, envir = export_env)
        }
        return(parallel::parLapply(cl, seq_len(n_units), run_unit))
      }
      res <- parallel::mclapply(seq_len(n_units), run_unit, mc.cores = n_workers)
      failed <- vapply(res, function(r) inherits(r, "try-error"), logical(1))
      if (any(failed)) stop("A parallel worker failed: ", as.character(res[[which(failed)[1]]]))
      return(res)
    }
  }
  lapply(seq_len(n_units), run_unit)
}

#' Load This Session's Copy of RobustLPA on PSOCK Workers
#'
#' PSOCK workers (used on Windows) are fresh R processes: by default they
#' would load RobustLPA from the standard library, which may hold another
#' version than the one loaded in this session (e.g. when a new version is
#' installed in a separate library). The workers therefore put the library of
#' the loaded copy first on their library path before loading the namespace,
#' and a warning is raised if their version still differs.
#' @param cl A cluster from \code{parallel::makeCluster()}.
#' @return \code{NULL}, invisibly.
#' @keywords internal
#' @noRd
.setup_psock_workers <- function(cl) {
  pkg_path <- getNamespaceInfo("RobustLPA", "path")
  installed <- file.exists(file.path(pkg_path, "Meta", "package.rds"))
  lib <- if (installed) dirname(pkg_path) else character(0)
  load_here <- function(lib) {
    if (length(lib) > 0) .libPaths(unique(c(lib, .libPaths())))
    if (!requireNamespace("RobustLPA", quietly = TRUE)) return(NA_character_)
    as.character(getNamespaceVersion("RobustLPA"))
  }
  # A function created here is bound to the RobustLPA namespace, and merely
  # unserializing it would make the worker load RobustLPA from its default
  # library before the library path is changed: detach it from the namespace.
  environment(load_here) <- globalenv()
  worker_versions <- parallel::clusterCall(cl, load_here, lib)
  here <- as.character(getNamespaceVersion("RobustLPA"))
  worker_versions <- unlist(worker_versions)
  if (any(is.na(worker_versions) | worker_versions != here)) {
    warning("Parallel workers could not load RobustLPA ", here, " (they found ",
            paste(unique(worker_versions), collapse = ", "), "). Install the package in a ",
            "library on the default library path, or use `cores = 1`.", call. = FALSE)
  }
  invisible(NULL)
}
