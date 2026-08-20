#' Dispatch a Loop of Independent Work Units in Parallel
#'
#' Shared internal parallel-backend helper used by both \code{\link{robust_lpa}}
#' (to parallelize EM random restarts / MCMC chains via its \code{cores}
#' argument) and \code{\link{blrt_robust}} (to parallelize its bootstrap
#' replicates via its own \code{cores} argument). Centralizing this logic in
#' one place keeps the parallel dispatch behavior (backend selection, RNG
#' handling, fallback conditions) identical and independently correct in only
#' one location.
#'
#' Falls back to sequential \code{lapply()} if \code{cores == 1}, the
#' \pkg{parallel} package is unavailable, or there is only one unit of work to
#' run anyway (parallelizing a single task is pure overhead). On Windows
#' (no \code{fork()}), uses a \code{parallel::makeCluster()} PSOCK cluster with
#' independent per-worker RNG streams (\code{parallel::clusterSetRNGStream()});
#' elsewhere uses \code{parallel::mclapply()} forking, whose per-worker RNG
#' independence is automatic (\code{mc.set.seed = TRUE} is the default).
#'
#' @param cores Integer, requested number of CPU cores.
#' @param n_units Integer, the number of independent units of work (e.g.
#'   EM restarts, MCMC chains, or bootstrap replicates).
#' @param FUN A function of one argument (the unit index, \code{1:n_units})
#'   that performs one unit of work and returns its result. Must be
#'   self-contained (a closure capturing everything it needs from its
#'   defining environment), so it can be dispatched unchanged to a worker
#'   process.
#' @param export_vars Character vector of variable names that \code{FUN}
#'   depends on, to be exported to PSOCK cluster workers via
#'   \code{parallel::clusterExport()} (ignored by the \code{mclapply()} and
#'   sequential code paths, where \code{FUN} already has direct access to
#'   these variables through its closure).
#' @param export_env The environment in which \code{export_vars} live
#'   (typically \code{environment()} of the caller).
#' @return A list of length \code{n_units}, the result of \code{FUN} applied
#'   to each unit, in order.
#' @keywords internal
#' @noRd
.run_parallel <- function(cores, n_units, FUN, export_vars, export_env) {
  if (cores > 1 && n_units > 1 && requireNamespace("parallel", quietly = TRUE)) {
    n_workers <- min(cores, n_units)
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, {
        if (requireNamespace("RobustLPA", quietly = TRUE)) library(RobustLPA)
      })
      parallel::clusterSetRNGStream(cl)
      parallel::clusterExport(cl, varlist = export_vars, envir = export_env)
      return(parallel::parLapply(cl, seq_len(n_units), FUN))
    }
    return(parallel::mclapply(seq_len(n_units), FUN, mc.cores = n_workers))
  }
  if (cores > 1 && n_units > 1 && !requireNamespace("parallel", quietly = TRUE)) {
    warning("`cores` > 1 was requested but the 'parallel' package is unavailable; running sequentially.")
  }
  lapply(seq_len(n_units), FUN)
}
