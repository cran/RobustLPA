#' Match One Set of Profile Labels to a Reference Set
#'
#' Mixture-model profile labels are arbitrary ("label switching"): a fresh
#' fit, or even a later iteration/chain of the same MCMC run, can converge
#' with what it calls "profile 1" actually corresponding to a different
#' profile in some reference labeling. This finds, for each reference
#' profile, the closest (by Euclidean distance between mean vectors)
#' not-yet-matched profile in \code{new_means}, via a simple greedy
#' nearest-mean assignment. Used by \code{\link{.relabel_mcmc_chains}} (to
#' align every MCMC iteration/chain to a common labeling before pooling) and
#' by \code{\link{bch_robust}} when \code{correction = "bootstrap"} (to align
#' each bootstrap refit's profiles to the original model's).
#'
#' This is a greedy heuristic, not the globally-optimal (Hungarian-algorithm)
#' assignment; it is adequate for the small numbers of profiles typical of
#' LPA (usually 2-6). It works best when variables are on comparable scales
#' (e.g. standardized data, as recommended throughout this package) and when
#' profiles are reasonably separated -- which is also when label-switching
#' correction matters most.
#'
#' @param orig_means A list of length \code{G} (the reference labeling's mean vectors).
#' @param new_means A list of length \code{G} (the labeling to align).
#' @return An integer vector of length \code{G}: element \code{g} is the index
#'   (in \code{new_means}) matched to reference profile \code{g}.
#' @keywords internal
#' @noRd
.match_profile_labels <- function(orig_means, new_means) {
  G <- length(orig_means)
  if (G <= 1) return(seq_len(G))

  orig_mat <- do.call(rbind, lapply(orig_means, as.numeric))
  new_mat <- do.call(rbind, lapply(new_means, as.numeric))

  assigned <- integer(G)
  used <- logical(G)
  for (g in seq_len(G)) {
    d <- sqrt(rowSums(sweep(new_mat, 2, orig_mat[g, ], "-")^2))
    d[used] <- Inf
    best <- which.min(d)
    assigned[g] <- best
    used[best] <- TRUE
  }
  assigned
}

#' Relabel MCMC Draws to a Common Profile Ordering (Label-Switching Correction)
#'
#' Raw \code{robust_mcmc_cpp()} draws carry no profile identity across
#' iterations or chains: because every profile's prior is exchangeable, the
#' Gibbs sampler can (and in practice does, especially across independently
#' initialized chains) settle on different orderings of the same underlying
#' profiles from one chain -- or even one iteration -- to the next. Naively
#' pooling/averaging raw draws across this "label switching" silently blends
#' together different real-world profiles: the telltale symptom is pooled
#' profile means that are nearly identical to each other and a classification
#' entropy near 0 (every observation looks ~equally likely to belong to every
#' profile), even when the underlying data has well-separated clusters.
#'
#' This relabels every iteration of every chain (in a single deterministic
#' pass: chain 1's iterations in order, then chain 2's, etc.) via a
#' sequential nearest-mean matching (\code{\link{.match_profile_labels}})
#' against a running reference that starts at chain 1's first iterate and is
#' incrementally updated with each newly-aligned iterate. It is a practical
#' approximation to the relabeling algorithms described in the MCMC mixture-
#' model literature (e.g. Stephens, 2000), not a full implementation of any
#' one of them, but directly targets the failure mode above. Called once on
#' the raw output of every chain, before burn-in is even discarded, so that
#' burn-in draws shown by \code{\link{plot_mcmc_chains}} are on the same
#' labeling as the post-burn-in draws used for posterior summaries and
#' \code{\link{.compute_mcmc_diagnostics}}.
#'
#' @param chains A list of length \code{n_chains}, each element as returned
#'   by \code{robust_mcmc_cpp()} (with \code{mu_chain}, \code{sigma_chain}, \code{pi_chain}).
#' @param G Integer, the number of latent profiles.
#' @return A list of the same shape as \code{chains}, relabeled.
#' @keywords internal
#' @noRd
.relabel_mcmc_chains <- function(chains, G) {
  if (G <= 1) return(chains)

  n_chains <- length(chains)
  n_iter <- length(chains[[1]]$mu_chain)
  relabeled <- chains

  ref_mu <- NULL
  ref_n <- 0L

  for (chain_id in seq_len(n_chains)) {
    for (iter in seq_len(n_iter)) {
      mu_iter <- chains[[chain_id]]$mu_chain[[iter]]

      if (is.null(ref_mu)) {
        # Seed the reference labeling with the very first draw processed
        # (chain 1, iteration 1); everything else is aligned to it. There is
        # no "true" absolute labeling in an exchangeable mixture -- only
        # internal consistency across all the stored draws is required.
        ref_mu <- mu_iter
        ref_n <- 1L
        next
      }

      perm <- .match_profile_labels(ref_mu, mu_iter)
      if (!identical(perm, seq_len(G))) {
        relabeled[[chain_id]]$mu_chain[[iter]] <- mu_iter[perm]
        relabeled[[chain_id]]$sigma_chain[[iter]] <- chains[[chain_id]]$sigma_chain[[iter]][perm]
        relabeled[[chain_id]]$pi_chain[iter, ] <- chains[[chain_id]]$pi_chain[iter, perm]
      }

      ref_n <- ref_n + 1L
      aligned_mu <- mu_iter[perm]
      for (g in seq_len(G)) {
        ref_mu[[g]] <- ref_mu[[g]] + (aligned_mu[[g]] - ref_mu[[g]]) / ref_n
      }
    }
  }

  relabeled
}
