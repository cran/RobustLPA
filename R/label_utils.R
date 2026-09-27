#' Match One Set of Profile Labels to a Reference Set
#'
#' Mixture-model profile labels are arbitrary ("label switching"): a fresh
#' fit, a bootstrap refit, or a later MCMC draw can call "profile 1" what
#' the reference calls "profile 2". This finds the permutation of
#' \code{new_means} that minimizes the total squared (optionally
#' standardized) distance to \code{orig_means}, solved \emph{exactly} as a
#' linear assignment problem with the Hungarian algorithm (a greedy
#' nearest-mean matching can mis-assign when two profiles are close to the
#' same reference profile).
#'
#' @param orig_means A list of length \code{G} (the reference labeling's mean vectors).
#' @param new_means A list of length \code{G} (the labeling to align).
#' @param scale Optional numeric vector of length \code{p} used to
#'   standardize each variable before computing distances (e.g. pooled
#'   within-profile standard deviations). Default \code{NULL}: no scaling.
#' @return An integer vector of length \code{G}: element \code{g} is the index
#'   (in \code{new_means}) matched to reference profile \code{g}.
#' @keywords internal
#' @noRd
.match_profile_labels <- function(orig_means, new_means, scale = NULL) {
  G <- length(orig_means)
  if (G <= 1) return(seq_len(G))
  orig_mat <- do.call(rbind, lapply(orig_means, as.numeric))
  new_mat <- do.call(rbind, lapply(new_means, as.numeric))
  if (!is.null(scale)) {
    orig_mat <- sweep(orig_mat, 2, scale, "/")
    new_mat <- sweep(new_mat, 2, scale, "/")
  }
  cost <- matrix(0, G, G)
  for (g in seq_len(G)) {
    cost[g, ] <- rowSums(sweep(new_mat, 2, orig_mat[g, ], "-")^2)
  }
  as.integer(solve_lsap_cpp(cost)) + 1L
}

#' Relabel MCMC Draws to a Common Profile Ordering (Label-Switching Correction)
#'
#' Every draw of every chain is relabeled to a fixed pivot -- the preliminary
#' EM solution used to initialize the chains -- by an exact optimal
#' assignment of its profile means to the pivot means on standardized
#' variables (pivotal reordering; Marin, Mengersen & Robert, 2005). Because
#' the pivot is fixed and does not depend on the order in which draws are
#' processed, the resulting labeling is the same for every chain and every
#' iteration, which is what pooling draws and computing between-chain
#' diagnostics require.
#'
#' @param chains A list of length \code{n_chains}, each element as returned
#'   by \code{mcmc_chain_cpp()} (with \code{mu_chain}, \code{sigma_chain},
#'   \code{pi_chain}).
#' @param G Integer, the number of latent profiles.
#' @param pivot_means A list of \code{G} reference mean vectors.
#' @param scale Numeric vector of length \code{p} (pooled within-profile
#'   standard deviations) used to standardize distances.
#' @return A list of the same shape as \code{chains}, relabeled.
#' @references
#'   Marin, J.-M., Mengersen, K., & Robert, C. P. (2005). Bayesian modelling
#'   and inference on mixtures of distributions. \emph{Handbook of
#'   Statistics}, 25, 459-507. \doi{10.1016/S0169-7161(05)25016-2}
#' @keywords internal
#' @noRd
.relabel_mcmc_chains <- function(chains, G, pivot_means, scale = NULL) {
  if (G <= 1) return(chains)
  identity_perm <- seq_len(G)
  for (chain_id in seq_along(chains)) {
    ch <- chains[[chain_id]]
    n_iter <- length(ch$mu_chain)
    for (iter in seq_len(n_iter)) {
      perm <- .match_profile_labels(pivot_means, ch$mu_chain[[iter]], scale = scale)
      if (!identical(perm, identity_perm)) {
        ch$mu_chain[[iter]] <- ch$mu_chain[[iter]][perm]
        ch$sigma_chain[[iter]] <- ch$sigma_chain[[iter]][perm]
        ch$pi_chain[iter, ] <- ch$pi_chain[iter, perm]
      }
    }
    chains[[chain_id]] <- ch
  }
  chains
}
