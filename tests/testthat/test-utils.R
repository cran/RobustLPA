test_that("the Hungarian solver finds the optimal assignment", {
  set.seed(1)
  perms <- function(v) {
    if (length(v) <= 1) return(list(v))
    out <- list()
    for (i in seq_along(v)) for (rest in perms(v[-i])) out[[length(out) + 1]] <- c(v[i], rest)
    out
  }
  all_p <- perms(1:5)
  for (rep in 1:20) {
    cost <- matrix(stats::runif(25), 5, 5)
    brute <- min(vapply(all_p, function(pp) sum(cost[cbind(1:5, pp)]), numeric(1)))
    sol <- as.integer(RobustLPA:::solve_lsap_cpp(cost)) + 1L
    expect_equal(sum(cost[cbind(1:5, sol)]), brute)
  }
})

test_that("profile labels are matched by the globally optimal assignment", {
  # A greedy nearest-mean matching would give reference 1 the new profile at
  # 0.6 (total cost 0.36 + 36); the optimal assignment swaps them (25 + 0.16).
  ref <- list(0, 1)
  new <- list(0.6, -5)
  expect_equal(RobustLPA:::.match_profile_labels(ref, new), c(2L, 1L))
  expect_equal(RobustLPA:::.match_profile_labels(list(c(0, 0), c(3, 3)), list(c(3.1, 2.9), c(0.1, 0))), c(2L, 1L))
})

test_that("log-sum-exp is stable for extreme log-densities", {
  M <- rbind(c(-1e4, -1e4 - 1), c(0, -Inf))
  lse <- RobustLPA:::.row_logsumexp(M)
  expect_equal(lse[1], -1e4 + log(1 + exp(-1)))
  expect_equal(lse[2], 0)
})

test_that("PSOCK workers load the same copy of the package as the session", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  pkg_path <- getNamespaceInfo("RobustLPA", "path")
  skip_if_not(file.exists(file.path(pkg_path, "Meta", "package.rds")), "package not installed (load_all)")
  cl <- parallel::makeCluster(1)
  on.exit(parallel::stopCluster(cl))
  expect_silent(RobustLPA:::.setup_psock_workers(cl))
  worker_path <- parallel::clusterCall(cl, function() getNamespaceInfo("RobustLPA", "path"))[[1]]
  expect_identical(normalizePath(worker_path), normalizePath(pkg_path))
})
