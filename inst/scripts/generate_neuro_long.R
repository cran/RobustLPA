# Generates the `neuro_long` dataset shipped with RobustLPA (data/neuro_long.rda).
# Run it with
#   source(system.file("scripts", "generate_neuro_long.R", package = "RobustLPA"))
# and compare with data(neuro_long): the result is identical.
# Synthetic longitudinal neuropsychological data: three latent classes of
# cognitive change, three outcomes on a T-score-like metric, random
# intercepts and slopes, missed visits, drop-out depending on the last
# observed Memory score (missing at random), and gross data-entry errors
# for a few persons.

set.seed(20260926)
N <- 400
years <- 0:5
cls <- sample(1:3, N, replace = TRUE, prob = c(0.55, 0.30, 0.15))
labels <- c("Stable", "Slow decline", "Fast decline")

# class mean trajectories (intercept, slope per year)
mu <- list(
  Memory    = rbind(c(53, 0.0), c(49, -2.0), c(45, -5.0)),
  Executive = rbind(c(52, 0.0), c(48, -1.5), c(44, -4.0)),
  Speed     = rbind(c(50, -0.5), c(50, -0.5), c(50, -0.5))   # no class differences
)
outcomes <- names(mu)

# random effects: (intercept, slope) for each outcome, correlated
sd_re <- rep(c(5, 0.4), 3)
R <- diag(6)
R[cbind(c(1, 1, 3), c(3, 5, 5))] <- 0.5       # intercepts across outcomes
R[cbind(c(2, 2, 4), c(4, 6, 6))] <- 0.4       # slopes across outcomes
R[cbind(c(1, 3, 5), c(2, 4, 6))] <- 0.2       # intercept-slope within outcome
R[lower.tri(R)] <- t(R)[lower.tri(R)]
Sigma_re <- diag(sd_re) %*% R %*% diag(sd_re)
L <- chol(Sigma_re)

rows <- list()
for (i in seq_len(N)) {
  b <- as.vector(stats::rnorm(6) %*% L)
  last_memory <- NA
  for (t in years) {
    if (t > 0) {
      p_drop <- stats::plogis(-3.2 + 0.12 * (50 - last_memory))
      if (stats::runif(1) < p_drop) break
    }
    y <- vapply(seq_along(outcomes), function(k) {
      m <- mu[[k]][cls[i], ]
      m[1] + b[2 * k - 1] + (m[2] + b[2 * k]) * t + stats::rnorm(1, 0, 2.5)
    }, numeric(1))
    if (t > 0 && stats::runif(1) < 0.08) y[] <- NA          # missed visit
    if (stats::runif(1) < 0.05) y[sample(3, 1)] <- NA          # single missing test
    if (!is.na(y[1])) last_memory <- y[1]
    if (is.na(last_memory)) last_memory <- 50
    rows[[length(rows) + 1]] <- data.frame(ID = i, Year = t, Memory = y[1], Executive = y[2], Speed = y[3])
  }
}
neuro_long <- do.call(rbind, rows)

# gross data-entry errors for 4% of the persons
bad <- sample(unique(neuro_long$ID), round(0.04 * N))
for (i in bad) {
  r <- sample(which(neuro_long$ID == i), 1)
  k <- sample(outcomes, 1)
  if (!is.na(neuro_long[r, k])) neuro_long[r, k] <- neuro_long[r, k] + sample(c(-1, 1), 1) * stats::runif(1, 25, 40)
}

baseline <- data.frame(
  ID = seq_len(N),
  True_Class = factor(labels[cls], levels = labels),
  Age = round(stats::rnorm(N, 64 + 4 * (cls == 3) + 1.5 * (cls == 2), 7), 1),
  Biomarker = round(stats::rnorm(N, 900 - 60 * (cls == 2) - 180 * (cls == 3), 150))
)
neuro_long <- merge(neuro_long, baseline, by = "ID", sort = TRUE)
neuro_long <- neuro_long[order(neuro_long$ID, neuro_long$Year),
                         c("ID", "Year", "Memory", "Executive", "Speed", "True_Class", "Age", "Biomarker")]
rownames(neuro_long) <- NULL
neuro_long[c("Memory", "Executive", "Speed")] <- round(neuro_long[c("Memory", "Executive", "Speed")], 1)

# To rebuild data/neuro_long.rda from the package sources:
# save(neuro_long, file = "data/neuro_long.rda", compress = "xz")
