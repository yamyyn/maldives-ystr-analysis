# Haploid method-of-moments AMOVA for squared repeat-size distances.
# No scaling by locus variance; negative variance components are retained.
amova_point <- function(X, population, region = NULL) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  assert(all(is.finite(X)) && all(X == round(X)) && nrow(X) == length(population),
         "AMOVA requires aligned complete integer repeat calls")
  pops <- unique(population); N <- nrow(X); K <- length(pops)
  assert(K >= 2 && N > K, "AMOVA needs >=2 populations and within-population replication")
  indices <- split(seq_len(N), factor(population, levels = pops))
  ni <- unname(lengths(indices)); grand <- colMeans(X)
  means <- lapply(indices, function(i) colMeans(X[i, , drop = FALSE]))
  within <- sum(vapply(seq_along(indices), function(j)
    sum(sweep(X[indices[[j]], , drop = FALSE], 2, means[[j]], "-")^2), 0.0))
  msw <- within / (N-K)
  if (is.null(region)) {
    between <- sum(vapply(seq_along(pops), function(j) ni[j] * sum((means[[j]]-grand)^2), 0.0))
    coef <- (N-sum(ni^2)/N)/(K-1)
    vb <- (between/(K-1)-msw)/coef; total <- vb+msw
    return(c(N = N, populations = K, regions = NA_real_, loci = ncol(X),
      SS_region = NA_real_, SS_population = between, SS_within = within,
      Va = NA_real_, Vb = vb, Vc = msw, RCT = NA_real_, RSC = NA_real_,
      RST = if (abs(total) > 1e-14) vb/total else NA_real_,
      percent_region = NA_real_, percent_population = if (abs(total)>1e-14) 100*vb/total else NA_real_,
      percent_within = if (abs(total)>1e-14) 100*msw/total else NA_real_))
  }
  assert(length(region) == N && !anyNA(region), "Invalid AMOVA region labels")
  reg_of_pop <- vapply(indices, function(i) {
    assert(length(unique(region[i])) == 1L, "Populations must nest within regions")
    region[i[1]]
  }, "")
  regs <- unique(region); G <- length(regs)
  assert(G >= 2 && K > G, "Nested AMOVA needs >=2 regions and more populations than regions")
  ng <- table(factor(region, levels = regs))
  sa <- sb <- 0
  for (g in regs) {
    gi <- which(region == g); gm <- colMeans(X[gi, , drop = FALSE])
    sa <- sa + length(gi)*sum((gm-grand)^2)
    for (j in which(reg_of_pop == g)) sb <- sb + ni[j]*sum((means[[j]]-gm)^2)
  }
  ng_pop <- as.numeric(ng[match(reg_of_pop, regs)])
  A <- (N-sum(ng^2)/N)/(G-1)
  B <- sum(ni^2*(1/ng_pop-1/N))/(G-1)
  C <- sum(ni-ni^2/ng_pop)/(K-G)
  vb <- (sb/(K-G)-msw)/C
  va <- (sa/(G-1)-msw-B*vb)/A
  total <- va+vb+msw
  ratio <- function(num, den) if (abs(den) > 1e-14) num/den else NA_real_
  c(N = N, populations = K, regions = G, loci = ncol(X),
    SS_region = sa, SS_population = sb, SS_within = within,
    Va = va, Vb = vb, Vc = msw, RCT = ratio(va,total), RSC = ratio(vb,vb+msw),
    RST = ratio(va+vb,total), percent_region = 100*ratio(va,total),
    percent_population = 100*ratio(vb,total), percent_within = 100*ratio(msw,total))
}

amova_permutations <- function(X, pop, reg, B, seed) {
  observed <- amova_point(X, pop, reg)[c("RCT", "RSC", "RST")]
  set.seed(seed); draws <- matrix(NA_real_, B, 3, dimnames = list(NULL, names(observed)))
  pop_names <- unique(pop); pop_regions <- vapply(pop_names, function(p) reg[match(p, pop)], "")
  for (b in seq_len(B)) {
    # Among-region test: populations are whole exchangeable units.
    reassigned <- pop_regions[sample.int(length(pop_names))]
    rg <- reassigned[match(pop, pop_names)]
    draws[b, "RCT"] <- amova_point(X, pop, rg)["RCT"]
    # Among-populations-within-region: permute complete haplotypes within regions.
    order <- seq_len(nrow(X))
    for (g in unique(reg)) {
      i <- which(reg == g); order[i] <- i[sample.int(length(i))]
    }
    draws[b, "RSC"] <- amova_point(X[order, , drop = FALSE], pop, reg)["RSC"]
    draws[b, "RST"] <- amova_point(X[sample.int(nrow(X)), , drop = FALSE], pop, reg)["RST"]
  }
  valid <- colSums(is.finite(draws))
  p <- vapply(seq_along(observed), function(j) {
    if (!is.finite(observed[j]) || valid[j] != B) return(NA_real_)
    (1+sum(draws[, j] >= observed[j]-1e-14))/(B+1)
  }, 0.0)
  data.frame(statistic = names(observed), estimate = unname(observed), p_greater = p,
    permutations = B, valid_permutations = valid, seed = seed)
}

amova_sensitivity <- function(d, populations, panels, cfg) {
  rows <- list()
  for (panel_name in names(panels)) {
    cohort <- structure_cohort(d, panels[[panel_name]], cfg)
    z <- cohort$data
    reg <- populations$region[match(z$population, populations$population)]
    for (drop in c("NONE", sort(unique(z$population)))) {
      keep <- if (drop == "NONE") rep(TRUE, nrow(z)) else z$population != drop
      r <- tryCatch(amova_point(z[keep, cohort$loci, drop = FALSE], z$population[keep], reg[keep]),
                    error = function(e) e)
      if (inherits(r, "error")) rows[[length(rows)+1L]] <- data.frame(panel = panel_name,
        excluded_population = drop, status = conditionMessage(r), statistic = "unavailable", value = NA_real_)
      else rows[[length(rows)+1L]] <- data.frame(panel = panel_name, excluded_population = drop,
        status = "estimated", statistic = names(r), value = unname(r))
    }
  }
  do.call(rbind, rows)
}

pairwise_rst_verification <- function(X, groups, B, seed) {
  set.seed(seed); pairs <- combn(sort(unique(groups)), 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    keep <- groups %in% pair; x <- X[keep, , drop = FALSE]; g <- groups[keep]
    observed <- amova_point(x, g)["RST"]
    draw <- replicate(B, amova_point(x[sample.int(nrow(x)), , drop = FALSE], g)["RST"])
    p <- if (all(is.finite(draw)) && is.finite(observed)) (1+sum(draw >= observed-1e-14))/(B+1) else NA_real_
    data.frame(population1 = pair[1], population2 = pair[2], RST = observed,
      p_greater = p, monte_carlo_se = if (is.finite(p)) sqrt(p*(1-p)/(B+1)) else NA_real_,
      p_resolution = 1/(B+1), permutations = B)
  })
  result <- do.call(rbind, rows)
  result$p_bonferroni <- p.adjust(result$p_greater, "bonferroni")
  result$p_BH <- p.adjust(result$p_greater, "BH")
  result
}
