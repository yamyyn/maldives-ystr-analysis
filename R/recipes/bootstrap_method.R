# Original mathematical/parsing helpers; top-level study scripts were not executed.
calc_haplotype_stats <- function(h) {
    h <- as.character(h)
    h <- h[!is.na(h) & h != ""]
    n <- length(h)
    if (n == 0L) {
        return(tibble(N = 0L, Distinct_Haplotypes = 0L, Singleton_Haplotypes = 0L, HMP = NA_real_, HD = NA_real_, DC = NA_real_))
    }
    tab <- table(h)
    p <- as.numeric(tab)/n
    hmp <- sum(p^2)
    tibble(N = as.integer(n), Distinct_Haplotypes = as.integer(length(tab)), Singleton_Haplotypes = as.integer(sum(tab == 
        1L)), HMP = hmp, HD = if (n > 1L) 
        (n/(n - 1L)) * (1 - hmp)
    else NA_real_, DC = length(tab)/n)
}
bootstrap_group <- function(h, B, seed_value) {
    h <- as.character(h)
    n <- length(h)
    lev <- unique(h)
    code <- match(h, lev)
    k <- length(lev)
    observed <- calc_haplotype_stats(h)
    set.seed(seed_value)
    distinct <- integer(B)
    singleton <- integer(B)
    hmp <- numeric(B)
    hd <- numeric(B)
    dc <- numeric(B)
    for (b in seq_len(B)) {
        sampled_codes <- sample(code, size = n, replace = TRUE)
        counts <- tabulate(sampled_codes, nbins = k)
        counts <- counts[counts > 0L]
        p <- counts/n
        hmp_b <- sum(p^2)
        distinct_b <- length(counts)
        distinct[b] <- distinct_b
        singleton[b] <- sum(counts == 1L)
        hmp[b] <- hmp_b
        hd[b] <- if (n > 1L) 
            (n/(n - 1L)) * (1 - hmp_b)
        else NA_real_
        dc[b] <- distinct_b/n
    }
    iterations <- tibble(Iteration = seq_len(B), Distinct_Haplotypes = distinct, Singleton_Haplotypes = singleton, HMP = hmp, 
        HD = hd, DC = dc)
    summarise_metric <- function(x, observed_value) {
        c(Observed = observed_value, Bootstrap_Mean = mean(x, na.rm = TRUE), Bootstrap_SE = sd(x, na.rm = TRUE), LCL_2.5 = as.numeric(quantile(x, 
            0.025, na.rm = TRUE, names = FALSE)), UCL_97.5 = as.numeric(quantile(x, 0.975, na.rm = TRUE, names = FALSE)))
    }
    metric_rows <- bind_rows(as_tibble_row(c(Metric = "Distinct_Haplotypes", summarise_metric(iterations$Distinct_Haplotypes, 
        observed$Distinct_Haplotypes))), as_tibble_row(c(Metric = "Singleton_Haplotypes", summarise_metric(iterations$Singleton_Haplotypes, 
        observed$Singleton_Haplotypes))), as_tibble_row(c(Metric = "HMP", summarise_metric(iterations$HMP, observed$HMP))), 
        as_tibble_row(c(Metric = "HD", summarise_metric(iterations$HD, observed$HD))), as_tibble_row(c(Metric = "DC", summarise_metric(iterations$DC, 
            observed$DC)))) %>% mutate(across(c(Observed, Bootstrap_Mean, Bootstrap_SE, LCL_2.5, UCL_97.5), as.numeric))
    list(observed = observed, summary = metric_rows, iterations = iterations)
}
