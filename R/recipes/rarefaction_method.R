# Original mathematical/parsing helpers; top-level study scripts were not executed.
expected_rarefied_distinct <- function(counts, m) {
    counts <- as.integer(counts)
    N <- sum(counts)
    if (m > N) 
        stop("Rarefaction size exceeds population sample size.", call. = FALSE)
    if (m == N) 
        return(as.numeric(length(counts)))
    absent_prob <- vapply(counts, function(f) {
        if ((N - f) < m) 
            return(0)
        exp(lchoose(N - f, m) - lchoose(N, m))
    }, FUN.VALUE = numeric(1))
    sum(1 - absent_prob)
}
expected_rarefied_singletons <- function(counts, m) {
    counts <- as.integer(counts)
    N <- sum(counts)
    if (m > N) 
        stop("Rarefaction size exceeds population sample size.", call. = FALSE)
    if (m == N) {
        return(as.numeric(sum(counts == 1L)))
    }
    probs <- vapply(counts, function(f) {
        if (m < 1L || (N - f) < (m - 1L)) 
            return(0)
        logp <- log(f) + lchoose(N - f, m - 1L) - lchoose(N, m)
        exp(logp)
    }, FUN.VALUE = numeric(1))
    sum(probs)
}
subsample_one_group <- function(h, m, B, seed_offset = 0L) {
    h <- as.character(h)
    N <- length(h)
    if (m > N) 
        stop("m > N in subsampling.", call. = FALSE)
    observed <- calc_haplotype_stats(h)
    counts <- table(h)
    expected_distinct <- expected_rarefied_distinct(counts, m)
    expected_singletons <- expected_rarefied_singletons(counts, m)
    if (m == N) {
        iter <- tibble(Iteration = seq_len(B), Distinct_Haplotypes = rep(observed$Distinct_Haplotypes, B), Singleton_Haplotypes = rep(observed$Singleton_Haplotypes, 
            B))
    }
    else {
        set.seed(SEED + as.integer(seed_offset))
        iter_list <- vector("list", B)
        for (b in seq_len(B)) {
            s <- sample.int(N, size = m, replace = FALSE)
            tab <- table(h[s])
            iter_list[[b]] <- c(Distinct_Haplotypes = length(tab), Singleton_Haplotypes = sum(tab == 1L))
        }
        mat <- do.call(rbind, iter_list)
        iter <- tibble(Iteration = seq_len(B), Distinct_Haplotypes = as.integer(mat[, "Distinct_Haplotypes"]), Singleton_Haplotypes = as.integer(mat[, 
            "Singleton_Haplotypes"]))
    }
    summary <- tibble(N_Observed = N, Rarefaction_N = m, Observed_Distinct = observed$Distinct_Haplotypes, Observed_Singletons = observed$Singleton_Haplotypes, 
        Expected_Rarefied_Distinct = expected_distinct, Distinct_LCL_2.5 = as.numeric(quantile(iter$Distinct_Haplotypes, 
            0.025, names = FALSE)), Distinct_UCL_97.5 = as.numeric(quantile(iter$Distinct_Haplotypes, 0.975, names = FALSE)), 
        Expected_Rarefied_Singletons = expected_singletons, Singleton_LCL_2.5 = as.numeric(quantile(iter$Singleton_Haplotypes, 
            0.025, names = FALSE)), Singleton_UCL_97.5 = as.numeric(quantile(iter$Singleton_Haplotypes, 0.975, names = FALSE)), 
        MC_Mean_Distinct = mean(iter$Distinct_Haplotypes), MC_Mean_Singletons = mean(iter$Singleton_Haplotypes))
    list(summary = summary, iterations = iter)
}
