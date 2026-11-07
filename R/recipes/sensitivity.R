# Original mathematical/parsing helpers; top-level study scripts were not executed.
compute_metrics_from_counts <- function(counts) {
    counts <- counts[counts > 0]
    n <- sum(counts)
    if (n == 0) {
        stop("Cannot calculate forensic metrics for an empty scope.", call. = FALSE)
    }
    distinct_haplotypes <- length(counts)
    unique_haplotypes <- sum(counts == 1)
    nonunique_distinct <- sum(counts > 1)
    profiles_repeated <- sum(counts[counts > 1])
    frequencies <- counts/n
    hmp <- sum(frequencies^2)
    hd <- if (n > 1) {
        n/(n - 1) * (1 - hmp)
    }
    else {
        NA_real_
    }
    dc <- distinct_haplotypes/n
    unique_profile_fraction <- unique_haplotypes/n
    unique_haplotype_fraction <- unique_haplotypes/distinct_haplotypes
    max_count <- max(counts)
    max_frequency <- max_count/n
    c(Distinct_Haplotypes = distinct_haplotypes, Unique_Haplotypes = unique_haplotypes, Nonunique_Distinct_Haplotypes = nonunique_distinct, 
        Profiles_in_Repeated_Haplotypes = profiles_repeated, HMP = hmp, HD = hd, DC = dc, Unique_Profile_Fraction = unique_profile_fraction, 
        Unique_Haplotype_Fraction = unique_haplotype_fraction, Most_Common_Haplotype_Count = max_count, Most_Common_Haplotype_Frequency = max_frequency)
}
compute_metrics_from_haps <- function(haps) {
    compute_metrics_from_counts(unname(table(haps)))
}
sample_vector_safe <- function(x, size, replace = FALSE) {
    if (length(x) == 0) {
        if (size == 0) 
            return(x)
        stop("Cannot sample a positive number of values from an empty vector.", call. = FALSE)
    }
    if (!replace && size > length(x)) {
        stop("Requested sample size exceeds vector length without replacement.", call. = FALSE)
    }
    x[sample.int(n = length(x), size = size, replace = replace)]
}
