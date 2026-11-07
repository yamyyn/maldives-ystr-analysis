# Original mathematical/parsing helpers; top-level study scripts were not executed.
evaluate_one_permutation <- function(haplotypes, sample_size_grid, iteration) {
    n <- length(haplotypes)
    random_order <- sample.int(n = n, size = n, replace = FALSE)
    ordered_haps <- haplotypes[random_order]
    distinct_values <- integer(length(sample_size_grid))
    singleton_values <- integer(length(sample_size_grid))
    hap_counts <- new.env(hash = TRUE, parent = emptyenv())
    distinct_n <- 0L
    singleton_n <- 0L
    grid_position <- 1L
    for (current_n in seq_len(n)) {
        h <- ordered_haps[current_n]
        if (!exists(h, envir = hap_counts, inherits = FALSE)) {
            assign(h, 1L, envir = hap_counts)
            distinct_n <- distinct_n + 1L
            singleton_n <- singleton_n + 1L
        }
        else {
            previous_count <- get(h, envir = hap_counts, inherits = FALSE)
            new_count <- previous_count + 1L
            assign(h, new_count, envir = hap_counts)
            if (previous_count == 1L) {
                singleton_n <- singleton_n - 1L
            }
        }
        if (grid_position <= length(sample_size_grid) && current_n == sample_size_grid[grid_position]) {
            distinct_values[grid_position] <- distinct_n
            singleton_values[grid_position] <- singleton_n
            grid_position <- grid_position + 1L
        }
    }
    if (grid_position != length(sample_size_grid) + 1L) {
        stop("Not all requested sample sizes were recorded.", call. = FALSE)
    }
    tibble(Iteration = iteration, Sample_Size = sample_size_grid, Distinct_Haplotypes = distinct_values, Singleton_Haplotypes = singleton_values)
}
