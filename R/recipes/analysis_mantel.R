# Adapted from the user's original analysis_mantel final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
lower_vector <- function(mat) {
    mat[lower.tri(mat)]
}
mantel_observed <- function(genetic_matrix, geographic_matrix) {
    stats::cor(lower_vector(genetic_matrix), lower_vector(geographic_matrix), method = "pearson")
}
mantel_permute <- function(genetic_matrix, geographic_matrix, permutations, seed) {
    set.seed(seed)
    n <- nrow(genetic_matrix)
    geo_vector <- lower_vector(geographic_matrix)
    output <- numeric(permutations)
    for (b in seq_len(permutations)) {
        perm <- sample.int(n, size = n, replace = FALSE)
        permuted_genetic <- genetic_matrix[perm, perm, drop = FALSE]
        output[b] <- stats::cor(lower_vector(permuted_genetic), geo_vector, method = "pearson")
    }
    output
}
empirical_p_greater <- function(observed, permutations) {
    (1 + sum(permutations >= observed))/(length(permutations) + 1)
}
empirical_p_two_sided <- function(observed, permutations) {
    (1 + sum(abs(permutations) >= abs(observed)))/(length(permutations) + 1)
}
