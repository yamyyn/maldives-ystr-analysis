# Adapted from the user's original analysis_correlogram final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_class_matrix <- function(geographic_matrix, class_label) {
    W <- matrix(0, nrow = nrow(geographic_matrix), ncol = ncol(geographic_matrix), dimnames = dimnames(geographic_matrix))
    if (class_label == "0-<100 km") {
        selected <- geographic_matrix > 0 & geographic_matrix < 100
    }
    else if (class_label == "100-<200 km") {
        selected <- geographic_matrix >= 100 & geographic_matrix < 200
    }
    else if (class_label == "200-<300 km") {
        selected <- geographic_matrix >= 200 & geographic_matrix < 300
    }
    else if (class_label == "300-<400 km") {
        selected <- geographic_matrix >= 300 & geographic_matrix < 400
    }
    else if (class_label == ">=400 km") {
        selected <- geographic_matrix >= 400
    }
    else {
        stop(paste("Unknown class label:", class_label), call. = FALSE)
    }
    W[selected] <- 1
    diag(W) <- 0
    W
}
compute_correlogram_coefficients <- function(genetic_vector, class_vectors) {
    -as.numeric(stats::cor(genetic_vector, class_vectors, method = "pearson"))
}
directional_p <- function(observed, permuted) {
    if (observed >= 0) {
        (1 + sum(permuted >= observed))/(length(permuted) + 1)
    }
    else {
        (1 + sum(permuted <= observed))/(length(permuted) + 1)
    }
}
two_sided_p <- function(observed, permuted) {
    (1 + sum(abs(permuted) >= abs(observed)))/(length(permuted) + 1)
}
progressive_holm <- function(p_values) {
    out <- numeric(length(p_values))
    for (j in seq_along(p_values)) {
        adjusted <- stats::p.adjust(p_values[seq_len(j)], method = "holm")
        out[j] <- adjusted[j]
    }
    out
}
