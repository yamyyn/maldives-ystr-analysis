# Original mathematical/parsing helpers; top-level study scripts were not executed.
gower_eigenvalues <- function(distance_matrix) {
    n <- nrow(distance_matrix)
    J <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    B <- -0.5 * J %*% (distance_matrix^2) %*% J
    eigen(B, symmetric = TRUE, only.values = TRUE)$values
}
orient_axis <- function(x) {
    idx <- which.max(abs(x))
    if (length(idx) == 0 || x[idx] >= 0) {
        return(x)
    }
    -x
}
edge_name <- function(node_number, tree) {
    n_tip <- ape::Ntip(tree)
    if (node_number <= n_tip) {
        tree$tip.label[node_number]
    }
    else {
        paste0("InternalNode_", node_number - n_tip)
    }
}
