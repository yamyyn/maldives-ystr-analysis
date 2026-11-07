# Original mathematical/parsing helpers; top-level study scripts were not executed.
pairwise_sharing <- function(data, group_var, group_order) {
    group_var <- rlang::ensym(group_var)
    sets <- data %>% mutate(.group = as.character(!!group_var)) %>% distinct(.group, Haplotype_ID) %>% group_by(.group) %>% 
        summarise(Haplotypes = list(Haplotype_ID), .groups = "drop")
    counts <- data %>% mutate(.group = as.character(!!group_var)) %>% count(.group, Haplotype_ID, name = "n")
    shared_distinct <- matrix(0, nrow = length(group_order), ncol = length(group_order), dimnames = list(group_order, group_order))
    matching_pairs <- shared_distinct
    jaccard <- matrix(NA_real_, nrow = length(group_order), ncol = length(group_order), dimnames = list(group_order, group_order))
    for (i in seq_along(group_order)) {
        g1 <- group_order[i]
        h1 <- sets$Haplotypes[match(g1, sets$.group)][[1]]
        for (j in seq_along(group_order)) {
            g2 <- group_order[j]
            h2 <- sets$Haplotypes[match(g2, sets$.group)][[1]]
            shared <- intersect(h1, h2)
            union_h <- union(h1, h2)
            shared_distinct[i, j] <- length(shared)
            jaccard[i, j] <- if (length(union_h) > 0L) {
                length(shared)/length(union_h)
            }
            else {
                NA_real_
            }
            if (i == j) {
                cg <- counts %>% filter(.group == g1)
                matching_pairs[i, j] <- sum(choose(cg$n, 2))
            }
            else {
                c1 <- counts %>% filter(.group == g1) %>% select(Haplotype_ID, n1 = n)
                c2 <- counts %>% filter(.group == g2) %>% select(Haplotype_ID, n2 = n)
                cross <- inner_join(c1, c2, by = "Haplotype_ID")
                matching_pairs[i, j] <- if (nrow(cross) == 0L) {
                  0
                }
                else {
                  sum(cross$n1 * cross$n2)
                }
            }
        }
    }
    list(shared_distinct = shared_distinct, matching_pairs = matching_pairs, jaccard = jaccard)
}
