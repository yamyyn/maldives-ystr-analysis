# Release adaptation: shared observed-state normalization, no inferred copies.
make_multicopy_state <- function(..., collapse_identical=FALSE) {
 observed_multicopy_state(c(...),collapse_identical=collapse_identical)
}
