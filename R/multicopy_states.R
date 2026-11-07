# Shared observed-state normalizer. Never infers unobserved copies.
observed_multicopy_state <- function(values, collapse_identical=FALSE) {
 values<-as.character(values);values<-values[!is.na(values)&nzchar(trimws(values))]
 if(!length(values))return(NA_character_)
 numbers<-suppressWarnings(as.numeric(values))
 assert(all(is.finite(numbers)),"Multicopy states require finite observed numeric alleles")
 values<-values[order(numbers,values)];numbers<-sort(numbers)
 if(collapse_identical)values<-values[!duplicated(numbers)]
 paste(values,collapse="-")
}
