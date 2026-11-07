frequency_analysis <- function(data,study) {
  analytic<-list(SINGLE_COPY_LOCI=study$single,MULTICOPY_STATE_COLUMNS=study$multicopy,
    LOCUS_DISPLAY_NAMES=study$displays[study$single],FORENSIC_MARKER_DISPLAY_NAMES=study$displays)
  e<-recipe_environment("analysis_frequency",analytic,list(),study)
  e$calculate_all_allele_frequency_tables(data)
}
locus_summary <- function(data,study,slots) {
  e<-method_env("analysis_locus")
  rows<-lapply(study$single,function(f){x<-data[[f]];x[x=="DEL"]<-NA_character_;cbind(Locus=study$displays[[f]],e$calc_single_locus_stats(x))})
  for(f in names(slots))rows[[length(rows)+1L]]<-cbind(Locus=study$displays[[f]],e$calc_multicopy_stats(data[[f]],data[slots[[f]]]))
  x<-as.data.frame(do.call(rbind,rows));names(x)[names(x)=="Na_Np"]<-"Na/Np";names(x)[names(x)=="Observed_Range"]<-"Observed range"
  x$Locus[x$Locus %in% study$rapidly_mutating]<-paste0(x$Locus[x$Locus %in% study$rapidly_mutating],"*")
  x
}
