# Additional descriptive verification. Run after the full pipeline.
source("R/pipeline.R")
cfg<-load_pipeline_config(Sys.getenv("YSTR_CONFIG","config/analysis_settings.yml"))
checkpoint<-readRDS(file.path(cfg$out,"pipeline_state.rds"))
assert(identical(checkpoint$input_hash,pipeline_input_hash(cfg)),"Inputs or code changed; start a new run")
cohort<-checkpoint$state$`02_build_datasets`
assert(!is.null(cohort),"Dataset preparation has not completed")
markers<-c(cohort$study$single,names(cohort$study$multicopy))
full<-cohort$full;complete<-cohort$complete
excluded<-full[!full$sample_id %in% complete$sample_id,,drop=FALSE]
re<-method_env("rarefaction_method")
groups<-sort(unique(full$atoll))
complete_sizes<-table(factor(complete$atoll,levels=groups))
assert(all(complete_sizes>0L),"A division has no complete profiles; define a different comparison scope explicitly")
common_n<-min(complete_sizes)
summarize<-function(data,label){
 h<-haplotype_keys(data,markers)
 do.call(rbind,lapply(groups,function(g){
  ii<-which(data$atoll==g);counts<-table(h[ii])
  data.frame(Dataset=label,Division=g,as.list(hap_metrics(h[ii])),
   Rarefaction_n=common_n,Expected_distinct=re$expected_rarefied_distinct(counts,common_n),
   Expected_singletons=re$expected_rarefied_singletons(counts,common_n),check.names=FALSE)
 }))
}
dest<-file.path(cfg$out,"verification");dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write_csv(rbind(summarize(full,"Full"),summarize(complete,"Complete")),file.path(dest,"Division_Common_Size_Rarefaction.csv"))
if(nrow(excluded)){
 keys<-haplotype_keys(excluded,markers)
 out<-do.call(rbind,lapply(split(seq_len(nrow(excluded)),excluded$atoll),function(ii){
  data.frame(Division=excluded$atoll[ii[1]],Profiles=length(ii),Distinct_haplotypes=length(unique(keys[ii])),Singleton_haplotypes=sum(table(keys[ii])==1L))
 }))
 write_csv(out,file.path(dest,"Excluded_Profile_Haplotypes_By_Division.csv"))
 write_csv(data.frame(Profiles=length(keys),Distinct_haplotypes=length(unique(keys)),Singleton_haplotypes=sum(table(keys)==1L)),file.path(dest,"Excluded_Profile_Haplotypes_Overall.csv"))
}
writeLines(c("Analytical hypergeometric rarefaction at the smallest complete-division sample size.","Descriptive comparison only; no new randomization test or inferred copy-number/deletion claim."),file.path(dest,"Local_Profile_Sensitivity_Methods.txt"))
message("Wrote local profile sensitivity to ",dest)
