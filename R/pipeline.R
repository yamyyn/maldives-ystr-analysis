source("R/runtime.R")
for(f in c("qc","multicopy_states","initial_qc","data_validation","original_audit_engine","forensic_statistics","allele_frequencies","population_structure",
  "amova_verification","spatial_statistics","dapc_analysis","haplogroup_analysis","rarefaction_sharing","run_outputs","output_bundle","plotting_functions","pipeline_stages"))source(file.path("R",paste0(f,".R")))
pipeline_order<-c("00_setup","01_import_validate","02_build_datasets","03_sampling_statistics","04_allele_frequencies",
  "05_forensic_statistics","06_sensitivity_analysis","07_marker_panels","08_haplotype_accumulation",
  "09a_prepare_structure","09b_export_arlequin","09c_import_arlequin","10_spatial_statistics","11_geographic_dapc",
  "12_haplogroup_analysis","13_rarefaction_sharing","14_bootstrap_uncertainty","15_generate_tables","16_generate_figures","17_final_qc","18_publication_package")
load_pipeline_config <- function(path) {
  cfg<-yaml::read_yaml(safe_path(path));cfg$config_path<-path
  assert(cfg$mode %in% c("synthetic","private"),"Choose synthetic or private mode")
  assert(grepl("^[A-Za-z0-9_-]+$",cfg$run_id),"Unsafe run ID")
  assert(identical(cfg$results_root,"results")&&identical(cfg$clean_root,"data/clean"),"Use repository-local result and clean-data roots")
  for(f in c("input","qc_config","panel_file","seed_file","prediction_leaves","parent_mapping","reported_parents"))safe_path(cfg[[f]])
  if(cfg$mode=="private"){
    assert(!any(grepl("example|synthetic",unlist(cfg[c("input","prediction_leaves","parent_mapping","reported_parents")]),ignore.case=TRUE)),"Private mode cannot use demonstration files")
  }
  cfg$out<-file.path("results",cfg$run_id);cfg$clean_out<-file.path("data/clean",cfg$run_id)
  cfg$seeds<-read.csv(cfg$seed_file);assert(!anyDuplicated(cfg$seeds$analysis),"Duplicate random seed key")
  assert(all(unlist(cfg$iterations)>=1)&&cfg$dapc$repetitions>=2,"Invalid iteration settings")
  cfg
}
stage_seed<-function(cfg,key){i<-match(key,cfg$seeds$analysis);assert(!is.na(i),paste("Missing seed",key));as.integer(cfg$seeds$seed[i])}
stage_write <- function(cfg,id,result) {
  out<-file.path(cfg$out,"intermediate");dir.create(out,recursive=TRUE,showWarnings=FALSE)
  saveRDS(result,file.path(out,paste0(id,".rds")))
  tables<-result[vapply(result,is.data.frame,TRUE)]
  dest<-file.path(cfg$out,"tables",id);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  if(length(tables)){
    for(n in names(tables))write_csv(tables[[n]],file.path(dest,paste0(n,".csv")))
    openxlsx::write.xlsx(tables,file.path(dest,"results.xlsx"),creator="Y-STR workflow",overwrite=TRUE)
  }
}
pipeline_input_hash <- function(cfg) {
  paths<-unique(c(cfg$config_path,cfg$input,cfg$qc_config,cfg$panel_file,cfg$seed_file,cfg$prediction_leaves,cfg$parent_mapping,cfg$reported_parents,"renv.lock",
    unlist(cfg[c("original_targets","original_observations","private_targets")],use.names=FALSE),
    list.files("R",recursive=TRUE,full.names=TRUE),list.files("config",pattern="\\.csv$",full.names=TRUE)))
  # Include private QC evidence and mappings referred to by the QC config.
  qe<-new.env(parent=globalenv());sys.source(cfg$qc_config,qe)
  assert(all(file.exists(paths)),"A frozen input, audit specification or source file is missing")
  digest::digest(list(files=setNames(unname(tools::md5sum(paths)),paths),qc=qe$qc_config),algo="sha256")
}
run_pipeline <- function(config_path="config/analysis_settings.yml",from=NULL,until=NULL,resume=FALSE) {
  require_packages()
  for(p in c("yaml","digest","adegenet"))assert(requireNamespace(p,quietly=TRUE),paste("Missing package",p))
  cfg<-load_pipeline_config(config_path);fingerprint<-pipeline_input_hash(cfg)
  grDevices::pdf(file=NULL);measure_device<-grDevices::dev.cur()
  on.exit(if(measure_device %in% grDevices::dev.list())grDevices::dev.off(measure_device),add=TRUE)
  path<-file.path(cfg$out,"pipeline_state.rds")
  if(resume){
    assert(file.exists(path),"No checkpoint to resume")
    saved<-readRDS(path);assert(identical(saved$input_hash,fingerprint),"Inputs/config/code changed; start a new run ID")
    state<-saved$state;done<-saved$completed
  }else{
    assert(!dir.exists(cfg$out),"Run output exists. Use --resume or choose a new run ID")
    dir.create(cfg$out,recursive=TRUE);state<-list();done<-character()
  }
  start<-if(is.null(from))1L else match(from,pipeline_order);end<-if(is.null(until))length(pipeline_order) else match(until,pipeline_order)
  assert(!is.na(start)&&!is.na(end)&&start<=end,"Unknown or reversed stage range")
  for(i in seq.int(start,end)){
    id<-pipeline_order[i]
    if(id %in% done)next
    assert(i==1L||all(pipeline_order[seq_len(i-1L)] %in% done),"Required preceding stages have not completed")
    message("Running ",id)
    warnings<-character();began<-Sys.time()
    value<-tryCatch(withCallingHandlers(dispatch_stage(id,state,cfg),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=identity)
    if(inherits(value,"error")){
      writeLines(conditionMessage(value),file.path(cfg$out,"STOPPED.txt"))
      saveRDS(list(state=state,completed=done,input_hash=fingerprint),path)
      stop(id,": ",conditionMessage(value),call.=FALSE)
    }
    state[[id]]<-value;stage_write(cfg,id,value);done<-c(done,id)
    saveRDS(list(state=state,completed=done,input_hash=fingerprint),path)
    log<-data.frame(Stage=id,Status="COMPLETED",Seconds=as.numeric(difftime(Sys.time(),began,units="secs")),Warnings=paste(unique(warnings),collapse=" | "))
    logpath<-file.path(cfg$out,"stage_log.csv")
    if(file.exists(logpath))log<-rbind(read.csv(logpath),log)
    write_csv(log,logpath)
  }
  message("Completed through ",tail(done,1),". See ",cfg$out,"/stage_log.csv")
  invisible(state)
}
run_stage_entry <- function(id) {
  args<-commandArgs(trailingOnly=TRUE);config<-if(length(args))args[1] else "config/analysis_settings.yml"
  cfg<-load_pipeline_config(config)
  run_pipeline(config,from=id,until=id,resume=file.exists(file.path(cfg$out,"pipeline_state.rds")))
}
