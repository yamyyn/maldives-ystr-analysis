if(!file.exists("R/pipeline.R"))stop("Run from the repository root")
source("R/pipeline.R")
args<-commandArgs(trailingOnly=TRUE)
get_arg<-function(flag,default=NULL){i<-match(flag,args);if(is.na(i))default else {if(i==length(args))stop("Missing value for ",flag);args[i+1L]}}
if("--list" %in% args){cat(paste(pipeline_order,collapse="\n"),"\n")}else{
  run_pipeline(get_arg("--config","config/analysis_settings.yml"),get_arg("--from"),get_arg("--until"),"--resume" %in% args)
}
