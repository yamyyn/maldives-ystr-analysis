r<-run_initial_qc(raw,config)
check(sum(r$provenance$Deletion_profile)>0L)
check(nrow(r$deletion)==sum(r$provenance$Deletion_profile))
check(all(r$deletion$sample_id %in% r$full$sample_id))
check(!any(r$complete$sample_id %in% r$deletion$sample_id))
check(setequal(c(r$complete$sample_id,r$deletion$sample_id),r$full$sample_id))
check(all(r$counts_by_population$Arithmetic_difference==0L)&&all(r$counts_by_region$Arithmetic_difference==0L))
check(sum(r$counts_by_population$Full_n)==nrow(raw))
check(all(r$recoding_detail$Valid)&&nrow(r$recoding_detail)==nrow(r$recoding))
for(i in seq_len(nrow(r$recoding_detail))){
 z<-r$recoding_detail[i,]
 check(is.na(r$deletion[r$deletion$sample_id==z$sample_id,z$Locus])||r$deletion[r$deletion$sample_id==z$sample_id,z$Locus]=="")
 check(r$full[r$full$sample_id==z$sample_id,z$Locus]=="DEL")
}
conf<-config;conf$confirmed_patterns$classification<-NULL
not_deletion<-run_initial_qc(raw,conf)
check(not_deletion$passed&&nrow(not_deletion$deletion)==0L)
check(all(not_deletion$provenance$Profile_status[not_deletion$provenance$Confirmed_no_call]=="CONFIRMED_OTHER_NO_CALL"))
check(all(not_deletion$counts_by_population$Arithmetic_difference==0L))
conf<-config;conf$confirmed_patterns$classification[1]<-"automatic_deletion"
check(reject(run_initial_qc(raw,conf)))
bad<-raw;bad[[config$single[1]]][1]<-"11.2"
micro<-run_initial_qc(bad,config)
check(any(micro$microvariants$Normalized=="11.2")&&any(micro$microvariant_summary$Normalized=="11.2"))
bad<-raw;bad[[config$single[1]]][1]<-"bad"
check(run_initial_qc(bad,config)$profile_details$Profile_status[1]=="INVALID_OR_UNRESOLVED")
dest<-file.path("outputs",paste0("prepared-test-",format(Sys.time(),"%Y%m%d%H%M%S")))
write_initial_qc(r,dest)
check(file.exists(file.path(dest,"00_project_info/03_YSTR_QC.xlsx")))
check(file.exists(file.path(dest,"04_clean_data/deletion/Confirmed_Deletion_Validated.rds")))
check(identical(readRDS(file.path(dest,"04_clean_data/full/Full_Archival_Unchanged.rds")),raw))
check(identical(readRDS(file.path(dest,"04_clean_data/complete/Complete_Analysis_Ready.rds")),r$complete))
check(nrow(openxlsx::read.xlsx(file.path(dest,"00_project_info/04_Analysis_Dataset_Creation_QC.xlsx"),sheet="Sample_Provenance"))==nrow(raw))
check(all(c("Microvariant_Observations","Multicopy_Exceptions","Deletion_Details") %in% openxlsx::getSheetNames(file.path(dest,"00_project_info/03_YSTR_QC.xlsx"))))
