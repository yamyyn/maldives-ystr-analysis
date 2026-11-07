source("R/initial_qc.R")
qe<-new.env(parent=globalenv());sys.source("config/initial_qc.R",qe)
raw<-read_qc_csv(qe$qc_input);config<-qe$qc_config
r<-run_initial_qc(raw,config)
check(r$passed)
check(identical(r$raw,raw))
check(nrow(r$full)==nrow(raw))
check(nrow(r$complete)+sum(r$provenance$Confirmed_no_call)==nrow(raw))
check(all(r$full[r$provenance$Confirmed_no_call,config$single][r$full[r$provenance$Confirmed_no_call,config$single]=="DEL"]=="DEL"))
check(nrow(r$recoding)>0L&&all(r$recoding$Derived=="DEL"))
check(all(r$provenance$Full_forensic_eligible))
bad<-raw;bad$sample_id[2]<-bad$sample_id[1]
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad$sample_id[1]<-" "
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad[[config$single[1]]][1]<-"12x"
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad[[config$single[1]]][1]<-"12.2"
micro<-run_initial_qc(bad,config)
check(micro$passed&&micro$full[[config$single[1]]][1]=="12.2")
bad<-raw;bad[[config$single[1]]][1]<-" 12.2 "
trimmed<-run_initial_qc(bad,config)
check(trimmed$passed&&nrow(trimmed$changes)==1L&&trimmed$raw[[config$single[1]]][1]==" 12.2 ")
bad<-raw;bad$dys385b[1]<-""
one<-run_initial_qc(bad,config)
check(one$passed&&one$copy_audit$Observed_slots[1]==1L)
bad<-raw;bad[1,c("dys385a","dys385b")]<-""
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad$dys385_state<-"invalid-derived-state"
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad$region[1]<-"Unexpected"
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad$latitude_centroid[1]<-"100"
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad$total_samples[1]<-"999"
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad$islands_sampled[1]<-"1.5"
check(!run_initial_qc(bad,config)$passed)
bad<-raw;bad[[config$single[1]]][1]<-""
blocked<-run_initial_qc(bad,config)
check(!blocked$passed&&is.null(blocked$full)&&is.null(blocked$complete))
conf<-config;conf$confirmed_patterns<-NULL
check(!run_initial_qc(raw,conf)$passed)
conf<-config;conf$confirmed_patterns$evidence[1]<-""
check(!run_initial_qc(raw,conf)$passed)
conf<-config;conf$confirmed_patterns$loci[1]<-config$single[1]
check(!run_initial_qc(raw,conf)$passed)
conf<-config;conf$expected_n<-nrow(raw)+1L
check(!run_initial_qc(raw,conf)$passed)
conf<-config;conf$allele_bounds<-data.frame(field=config$single[1],min=100,max=200)
check(!run_initial_qc(raw,conf)$passed)
conf<-config;conf$id_pattern<-"^NONMATCH"
check(!run_initial_qc(raw,conf)$passed)
bad<-raw;bad[[config$single[1]]]<-NULL
check(!run_initial_qc(bad,config)$passed)
temp<-file.path("outputs",paste0("qc-test-",format(Sys.time(),"%Y%m%d%H%M%S")))
write_initial_qc(blocked,temp)
check(file.exists(file.path(temp,"Initial_QC.xlsx"))&&!file.exists(file.path(temp,"forensic_full.rds")))
check(reject(write_initial_qc(r,temp)))
renamed<-raw;names(renamed)[names(renamed)=="sample_id"]<-"Raw Sample ID"
conf<-config;conf$column_map<-data.frame(source="Raw Sample ID",target="sample_id")
mapped<-run_initial_qc(renamed,conf)
check(mapped$passed&&"Raw Sample ID" %in% names(mapped$raw)&&"sample_id" %in% names(mapped$full))
conf$column_map<-data.frame(source="Raw Sample ID",target="atoll")
check(reject(run_initial_qc(renamed,conf)))
swapped<-raw;swapped$dys385a<-raw$dys385b;swapped$dys385b<-raw$dys385a
check(identical(run_initial_qc(swapped,config)$full$dys385_state,r$full$dys385_state))

# One observed size and repeated equal-size slots have the same analytical state.
homo<-raw;homo$dys385b[1]<-homo$dys385a[1]
onecopy<-homo;onecopy$dys385b[1]<-""
homo_qc<-run_initial_qc(homo,config);one_qc<-run_initial_qc(onecopy,config)
check(homo_qc$passed && one_qc$passed)
check(homo_qc$full$dys385_state[1]==homo$dys385a[1])
check(identical(homo_qc$full$dys385_state[1],one_qc$full$dys385_state[1]))
check(identical(homo_qc$raw,homo) && identical(one_qc$raw,onecopy))
check(identical(homo_qc$full$dyf387s1_state,r$full$dyf387s1_state))

# Downstream table preparation must use the same DYS385 analytical states.
slots<-config$multicopy$dys385_state
table_helper<-method_env("table_alleles",list(DYS385_COLUMNS=slots))
check(identical(table_helper$make_pattern(homo_qc$full,slots),homo_qc$full$dys385_state))
check(identical(table_helper$make_pattern(one_qc$full,slots),one_qc$full$dys385_state))
locus_helper<-method_env("analysis_locus")
check(identical(unname(locus_helper$make_unordered_pattern(homo_qc$full,slots)),unname(homo_qc$full$dys385_state)))
