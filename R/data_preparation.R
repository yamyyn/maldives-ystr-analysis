# Count-neutral preparation matching the roles of original project-info and
# clean-data outputs. Confirmation evidence is supplied, never inferred.
enrich_initial_qc <- function(r,config) {
  c<-config;d<-r$clean;n<-nrow(d);p<-r$provenance
  conf<-c$confirmed_patterns
  classification<-rep("",n);evidence<-rep("",n)
  if(!is.null(conf)&&nrow(conf)){
    kind<-if("classification" %in% names(conf))conf$classification else rep("confirmed_no_call",nrow(conf))
    assert(all(kind %in% c("confirmed_deletion","confirmed_no_call")),"Confirmation classification must be confirmed_deletion or confirmed_no_call")
    ii<-match(d[[c$id]],conf$sample_id);okay<-!is.na(ii)&p$Confirmed_no_call
    classification[okay]<-kind[ii[okay]];evidence[okay]<-conf$evidence[ii[okay]]
  }
  long<-r$long
  missing_n<-vapply(strsplit(p$No_call_pattern,"|",fixed=TRUE),function(x)sum(nzchar(x)),0L)
  invalid_rows<-unique(long$Row[long$State=="invalid"])
  zero_rows<-unique(r$copy_audit$Row[r$copy_audit$Observed_slots==0L])
  no_multi<-unique(long$Row[long$Field %in% unlist(c$multicopy,use.names=FALSE)&long$State=="declared_no_call"])
  invalid<-seq_len(n) %in% c(invalid_rows,zero_rows,no_multi)
  status<-ifelse(invalid,"INVALID_OR_UNRESOLVED",ifelse(missing_n==0L,"COMPLETE",
    ifelse(classification=="confirmed_deletion","CONFIRMED_DELETION",
      ifelse(classification=="confirmed_no_call","CONFIRMED_OTHER_NO_CALL","PARTIAL_UNREVIEWED"))))
  p$atoll<-d[[c$population]];p$region<-d[[c$region]]
  p$Single_copy_missing_n<-missing_n;p$Profile_status<-status
  known<-if(!is.null(conf)&&nrow(conf))lapply(conf$loci,function(x)strsplit(x,"|",fixed=TRUE)[[1]]) else list()
  p$Pattern_relationship<-vapply(p$No_call_pattern,function(x){
    if(!nzchar(x))return("NO_SINGLE_COPY_NO_CALL")
    ff<-strsplit(x,"|",fixed=TRUE)[[1]]
    if(any(vapply(known,function(k)setequal(ff,k),TRUE)))"EXACT_CONFIGURED_PATTERN" else
      if(any(vapply(known,function(k)length(intersect(ff,k))>0L,TRUE)))"PARTIAL_OR_DISCORDANT_CONFIGURED_PATTERN" else "OTHER_UNREVIEWED_PATTERN"
  },"")
  p$Deletion_profile<-status=="CONFIRMED_DELETION";p$Confirmation_evidence<-evidence
  p$Microvariant_observations<-tabulate(long$Row[long$State=="microvariant"],nbins=n)
  p$Multicopy_systems_with_over_two_observations<-tabulate(r$copy_audit$Row[r$copy_audit$Observed_slots>2L],nbins=n)
  p$Invalid_observations<-tabulate(long$Row[long$State=="invalid"],nbins=n)
  # This describes observed completeness separately from release eligibility.
  p$Genotype_complete<-status=="COMPLETE"
  p$Complete_structure_eligible<-r$passed&p$Genotype_complete
  r$provenance<-p;r$profile_details<-p
  r$microvariants<-long[long$State=="microvariant",,drop=FALSE]
  r$microvariants$sample_id<-d[[c$id]][r$microvariants$Row]
  r$microvariant_summary<-if(nrow(r$microvariants))aggregate(list(Count=rep(1L,nrow(r$microvariants))),r$microvariants[c("Field","Normalized")],sum) else data.frame(Field=character(),Normalized=character(),Count=integer())
  r$nonstandard<-long[long$State=="invalid",,drop=FALSE]
  r$partial_details<-p[missing_n>0L|invalid,,drop=FALSE]
  r$deletion_details<-p[p$Deletion_profile,,drop=FALSE]
  r$blank_by_field<-do.call(rbind,lapply(unique(long$Field),function(f){
    z<-long[long$Field==f,];data.frame(Field=f,Profiles=nrow(z),Numeric_n=sum(z$State %in% c("integer","microvariant")),
      Blank_n=sum(z$State=="blank"),Declared_no_call_n=sum(z$State=="declared_no_call"),Invalid_n=sum(z$State=="invalid"))
  }))
  pattern_counts<-as.data.frame(table(Pattern=ifelse(nzchar(p$No_call_pattern),p$No_call_pattern,"<COMPLETE>")),stringsAsFactors=FALSE)
  names(pattern_counts)[2]<-"Profiles";r$pattern_counts<-pattern_counts
  r$copy_counts<-as.data.frame(table(System=r$copy_audit$System,Observed_slots=r$copy_audit$Observed_slots),stringsAsFactors=FALSE)
  valid_states<-r$copy_audit[!is.na(r$copy_audit$Unordered_state),,drop=FALSE]
  r$copy_states<-if(nrow(valid_states))aggregate(list(Profiles=rep(1L,nrow(valid_states))),valid_states[c("System","Observed_slots","Unordered_state")],sum) else
    data.frame(System=character(),Observed_slots=integer(),Unordered_state=character(),Profiles=integer())
  r$copy_exceptions<-r$copy_audit[r$copy_audit$Observed_slots!=2L,,drop=FALSE]
  r$locus_metadata<-data.frame(Field=c(c$single,unlist(c$multicopy,use.names=FALSE)),
    System=c(c$single,rep(names(c$multicopy),lengths(c$multicopy))),
    Class=c(rep("Single-copy",length(c$single)),rep("Multicopy observed slot",sum(lengths(c$multicopy)))))
  r$marker_definition<-data.frame(Marker_order=seq_len(length(c$single)+length(c$multicopy)),
    Analysis_column=c(c$single,names(c$multicopy)),Class=c(rep("Single-copy",length(c$single)),rep("Unordered observed state",length(c$multicopy))))
  summarize<-function(groups){
    partition<-split(seq_len(n),groups,drop=TRUE)
    if(!length(partition))return(data.frame(Group=character(),Full_n=integer(),Complete_n=integer(),Confirmed_deletion_n=integer(),Other_partial_or_invalid_n=integer(),Arithmetic_difference=integer()))
    do.call(rbind,lapply(partition,function(ii){
      data.frame(Group=as.character(groups[ii[1]]),Full_n=length(ii),Complete_n=sum(p$Genotype_complete[ii]),
        Confirmed_deletion_n=sum(p$Deletion_profile[ii]),Other_partial_or_invalid_n=sum(!p$Genotype_complete[ii]&!p$Deletion_profile[ii]),
        Arithmetic_difference=length(ii)-sum(p$Genotype_complete[ii])-sum(p$Deletion_profile[ii])-sum(!p$Genotype_complete[ii]&!p$Deletion_profile[ii]))
    }))
  }
  r$counts_by_population<-summarize(d[[c$population]]);r$counts_by_region<-summarize(d[[c$region]])
  r$dataset_summary<-data.frame(Dataset=c("Full archival","Full forensic","Complete","Confirmed deletion","Other partial or invalid"),
    N=c(n,if(r$passed)n else 0L,sum(p$Genotype_complete),sum(p$Deletion_profile),sum(!p$Genotype_complete&!p$Deletion_profile)),
    Released=c(TRUE,r$passed,r$passed,r$passed,FALSE))
  r$deletion<-if(r$passed) {
    z<-d[p$Deletion_profile,,drop=FALSE]
    for(f in names(c$multicopy))z[[f]]<-r$full[[f]][p$Deletion_profile]
    z
  } else NULL
  r$complete<-if(r$passed)r$full[p$Genotype_complete,,drop=FALSE] else NULL
  # Every confirmed no-call observation, including pre-coded tokens, is audited.
  rr<-which(p$Confirmed_no_call)
  r$recoding_detail<-if(length(rr))do.call(rbind,lapply(rr,function(i){
    ff<-strsplit(p$No_call_pattern[i],"|",fixed=TRUE)[[1]]
    original<-unlist(d[i,ff,drop=FALSE],use.names=FALSE)
    data.frame(Row=i,sample_id=d[[c$id]][i],Locus=ff,Source_value=original,
      Derived_value=c$no_call_token,Change=ifelse(is.na(original)|original=="","BLANK_TO_TOKEN","PREEXISTING_TOKEN"),
      Classification=classification[i],Evidence=evidence[i],Valid=(!is.na(original)&original==c$no_call_token)|is.na(original)|original=="")
  })) else data.frame(Row=integer(),sample_id=character(),Locus=character(),Source_value=character(),Derived_value=character(),Change=character(),Classification=character(),Evidence=character(),Valid=logical())
  r
}

write_prepared_project <- function(r,out) {
  if(is.null(r$profile_details))return(invisible(NULL))
  info<-file.path(out,"00_project_info");clean<-file.path(out,"04_clean_data")
  dir.create(info,recursive=TRUE);dir.create(clean,recursive=TRUE)
  writebook<-function(name,tables)openxlsx::write.xlsx(tables,file.path(info,name),overwrite=FALSE)
  writebook("01_Import_QC_Summary.xlsx",list(Import_Summary=r$checks,Column_Mapping=r$column_map,Normalization_Changes=r$changes))
  writebook("02_Metadata_QC.xlsx",list(QC_Summary=r$checks,Counts_By_Population=r$counts_by_population,Counts_By_Region=r$counts_by_region,Coordinate_QC=r$coordinates,Issues=r$issues))
  writebook("03_YSTR_QC.xlsx",list(Locus_Metadata=r$locus_metadata,Blank_By_Field=r$blank_by_field,
    Profile_Details=r$profile_details,SingleCopy_BlankPatterns=r$pattern_counts,Deletion_Details=r$deletion_details,
    Partial_Profile_Details=r$partial_details,Microvariant_Observations=r$microvariants,Nonstandard_Values=r$nonstandard,
    Microvariant_Summary=r$microvariant_summary,
    Allele_Range_Summary=r$observed_ranges,Multicopy_Counts=r$copy_counts,Multicopy_States=r$copy_states,
    Multicopy_Exceptions=r$copy_exceptions,Exact_Value_Inventory=r$long))
  writebook("04_Analysis_Dataset_Creation_QC.xlsx",list(Whole_Dataset_Summary=r$dataset_summary,
    Counts_By_Population=r$counts_by_population,Counts_By_Region=r$counts_by_region,Sample_Provenance=r$provenance,
    DEL_Recoding_Detail=r$recoding_detail,Forensic_Marker_Definition=r$marker_definition))
  for(n in c("dataset_summary","counts_by_population","counts_by_region","recoding_detail","profile_details"))write_csv(r[[n]],file.path(info,paste0(n,".csv")))
  export<-function(x,folder,stem){
    dest<-file.path(clean,folder);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    write_csv(x,file.path(dest,paste0(stem,".csv")));saveRDS(x,file.path(dest,paste0(stem,".rds")))
    openxlsx::write.xlsx(list(Data=x),file.path(dest,paste0(stem,".xlsx")),overwrite=FALSE)
  }
  export(r$raw,"full","Full_Archival_Unchanged")
  export(r$clean,"full","Master_Normalized_For_Review")
  # Review material is exported even on failure, explicitly outside analysis cohorts.
  export(r$partial_details,"review","Partial_And_Invalid_Profile_Details")
  export(r$provenance,"review","Sample_Analysis_Provenance")
  if(r$passed){
    export(r$full,"full","Forensic_Full_Analysis_Ready")
    export(r$complete,"complete","Complete_Analysis_Ready")
    export(r$deletion,"deletion","Confirmed_Deletion_Validated")
    openxlsx::write.xlsx(list(Full_Archival=r$raw,Forensic_Full=r$full,Complete=r$complete,
      Confirmed_Deletion=r$deletion,Sample_Provenance=r$provenance),file.path(clean,"YSTR_Analysis_Datasets.xlsx"),overwrite=FALSE)
  }
  files<-list.files(out,recursive=TRUE,full.names=TRUE)
  files<-files[!dir.exists(files)]
  write_csv(data.frame(File=substring(files,nchar(out)+2L),Bytes=file.info(files)$size,MD5=unname(tools::md5sum(files))),file.path(info,"Preparation_File_Manifest.csv"))
  writeLines(c("Archival and review files are not certified analysis cohorts.",
    "Analysis-ready files exist only when configured QC passes.",
    "Confirmed deletion requires explicit classification and independent evidence.",
    "Multicopy counts describe observed slots, not inferred genomic copies.",
    "Structure-specific DYS389II transformation and Arlequin integer recoding remain downstream."),file.path(clean,"README.txt"))
}
