# Configurable adaptation of the import, metadata, genotype and cohort checks
# in original scripts 01-04. No original cohort counts, ID patterns or records.
source("R/data_preparation.R")
initial_qc_defaults <- function(single, multicopy) {
  list(single=single,multicopy=multicopy,id="sample_id",population="atoll",region="region",
    trim_whitespace=TRUE,id_pattern=NULL,expected_ids=NULL,expected_n=NULL,
    expected_populations=NULL,expected_regions=NULL,population_mapping=NULL,
    metadata_required=TRUE,metadata_constant=c("latitude_centroid","longitude_centroid",
      "latitude_weighted","longitude_weighted","islands_sampled","total_samples"),
    allele_bounds=NULL,confirmed_patterns=NULL,no_call_token="DEL",column_map=NULL)
}

read_qc_csv <- function(path) {
  # Import all fields as character; no automatic allele/date/ID conversion.
  x<-read.csv(safe_path(path),colClasses="character",check.names=FALSE,
    na.strings=character(),strip.white=FALSE,stringsAsFactors=FALSE)
  assert(nrow(x)>0L,"Input contains no profiles")
  assert(!anyDuplicated(names(x)),"Duplicate input column names")
  x
}

run_initial_qc <- function(raw,config) {
  c<-config;d<-as.data.frame(raw,stringsAsFactors=FALSE,check.names=FALSE)
  assert(nrow(d)>0L&&!anyDuplicated(names(d)),"Empty input or duplicate columns")
  column_map<-data.frame(Original=names(d),Canonical=names(d))
  if(!is.null(c$column_map)){
    m<-c$column_map
    assert(all(c("source","target") %in% names(m))&&!anyDuplicated(m$source)&&!anyDuplicated(m$target)&&
      all(m$source %in% names(d))&&all(!is.na(m$target)&nzchar(m$target)),"Invalid explicit column mapping")
    names(d)[match(m$source,names(d))]<-m$target
    assert(!anyDuplicated(names(d)),"Column mapping creates duplicate canonical names")
    column_map$Canonical<-names(d)
  }
  assert(length(c$single)>0L&&length(c$multicopy)>0L&&!is.null(names(c$multicopy)),"Define single-copy fields and named multicopy slot lists")
  slots<-unlist(c$multicopy,use.names=FALSE);fields<-c(c$single,slots)
  assert(!anyDuplicated(fields),"Allele fields must belong to exactly one system")
  assert(length(c$no_call_token)==1L&&nzchar(c$no_call_token)&&!grepl("^[0-9]",c$no_call_token),"Invalid no-call token")
  issues<-list();changes<-list();checks<-list()
  add<-function(severity,check,rows=integer(),field="",detail="") {
    if(!length(rows))rows<-NA_integer_
    issues[[length(issues)+1L]]<<-data.frame(Severity=severity,Check=check,Row=rows,Field=field,Detail=detail)
  }
  record<-function(check,pass,detail="",severity="ERROR") {
    checks[[length(checks)+1L]]<<-data.frame(Check=check,Status=if(is.na(pass))"NOT_CONFIGURED" else if(pass)"PASS" else "FAIL",Detail=detail)
    if(!is.na(pass)&&!pass)add(severity,check,detail=detail)
  }
  blank<-function(x)is.na(x)|x==""
  numeric_allele<-function(x)!blank(x)&grepl("^[0-9]+(?:\\.[0-9]+)?$",x,perl=TRUE)
  required<-c(c$id,c$population,c$region,fields)
  if(c$metadata_required)required<-c(required,c$metadata_constant)
  missing<-setdiff(required,names(d));record("Required columns",!length(missing),paste(missing,collapse=";"))
  if(length(missing))return(list(raw=raw,clean=d,checks=do.call(rbind,checks),issues=do.call(rbind,issues),passed=FALSE))
  for(f in names(d)){
    old<-as.character(d[[f]]);value<-if(c$trim_whitespace)trimws(old) else old
    changed<-which(!is.na(old)&old!=value)
    if(length(changed))changes[[length(changes)+1L]]<-data.frame(Row=changed,Field=f,Original=old[changed],Normalized=value[changed])
    d[[f]]<-value
  }
  ids<-d[[c$id]]
  bad_id<-which(blank(ids)|duplicated(ids)|duplicated(ids,fromLast=TRUE))
  record("Unique nonblank identifiers",!length(bad_id));if(length(bad_id))add("ERROR","Invalid identifier",bad_id,c$id)
  for(f in c(c$population,c$region))if(any(blank(d[[f]])))add("ERROR","Missing grouping metadata",which(blank(d[[f]])),f)
  record("Configured identifier format",if(is.null(c$id_pattern))NA else all(grepl(c$id_pattern,ids)),"No study-specific ID format inferred")
  record("Configured identifier set",if(is.null(c$expected_ids))NA else setequal(ids,c$expected_ids))
  record("Configured profile count",if(is.null(c$expected_n))NA else nrow(d)==c$expected_n)
  record("Configured population set",if(is.null(c$expected_populations))NA else setequal(unique(d[[c$population]]),c$expected_populations))
  record("Configured region set",if(is.null(c$expected_regions))NA else setequal(unique(d[[c$region]]),c$expected_regions))
  groups<-split(seq_len(nrow(d)),d[[c$population]],drop=TRUE)
  for(pop in names(groups)){
    ii<-groups[[pop]]
    if(length(unique(d[[c$region]][ii]))!=1L)add("ERROR","Population maps to multiple regions",ii,c$region,pop)
    for(f in intersect(c$metadata_constant,names(d)))if(any(blank(d[[f]][ii]))||length(unique(d[[f]][ii]))!=1L)
      add("ERROR","Population metadata not constant or missing",ii,f,pop)
    if("total_samples" %in% names(d)){
      v<-suppressWarnings(as.numeric(d$total_samples[ii]))
      if(any(!is.finite(v)|v!=length(ii)))add("ERROR","Recorded population total differs from rows",ii,"total_samples",pop)
    }
    if("islands_sampled" %in% names(d)){
      v<-suppressWarnings(as.numeric(d$islands_sampled[ii]))
      if(any(!is.finite(v)|v<1|v>length(ii)|v!=floor(v)))add("ERROR","Invalid islands sampled",ii,"islands_sampled",pop)
    }
  }
  if(!is.null(c$population_mapping)){
    m<-c$population_mapping
    assert(all(c("population","region") %in% names(m))&&!anyDuplicated(m$population),"Invalid configured population mapping")
    expected<-m$region[match(d[[c$population]],m$population)]
    bad<-which(is.na(expected)|expected!=d[[c$region]])
    record("Configured population-region mapping",!length(bad));if(length(bad))add("ERROR","Population-region mismatch",bad,c$region)
  }else record("Configured population-region mapping",NA)
  coords<-c("latitude_centroid","longitude_centroid","latitude_weighted","longitude_weighted")
  for(f in intersect(coords,names(d))){
    v<-suppressWarnings(as.numeric(d[[f]]));lim<-if(grepl("^latitude",f))90 else 180
    bad<-which(!is.finite(v)|abs(v)>lim)
    if(length(bad))add("ERROR","Invalid geographic coordinate",bad,f)
  }
  distances<-data.frame()
  if(all(coords %in% names(d))){
    e<-recipe_environment("metadata_qc_helpers",list(),list(),list())
    v<-lapply(d[coords],function(x)suppressWarnings(as.numeric(x)))
    distances<-data.frame(Row=seq_len(nrow(d)),Centroid_to_weighted_km=do.call(e$haversine_km,unname(v)))
  }
  long<-do.call(rbind,lapply(fields,function(f){
    x<-d[[f]];num<-numeric_allele(x);no<-!blank(x)&x==c$no_call_token
    data.frame(Row=seq_len(nrow(d)),Field=f,Raw=as.character(raw[[column_map$Original[match(f,column_map$Canonical)]]]),Normalized=x,
      State=ifelse(blank(x),"blank",ifelse(no,"declared_no_call",ifelse(num,ifelse(grepl(".",x,fixed=TRUE),"microvariant","integer"),"invalid"))))
  }))
  invalid<-long[long$State=="invalid",]
  if(nrow(invalid))for(i in seq_len(nrow(invalid)))add("ERROR","Non-numeric allele string",invalid$Row[i],invalid$Field[i])
  if(!is.null(c$allele_bounds)){
    bounds<-c$allele_bounds
    assert(all(c("field","min","max") %in% names(bounds))&&!anyDuplicated(bounds$field),"Invalid allele bounds")
    assert(all(is.finite(bounds$min)&is.finite(bounds$max)&bounds$min<=bounds$max),"Invalid allele-bound limits")
    for(i in seq_len(nrow(bounds))){
      f<-bounds$field[i];assert(f %in% fields,"Unknown bounded allele field")
      v<-suppressWarnings(as.numeric(d[[f]]));bad<-which(is.finite(v)&(v<bounds$min[i]|v>bounds$max[i]))
      if(length(bad))add("ERROR","Allele outside configured range",bad,f)
    }
  }else record("Laboratory-specific allele bounds",NA,"Observed ranges are descriptive, not acceptance limits")
  # All no-call patterns are inventoried. Biological confirmation is independent.
  missing_fields<-lapply(seq_len(nrow(d)),function(i)c$single[vapply(d[i,c$single,drop=FALSE],function(x)blank(x)||x==c$no_call_token,TRUE)])
  pattern<-vapply(missing_fields,function(x)paste(sort(x),collapse="|"),"")
  confirmed<-rep(FALSE,nrow(d));recoding<-list()
  if(!is.null(c$confirmed_patterns)&&nrow(c$confirmed_patterns)){
    conf<-c$confirmed_patterns
    assert(all(c("sample_id","loci","evidence") %in% names(conf))&&!anyDuplicated(conf$sample_id),"Invalid confirmed-pattern specification")
    for(i in seq_len(nrow(conf))){
      ii<-which(ids==conf$sample_id[i]);loc<-strsplit(conf$loci[i],"|",fixed=TRUE)[[1]]
      evidence<-!is.na(conf$evidence[i])&&nzchar(trimws(conf$evidence[i]))
      okay<-length(ii)==1L&&evidence&&length(loc)>0L&&!anyDuplicated(loc)&&all(loc %in% c$single)&&setequal(loc,missing_fields[[ii]])
      if(!okay){add("ERROR","Confirmation does not match exact observed no-call pattern",field="confirmation",detail=paste("Specification row",i));next}
      confirmed[ii]<-TRUE
      for(f in loc)if(blank(d[[f]][ii])){
        recoding[[length(recoding)+1L]]<-data.frame(Row=ii,Field=f,Original=d[[f]][ii],Derived=c$no_call_token,Evidence=conf$evidence[i])
      }
    }
  }
  for(i in which(nzchar(pattern)&!confirmed))add("ERROR","Unreviewed single-copy no-call pattern",i,detail=pattern[i])
  e<-recipe_environment("initial_qc_helpers",list(),list(),list())
  copy_audit<-list();states<-list()
  for(system in names(c$multicopy)){
    ff<-c$multicopy[[system]]
    observed<-vapply(seq_len(nrow(d)),function(i)sum(numeric_allele(unlist(d[i,ff,drop=FALSE],use.names=FALSE))),0L)
    for(i in which(observed==0L))add("ERROR","Multicopy system has no numeric observation",i,system)
    no<-which(vapply(seq_len(nrow(d)),function(i)any(unlist(d[i,ff,drop=FALSE],use.names=FALSE)==c$no_call_token,na.rm=TRUE),TRUE))
    if(length(no))add("ERROR","No-call token in multicopy observation slots requires review",no,system)
    states[[system]]<-vapply(seq_len(nrow(d)),function(i)observed_multicopy_state(unlist(d[i,ff,drop=FALSE],use.names=FALSE),collapse_identical=identical(system,"dys385_state")),"")
    # Optional slots can be blank: observed peaks are not inferred chromosome copies.
    copy_audit[[system]]<-data.frame(Row=seq_len(nrow(d)),System=system,Observed_slots=observed,Unordered_state=states[[system]])
    if(system %in% names(d)){
      bad<-which(is.na(d[[system]])!=is.na(states[[system]])|(!is.na(d[[system]])&!is.na(states[[system]])&d[[system]]!=states[[system]]))
      if(length(bad))add("ERROR","Existing multicopy state disagrees with observed slots",bad,system)
    }
  }
  issue_table<-if(length(issues))do.call(rbind,issues) else data.frame(Severity=character(),Check=character(),Row=integer(),Field=character(),Detail=character())
  passed<-!any(issue_table$Severity=="ERROR")
  record("Configured metadata coverage",if(c$metadata_required)TRUE else NA,"Metadata coverage must be reviewed if disabled")
  record("Initial QC error count",passed,paste(sum(issue_table$Severity=="ERROR"),"error records"))
  forensic<-d
  for(i in which(confirmed))for(f in missing_fields[[i]])forensic[[f]][i]<-c$no_call_token
  for(f in names(states))forensic[[f]]<-states[[f]]
  provenance<-data.frame(Row=seq_len(nrow(d)),sample_id=ids,No_call_pattern=pattern,
    Confirmed_no_call=confirmed,Full_forensic_eligible=passed,Complete_structure_eligible=passed&!nzchar(pattern),
    Decision=if(!passed)"BLOCKED_QC_FAILURE" else ifelse(nzchar(pattern),"FULL_ONLY_CONFIRMED_NO_CALL","FULL_AND_COMPLETE"))
  inventory<-as.data.frame(table(Field=long$Field,State=long$State),stringsAsFactors=FALSE)
  observed_ranges<-do.call(rbind,lapply(fields,function(f){v<-suppressWarnings(as.numeric(d[[f]]));v<-v[is.finite(v)];data.frame(Field=f,Typed_n=length(v),Min=if(length(v))min(v) else NA_real_,Max=if(length(v))max(v) else NA_real_)}))
  fingerprint<-apply(d[,fields,drop=FALSE],1,function(x)paste(ifelse(is.na(x),"<NA>",x),collapse="|"))
  repeated<-data.frame(Row=which(duplicated(fingerprint)|duplicated(fingerprint,fromLast=TRUE)))
  result<-list(raw=raw,clean=d,column_map=column_map,repeated_profiles=repeated,checks=do.call(rbind,checks),issues=issue_table,passed=passed,long=long,
    inventory=inventory,observed_ranges=observed_ranges,copy_audit=do.call(rbind,copy_audit),coordinates=distances,provenance=provenance,
    changes=if(length(changes))do.call(rbind,changes) else data.frame(Row=integer(),Field=character(),Original=character(),Normalized=character()),
    recoding=if(length(recoding))do.call(rbind,recoding) else data.frame(Row=integer(),Field=character(),Original=character(),Derived=character(),Evidence=character()),
    full=if(passed)forensic else NULL,complete=if(passed)forensic[!nzchar(pattern),,drop=FALSE] else NULL)
  enrich_initial_qc(result,c)
}

write_initial_qc <- function(result,out,input_path=NULL) {
  assert(!file.exists(out),"QC output already exists; choose a fresh folder")
  dir.create(out,recursive=TRUE)
  report_names<-intersect(c("checks","issues","column_map","long","inventory","observed_ranges","copy_audit","coordinates","provenance","changes","recoding","repeated_profiles"),names(result))
  tables<-result[report_names]
  for(n in names(tables))write_csv(tables[[n]],file.path(out,paste0(n,".csv")))
  openxlsx::write.xlsx(tables,file.path(out,"Initial_QC.xlsx"),overwrite=FALSE)
  saveRDS(result$raw,file.path(out,"preserved_import.rds"))
  if(result$passed){saveRDS(result$full,file.path(out,"forensic_full.rds"));saveRDS(result$complete,file.path(out,"structure_complete.rds"))}
  writeLines(c(if(result$passed)"PASS: configured initial QC completed" else "BLOCKED: review issues; analytical datasets not released",
    "Optional unconfigured rules are not certified. Review all reports before downstream analysis.",
    if(!is.null(input_path))paste("Input MD5:",unname(tools::md5sum(input_path))) else "Input provided in memory"),file.path(out,"STATUS.txt"))
  write_prepared_project(result,out)
  invisible(result$passed)
}
