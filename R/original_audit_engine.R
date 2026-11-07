# Runs the original project's equality/near/file helpers against explicitly
# supplied private observations and an independent private target specification.
run_original_audit_plan <- function(plan,observations=NULL,targets=NULL,style,study) {
  if(!is.null(targets))assert(!anyDuplicated(targets$key),"Duplicate private target key")
  if(!is.null(observations))assert(!anyDuplicated(observations$key),"Duplicate private observation key")
  e<-recipe_environment("project_audit",list(audit_rows=list()),style,study)
  absent<-list()
  for(i in seq_len(nrow(plan))) {
    item<-plan[i,];key<-item$key
    ti<-if(is.null(targets))integer() else which(targets$key==key)
    oi<-if(is.null(observations))integer() else which(observations$key==key)
    missing_target<-!length(ti)||is.na(targets$expected[ti])||!nzchar(targets$expected[ti])
    missing_observed<-!length(oi)||is.na(observations$value[oi])
    if(missing_target||missing_observed) {
      absent[[length(absent)+1L]]<-data.frame(Section=item$section,Check=item$check,Expected="private target required",
        Observed=if(missing_observed)"not observed" else observations$value[oi],
        Status=if(missing_target)"NOT_CONFIGURED" else "NOT_OBSERVED",Critical=item$critical,
        Evidence=if(missing_observed)NA_character_ else observations$evidence[oi],Notes=paste("Source check key",key))
      next
    }
    t<-targets[ti,];o<-observations[oi,]
    assert(length(ti)==1L&&length(oi)==1L,"Duplicate private rule or observation key")
    rule<-t$comparison
    if(rule=="equal")e$check_equal(item$section,item$check,o$value,t$expected,item$critical,o$evidence)
    else if(rule=="near") {
      value<-suppressWarnings(as.numeric(o$value));expected<-suppressWarnings(as.numeric(t$expected));tol<-suppressWarnings(as.numeric(t$tolerance))
      assert(is.finite(expected)&&is.finite(tol)&&tol>=0,"Invalid private near target/tolerance")
      e$check_near(item$section,item$check,value,expected,tol,item$critical,o$evidence)
    }else if(rule=="file_exists")e$check_file(item$section,item$check,o$value,item$critical)
    else if(rule %in% c("at_least","at_most")) {
      value<-suppressWarnings(as.numeric(o$value));expected<-suppressWarnings(as.numeric(t$expected))
      pass<-is.finite(value)&&is.finite(expected)&&if(rule=="at_least")value>=expected else value<=expected
      e$check_equal(item$section,item$check,pass,TRUE,item$critical,o$evidence,
        paste("Explicit private rule",rule,"observed",o$value,"threshold",t$expected))
    }else stop("Unsupported private comparison for key ",key,call.=FALSE)
  }
  rows<-c(e$audit_rows,absent)
  if(length(rows))do.call(rbind,rows) else data.frame()
}

write_original_project_audit <- function(path,audit,targets=NULL,style,study) {
  # Preserve the original workbook's Audit/Status_Summary/Section_Summary/
  # Locked_Values/Run_Info architecture and formatting. Incomplete states are
  # additional fail-closed states; they are never counted as successful checks.
  if(!"Evidence" %in% names(audit))audit$Evidence<-NA_character_
  if(!"Notes" %in% names(audit))audit$Notes<-NA_character_
  status_summary<-as.data.frame(table(Status=audit$Status));names(status_summary)[2]<-"N"
  sections<-unique(audit$Section)
  section_summary<-do.call(rbind,lapply(sections,function(section){
    x<-audit[audit$Section==section,]
    data.frame(Section=section,PASS=sum(x$Status=="PASS"),WARN=sum(x$Status=="WARN"),FAIL=sum(x$Status=="FAIL"),
      NOT_CONFIGURED=sum(x$Status=="NOT_CONFIGURED"),NOT_OBSERVED=sum(x$Status=="NOT_OBSERVED"),
      Critical_FAIL=sum(x$Status=="FAIL"&x$Critical))
  }))
  critical_fail_n<-sum(audit$Status=="FAIL"&audit$Critical)
  incomplete<-sum(audit$Status %in% c("NOT_CONFIGURED","NOT_OBSERVED")&audit$Critical)
  overall_status<-if(critical_fail_n)"FAIL - CRITICAL ISSUES REQUIRE RESOLUTION" else
    if(incomplete)"INCOMPLETE - PRIVATE VALIDATION REQUIRED" else if(any(audit$Status=="WARN"))"PASS WITH WARNINGS" else "PASS"
  locked_values<-if(is.null(targets))data.frame(Domain="Private baseline",Metric="No target values included",Locked_Value="NOT_CONFIGURED") else targets
  e<-recipe_environment("project_audit",list(audit=audit,status_summary=status_summary,section_summary=section_summary,
    locked_values=locked_values,overall_status=overall_status,critical_fail_n=critical_fail_n,
    SCRIPT_NAME="source-faithful validation",root="."),style,study)
  e$write_workbook(path)
  invisible(overall_status)
}
