run_source_outputs <- function(bundle,cfg,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(folder in c("figures","tables","qc","inventory","captions"))dir.create(file.path(out,folder),showWarnings=FALSE)
  assert(!is.null(cfg$initial_qc),"Configure initial QC before downstream output generation")
  if(!isTRUE(cfg$skip_initial_qc)){
    qe<-new.env(parent=globalenv());sys.source(safe_path(cfg$initial_qc),qe)
    initial<-run_initial_qc(bundle$raw,qe$qc_config)
    write_initial_qc(initial,file.path(out,"qc","initial"))
    assert(initial$passed,"Initial QC failed; see qc/initial/Initial_QC.xlsx. No figures generated.")
    assert(identical(as.character(initial$complete$sample_id),as.character(bundle$raw$sample_id[bundle$complete])),"Prepared complete cohort disagrees with initial QC")
  }
  bundle$style$font<-cfg$font;bundle$style$dpi<-cfg$dpi
  registry<-read.csv("config/items.csv",stringsAsFactors=FALSE,na.strings="",check.names=FALSE)
  provenance<-read.csv("source_recipe_provenance.csv",stringsAsFactors=FALSE)
  analytic<-list(SINGLE_COPY_LOCI=bundle$study$single,MULTICOPY_STATE_COLUMNS=bundle$study$multicopy,
    LOCUS_DISPLAY_NAMES=bundle$study$displays[bundle$study$single],FORENSIC_MARKER_DISPLAY_NAMES=bundle$study$displays)
  af<-recipe_environment("analysis_frequency",analytic,bundle$style,bundle$study)
  hf<-recipe_environment("analysis_forensic",analytic,bundle$style,bundle$study)
  for(name in names(bundle$frequency)) {
    if(grepl("Typing",name))next
    if(grepl("SingleCopy",name))af$validate_single_copy_frequency_table(bundle$frequency[[name]])
    else af$validate_multicopy_frequency_table(bundle$frequency[[name]])
  }
  for(name in names(bundle$forensic))if(grepl("Frequencies",name))hf$validate_haplotype_frequency_table(bundle$forensic[[name]])else hf$validate_forensic_parameter_table(bundle$forensic[[name]])
  figures<-list();allqc<-list()
  if(!isFALSE(cfg$render_figures))for(i in which(registry$kind=="figure"&!is.na(registry$family)&!registry$id %in% bundle$skip_figures)) {
    item<-registry[i,];message("Original figure builder: ",item$id)
    bindings<-bundle$contexts[[item$id]]
    assert(!is.null(bindings),paste("Missing original figure input contract",item$id))
    manifest<-export_source_figure(item,bindings,bundle$style,bundle$study,cfg,file.path(out,"figures"))
    figures[[item$id]]<-manifest
    qc<-rbind(qc_row(item$id,"Color outputs","complete requested formats","created",all(cfg$formats %in% manifest$format[manifest$palette=="color"])),
      qc_row(item$id,"Grayscale outputs","complete requested formats","created",all(cfg$formats %in% manifest$format[manifest$palette=="grayscale"])),
      qc_row(item$id,"Source builder provenance","registered source recipe","registered",item$family %in% provenance$family||item$family=="pcoa"))
    if(item$id=="locus_diversity")qc<-rbind(qc,check_locus_table(bundle$tables$table_locus))
    if(item$id=="accumulation")qc<-rbind(qc,check_accumulation(bindings$plot_data))
    if(item$id=="nj")qc<-rbind(qc,check_nj(bindings$raw_tree,bindings$display_tree))
    data_bindings<-bindings[vapply(bindings,function(x)is.data.frame(x)||is.matrix(x),TRUE)]
    if(length(data_bindings)==1L)names(data_bindings)<-"Plot_Data"
    write_qc_workbook(file.path(out,"qc",paste0(item$id,"_QC.xlsx")),qc,data_bindings,
      data.frame(Source=item$family,Revision=if(is.na(item$revision_label))"locked source" else "final manuscript numbering"))
    allqc[[item$id]]<-qc
  }
  if(!isFALSE(cfg$render_tables))for(i in which(registry$kind=="table"&!is.na(registry$family))) {
    item<-registry[i,];message("Original table writer: ",item$id)
    data<-bundle$tables[[item$id]]
    custom<-item$id %in% c("table_structure","allele_workbook","recoding_audit","pairwise_statistics")
    if(!custom)assert(!is.null(data),paste("Missing table contract",item$id))
    title<-paste0(if(is.na(item$revision_label))item$locked_label else item$revision_label,". ",gsub("_"," ",sub("^table_","",item$id)))
    sheet<-paste0("Table_",sprintf("%02d",match(item$id,c("table_sampling","table_forensic","table_locus","table_panels","table_structure","table_predicted"))))
    if(item$id=="sensitivity_table")sheet<-"Sensitivity_Table"
    table_meta<-list(sheet=sheet,title=title)
    qc<-qc_row(item$id,"Prepared source contract supplied","nonempty input","supplied",if(custom)TRUE else nrow(data)>0)
    footnotes<-data.frame(Symbol="Note",Footnote=if(cfg$synthetic)"SYNTHETIC EXAMPLE: invented inputs. Not manuscript results." else "See approved study-specific caption and footnote specification.")
    source_audit<-data.frame(Source="Explicit prepared-input bundle",Synthetic=cfg$synthetic)
    e<-recipe_environment(item$family,list(),bundle$style,bundle$study,table_meta)
    path<-file.path(out,"tables",paste0(item$id,".xlsx"))
    assert(!file.exists(path),"Refusing to overwrite a table")
    if(item$id=="haplotype_catalogue") {
      e$tab_dir<-file.path(out,"tables");e$write_table(data,paste0(item$id,".csv"))
    }else if(item$id=="allele_workbook") {
      e$SINGLE_COPY_LOCI<-bundle$study$single;e$SINGLE_COPY_LABELS<-unname(bundle$study$displays[bundle$study$single])
      e$DYS385_COLUMNS<-c("dys385a","dys385b");e$DYF387S1_COLUMNS<-c("dyf387s1a","dyf387s1b","dyf387s1c")
      wb<-openxlsx::createWorkbook()
      for(pop in c(bundle$study$population_label,bundle$study$populations)) {
        z<-if(pop==bundle$study$population_label)bundle$raw else bundle$raw[bundle$raw$atoll==pop,]
        for(locus in bundle$study$single)z[[locus]][z[[locus]]=="DEL"]<-NA_character_
        object<-e$build_population_table(z)
        e$write_population_sheet(wb,if(pop==bundle$study$population_label)"Overall" else pop,pop,z,object)
      }
      openxlsx::saveWorkbook(wb,path,overwrite=FALSE)
    }else if(item$id=="population_forensic_table") {
      e$tab<-data;e$notes<-footnotes;e$write_workbook(path)
    }else if(item$id=="recoding_audit") {
      list2env(bundle$supplemental$recoding,e);e$write_workbook(path)
    }else if(item$id=="pairwise_statistics") {
      list2env(bundle$supplemental$pairwise,e);e$write_workbook(path)
    }else if(item$id=="table_structure") e$write_table_workbook(path,bundle$structure$amova,bundle$structure$pairwise,bundle$structure$mantel,footnotes,qc,source_audit)
    else if(item$id=="table_panels")e$write_table_workbook(path,data,bundle$panel_definitions,footnotes,qc,source_audit)
    else get(item$builder,e)(path,data,footnotes,qc,source_audit)
    if(!custom&&item$id!="haplotype_catalogue")write_csv(data,file.path(out,"tables",paste0(item$id,".csv")))
    write_qc_workbook(file.path(out,"qc",paste0(item$id,"_QC.xlsx")),qc,
      if(is.data.frame(data))list(Table_Data=data) else NULL,source_audit)
    allqc[[item$id]]<-qc
  }
  # Source frequency/haplotype validators are retained as executable originals.
  core_qc<-qc_row("Frequency/forensic core","Original validators completed","TRUE",TRUE,TRUE)
  allqc$core<-core_qc
  observed<-list(full.n=bundle$forensic$Overall_Forensic_Parameters$N,
    complete.n=sum(bundle$complete),full.distinct_haplotypes=bundle$forensic$Overall_Forensic_Parameters$Distinct_Haplotypes,
    full.singleton_haplotypes=bundle$forensic$Overall_Forensic_Parameters$Unique_Haplotypes,
    full.HMP=bundle$forensic$Overall_Forensic_Parameters$HMP,full.HD=bundle$forensic$Overall_Forensic_Parameters$HD,
    full.DC=bundle$forensic$Overall_Forensic_Parameters$DC,locus.system_count=nrow(bundle$tables$table_locus),
    nj.negative_raw_edges=sum(bundle$contexts$nj$raw_tree$edge.length<0))
  private_qc<-check_private_targets(observed,cfg$private_targets)
  if(is.null(cfg$private_targets))private_qc$Critical<-FALSE
  allqc$private<-private_qc
  qc<-do.call(rbind,allqc);write_csv(qc,file.path(out,"qc","Project_Audit.csv"))
  original_plan<-read.csv("config/original_audit_plan.csv",stringsAsFactors=FALSE)
  targets<-if(is.null(cfg$original_targets))NULL else read.csv(safe_path(cfg$original_targets),colClasses="character",na.strings="")
  original_audit<-run_original_audit_plan(original_plan,bundle$original_observations,targets,bundle$style,bundle$study)
  write_csv(original_audit,file.path(out,"qc","Original_Project_Checks.csv"))
  qc_full<-rbind(transform(qc,Evidence=NA_character_,Notes=NA_character_),original_audit)
  original_status<-write_original_project_audit(file.path(out,"qc","Project_Audit.xlsx"),qc_full,targets,bundle$style,bundle$study)
  if(length(figures))write_csv(do.call(rbind,figures),file.path(out,"inventory","Figure_Manifest.csv"))
  registry$Rendered<-registry$id %in% c(names(figures),names(allqc))
  registry$Final_label<-ifelse(is.na(registry$revision_label),registry$locked_label,registry$revision_label)
  registry$Publication_status<-ifelse(registry$Rendered,"PRIVATE_VALIDATION_PENDING",registry$state)
  write_csv(registry,file.path(out,"inventory","Migration_Status.csv"))
  writeLines(c("This is an output-builder migration candidate, not a publication-certified reproduction.",
    "Private targets not supplied are NOT_CONFIGURED, never PASS.",
    paste("Original project audit status:",original_status),
    "The source sampling map was pending; no replacement map is invented.",
    "The revised Figure 5 label is retained with the old Figure 4 label recorded in the registry."),file.path(out,"inventory","RELEASE_GATE.txt"))
  invisible(list(manifest=if(length(figures))do.call(rbind,figures) else NULL,qc=qc,registry=registry))
}
