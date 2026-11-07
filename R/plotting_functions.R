output_config <- function(cfg,figures=TRUE,tables=TRUE) {
  list(initial_qc=cfg$qc_config,skip_initial_qc=TRUE,font=cfg$graphics$font,dpi=cfg$graphics$dpi,
    formats=unlist(cfg$graphics$formats),palettes=unlist(cfg$graphics$palettes),synthetic=cfg$mode=="synthetic",
    private_targets=cfg$private_targets,original_targets=cfg$original_targets,render_figures=figures,render_tables=tables)
}
generate_analysis_tables <- function(s,cfg) {
  bundle<-make_analysis_bundle(s,cfg)
  out<-file.path(cfg$out,"publication_outputs")
  run_source_outputs(bundle,output_config(cfg,FALSE,TRUE),out)
  # Additional inference tables remain distinct from the compact manuscript table.
  openxlsx::write.xlsx(list(Stratified=s$`06_sensitivity_analysis`$stratified_summary,
    Unstratified=s$`06_sensitivity_analysis`$unstratified_summary,Paired_Bootstrap=s$`06_sensitivity_analysis`$bootstrap_summary),
    file.path(out,"tables/Sensitivity_Inference.xlsx"),creator="Y-STR workflow",overwrite=FALSE)
  cv<-rbind(transform(s$`11_geographic_dapc`$region$metadata,Scheme="Regions"),
    transform(s$`11_geographic_dapc`$atoll$metadata,Scheme="Populations"),transform(s$`12_haplogroup_analysis`$dapc$metadata,Scheme="Predicted parents"))
  write_csv(cv,file.path(out,"tables/DAPC_Assignment_Comparison.csv"))
  openxlsx::write.xlsx(list(Assignment_Comparison=cv,Region_CV=s$`11_geographic_dapc`$region$cv_summary,
    Population_CV=s$`11_geographic_dapc`$atoll$cv_summary,Parent_CV=s$`12_haplogroup_analysis`$dapc$cv_summary,
    Parent_Rare_Group_CV=s$`12_haplogroup_analysis`$rare_group_sensitivity$cv_summary),file.path(out,"tables/DAPC_Assignment_Comparison.xlsx"),creator="Y-STR workflow",overwrite=FALSE)
  list(bundle=bundle,assignment_comparison=cv,output=out)
}
export_extra_plot <- function(p,id,cfg,width=174,height=140) {
  rows<-list();out<-file.path(cfg$out,"figures/additional");dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(mode in unlist(cfg$graphics$palettes))for(fmt in unlist(cfg$graphics$formats)){
    plot<-if(is.function(p))p(mode) else p
    path<-file.path(out,paste0(id,"_",mode,".",fmt))
    assert(!file.exists(path),"Refusing to overwrite extra figure")
    if(fmt=="pdf")grDevices::cairo_pdf(path,width=width/25.4,height=height/25.4,family=cfg$graphics$font)
    else if(fmt=="svg")grDevices::svg(path,width=width/25.4,height=height/25.4,family=cfg$graphics$font,bg="white")
    else if(fmt=="png")grDevices::png(path,width=width,height=height,units="mm",res=cfg$graphics$dpi,type="cairo",bg="white")
    else grDevices::tiff(path,width=width,height=height,units="mm",res=cfg$graphics$dpi,type="cairo",compression="lzw",bg="white")
    tryCatch(print(plot),finally=dev.off())
    rows[[length(rows)+1L]]<-data.frame(ID=id,Palette=mode,Format=fmt,File=path,Width_mm=width,Height_mm=height,DPI=if(fmt %in% c("pdf","svg"))NA else cfg$graphics$dpi,MD5=unname(tools::md5sum(path)))
  }
  do.call(rbind,rows)
}
generate_analysis_figures <- function(s,cfg) {
  bundle<-s$`15_generate_tables`$bundle
  out<-file.path(cfg$out,"publication_outputs")
  source_result<-run_source_outputs(bundle,output_config(cfg,TRUE,FALSE),out)
  manifest<-list();add<-function(p,id,width=174,height=140){manifest[[id]]<<-export_extra_plot(p,id,cfg,width,height)}
  # Original final allele-frequency recipes, one combined plate per population.
  study<-s$`02_build_datasets`$study;freq<-s$`04_allele_frequencies`$full
  fe<-method_env("frequency_plots",list(FIGURE_FONT=cfg$graphics$font,LOCUS_ORDER=unname(study$displays[study$single]),
    MULTICOPY_ORDER=unname(study$multicopy),COLOUR_SINGLE="#0072B2",COLOUR_DYS385="#0072B2",COLOUR_DYF387S1="#D55E00",Y_HEADROOM_MULTIPLIER=1.08,Y_ROUND_STEP=.05))
  for(level in c("Overall","Regions","Atolls")){
    a<-freq[[paste0(level,"_SingleCopy")]];b<-freq[[paste0(level,"_Multicopy")]]
    assert(!is.null(a)&&!is.null(b),paste("Missing frequency plate",level))
    unit_column<-switch(level,Overall="Population",Regions="Region",Atolls="Administrative_Division")
    for(unit in unique(a[[unit_column]])){
      aa<-a[a[[unit_column]]==unit,];bb<-b[b[[unit_column]]==unit,]
      id<-paste0("allele_frequency_",level,"_",gsub("[^A-Za-z0-9_-]","_",unit))
      add(function(mode)fe$make_combined_plot(aa,bb,style=if(mode=="color")"colour" else "grayscale"),paste0(id,"_combined"),174,234)
      add(function(mode)fe$make_single_copy_plot(aa,descriptive_title=paste(unit,"single-copy frequencies"),publication_mode=level=="Overall",style=if(mode=="color")"colour" else "grayscale"),paste0(id,"_single"),174,220)
      add(function(mode)fe$make_multicopy_plot(bb,descriptive_title=paste(unit,"multicopy frequencies"),publication_mode=level=="Overall",style=if(mode=="color")"colour" else "grayscale"),paste0(id,"_multicopy"),174,95)
      dest<-file.path(cfg$out,"figures/source_data");dir.create(dest,showWarnings=FALSE)
      write_csv(aa,file.path(dest,paste0(id,"_single.csv")));write_csv(bb,file.path(dest,paste0(id,"_multicopy.csv")))
    }
  }
  theme<-ggplot2::theme_classic(base_family=cfg$graphics$font,base_size=10)
  contribution<-s$`07_marker_panels`$contribution;z<-contribution[contribution$Scope=="Overall",]
  add(function(mode)ggplot(z,aes(x=reorder(Marker,Delta_HD),y=Delta_HD))+geom_col(fill=if(mode=="color")"#2C7FB8" else "grey50")+coord_flip()+labs(x=NULL,y="Full-panel HD minus leave-one-system-out HD")+theme,"individual_marker_contribution")
  compare<-s$`15_generate_tables`$assignment_comparison
  add(function(mode)ggplot(compare,aes(Scheme,100*CV_mean_success))+geom_col(fill=if(mode=="color")"#2C7FB8" else "grey50")+geom_text(aes(label=sprintf("%.1f%%",100*CV_mean_success)),vjust=-.4)+scale_y_continuous(limits=c(0,110))+labs(x=NULL,y="Mean group assignment success (%)",caption="Repeated stratified holdout; schemes may use different eligible cohorts")+theme,"dapc_cv_assignment_comparison")
  for(scheme in c("region","atoll")){
    cv<-s$`11_geographic_dapc`[[scheme]]$cv_summary
    add(function(mode)ggplot(cv,aes(n.pca,Mean_Group_Assignment_Percent))+geom_line()+geom_point()+labs(x="Retained PCs",y="Mean group assignment success (%)")+theme,paste0("geographic_dapc_cv_",scheme))
  }
  hg<-s$`12_haplogroup_analysis`
  for(name in c("dapc","rare_group_sensitivity")){
    cv<-hg[[name]]$cv_summary
    add(function(mode)ggplot(cv,aes(n.pca,Mean_Group_Assignment_Percent))+geom_line()+geom_point()+labs(x="Retained PCs",y="Mean group assignment success (%)",caption="CV eligible groups only")+theme,paste0("haplogroup_",name,"_cv"))
  }
  a<-hg$rare_group_sensitivity$scores
  add(function(mode)ggplot(a,aes(DA1,DA2))+geom_point(shape=21,fill=if(mode=="color")"#2C7FB8" else "grey60")+facet_wrap(~Group)+theme,"haplogroup_rare_group_ordination")
  hf<-hg$atoll_frequency
  add(function(mode)ggplot(hf,aes(Group,Parent,fill=Percent))+geom_tile()+
    (if(mode=="color")scale_fill_gradient(low="white",high="#54278F") else scale_fill_gradient(low="white",high="grey15"))+
    labs(x="Population",y="Predicted parent",fill="Within-population %")+theme+ggplot2::theme(axis.text.x=element_text(angle=45,hjust=1)),"administrative_haplogroup_heatmap")
  hc<-hg$confidence
  add(function(mode)ggplot(hc,aes(parent,parent_probability_pct))+geom_boxplot()+coord_flip()+labs(x="Reported predicted parent",y="Aggregated parent posterior (%)")+theme,"parent_prediction_probabilities")
  add(function(mode)ggplot(hc,aes(atoll,fill=Probability_stratum))+geom_bar(position="fill")+
    (if(mode=="color")scale_fill_brewer(palette="Set2") else scale_fill_grey())+labs(x="Population",y="Profile proportion")+theme+ggplot2::theme(axis.text.x=element_text(angle=45,hjust=1)),"administrative_probability_strata")
  boot<-s$`14_bootstrap_uncertainty`$summary;boot<-boot[boot$Metric=="HD",]
  add(function(mode)ggplot(boot,aes(reorder(Scope,Observed),Observed,ymin=LCL_2.5,ymax=UCL_97.5))+geom_pointrange()+coord_flip()+labs(x=NULL,y="HD with bootstrap percentile interval",caption="Internal uncertainty diagnostic, not independent validation")+theme,"internal_bootstrap_hd",174,190)
  extra<-do.call(rbind,manifest);write_csv(extra,file.path(cfg$out,"figures/Additional_Figure_Manifest.csv"))
  captions<-data.frame(Content_ID=c(names(bundle$contexts),names(manifest)),Caption=c(paste("Source-derived",gsub("_"," ",names(bundle$contexts)),"from current-run calculations."),paste("Additional analysis",gsub("_"," ",names(manifest)),"from current-run calculations.")),Status="DRAFT_AUTHOR_REVIEW_REQUIRED")
  write_csv(captions,file.path(cfg$out,"captions/Caption_Review.csv"))
  list(source_manifest=source_result$manifest,additional_manifest=extra,captions=captions)
}
final_pipeline_qc <- function(s,cfg) {
  result<-s$`09c_import_arlequin`;full<-s$`02_build_datasets`$full;complete<-s$`02_build_datasets`$complete
  q<-rbind(qc_row("Cohorts","Unique full IDs",nrow(full),length(unique(full$sample_id)),!anyDuplicated(full$sample_id)),
    qc_row("Cohorts","Complete subset of full",TRUE,all(complete$sample_id %in% full$sample_id),all(complete$sample_id %in% full$sample_id)),
    check_locus_table(s$`04_allele_frequencies`$locus),check_accumulation(s$`08_haplotype_accumulation`$summary),
    check_nj(result$ordination$raw_tree,result$ordination$display_tree))
  original<-read.csv(file.path(cfg$out,"publication_outputs/qc/Original_Project_Checks.csv"))
  missing<-original$Status %in% c("NOT_CONFIGURED","NOT_OBSERVED")
  mapfiles<-unlist(lapply(c("color","grayscale"),function(p)paste0("sampling_map_",p,".",c("pdf","svg","png","tiff"))))
  map_ok<-isTRUE(cfg$sampling_map_approved)&&!is.null(cfg$sampling_map_directory)&&all(file.exists(file.path(safe_path(cfg$sampling_map_directory),mapfiles)))
  gate<-data.frame(Requirement=c("Private study run","External Arlequin results","Original audit complete with no FAIL","Sampling map approved","License approved","Captions approved","Private visual/numerical comparison"),
    Status=c(if(cfg$mode=="private")"PASS" else "SYNTHETIC_ONLY",if(result$external)"PASS" else "NOT_EXTERNAL",
      if(!any(missing)&&!any(original$Status=="FAIL"))"PASS" else "INCOMPLETE",if(map_ok)"PASS" else "PENDING",if(isTRUE(cfg$license_approved))"PASS" else "PENDING",if(isTRUE(cfg$captions_approved))"PASS" else "PENDING",if(isTRUE(cfg$private_comparison_approved))"PASS" else "PENDING"))
  write_csv(q,file.path(cfg$out,"qc/Final_Engineering_Checks.csv"));write_csv(gate,file.path(cfg$out,"qc/Publication_Gate.csv"))
  assert(!any(q$Status=="FAIL"),"Final engineering QC failed")
  # Reconcile actual files from both output stages rather than overwriting status.
  reg<-read.csv("config/items.csv",na.strings="",stringsAsFactors=FALSE)
  reg$Rendered<-vapply(seq_len(nrow(reg)),function(i){id<-reg$id[i];if(reg$kind[i]=="figure")file.exists(file.path(cfg$out,"publication_outputs/figures",paste0(id,"_color.pdf"))) else any(file.exists(file.path(cfg$out,"publication_outputs/tables",paste0(id,c(".xlsx",".csv")))))},TRUE)
  write_csv(reg,file.path(cfg$out,"qc/Final_Output_Inventory.csv"))
  trace<-read.csv("docs/output_traceability.csv",check.names=FALSE)
  trace$Run_mode<-cfg$mode;trace$Verification<-"MODERN_STAGE_IMPLEMENTED; PRIVATE_OUTPUT_EQUIVALENCE_REQUIRES_REVIEW"
  write_csv(trace,file.path(cfg$out,"qc/Output_Traceability_Review.csv"))
  list(checks=q,publication_gate=gate,inventory=reg,traceability=trace,publication_ready=all(gate$Status=="PASS"))
}
publication_package <- function(s,cfg) {
  ready<-s$`17_final_qc`$publication_ready
  if(!ready){
    writeLines("NOT CERTIFIED FOR PUBLICATION. Review qc/Publication_Gate.csv. No locked manuscript package was produced.",file.path(cfg$out,"PUBLICATION_PENDING.txt"))
    return(list(status=data.frame(Status="WITHHELD",Reason="Publication requirements incomplete")))
  }
  out<-file.path(cfg$out,"locked_package");assert(!dir.exists(out),"Locked package already exists")
  dir.create(out,recursive=TRUE)
  final<-file.path(cfg$out,"publication_outputs");reg<-read.csv("config/items.csv",stringsAsFactors=FALSE,na.strings="")
  rules<-method_env("inventory_rules");manifest<-list()
  for(i in seq_len(nrow(reg))){
    item<-reg[i,]
    if(item$kind=="figure"){
      src<-if(item$id=="sampling_map")safe_path(cfg$sampling_map_directory) else file.path(final,"figures")
      files<-list.files(src,pattern=paste0("^",item$id,"_(color|grayscale)\\.(pdf|svg|png|tiff)$"),full.names=TRUE)
      for(p in c("color","grayscale"))assert(rules$check_format_rule(data.frame(Extension=tools::file_ext(files[grepl(p,basename(files))])),"pdf|tiff|png","ALL"),paste("Incomplete required figure formats",item$id,p))
    }else{
      files<-list.files(file.path(final,"tables"),pattern=paste0("^",item$id,"\\.(xlsx|csv)$"),full.names=TRUE)
      assert(rules$check_format_rule(data.frame(Extension=tools::file_ext(files)),"xlsx|csv","ANY"),paste("Missing table",item$id))
    }
    dest<-file.path(out,item$kind);dir.create(dest,showWarnings=FALSE)
    assert(all(file.copy(files,dest,overwrite=FALSE)),"Package copy failed")
    manifest[[item$id]]<-data.frame(Content_ID=item$id,File=file.path(item$kind,basename(files)),MD5=unname(tools::md5sum(files)))
  }
  file.copy(file.path(cfg$out,"captions/Caption_Review.csv"),out)
  write_csv(do.call(rbind,manifest),file.path(out,"Locked_File_Manifest.csv"))
  list(status=data.frame(Status="ASSEMBLED",Path=out),manifest=do.call(rbind,manifest))
}
