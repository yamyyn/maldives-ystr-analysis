dispatch_stage <- function(id,s,cfg) {
  full<-if(!is.null(s$`02_build_datasets`))s$`02_build_datasets`$full else NULL
  complete<-if(!is.null(s$`02_build_datasets`))s$`02_build_datasets`$complete else NULL
  study<-if(!is.null(s$`02_build_datasets`))s$`02_build_datasets`$study else NULL
  markers<-if(!is.null(study))c(study$single,names(study$multicopy)) else NULL
  switch(id,
    `00_setup`={
      for(d in c("tables","figures","qc","logs","intermediate","captions"))dir.create(file.path(cfg$out,d),recursive=TRUE,showWarnings=FALSE)
      yaml::write_yaml(cfg,file.path(cfg$out,"logs/effective_configuration.yml"))
      write_csv(cfg$seeds,file.path(cfg$out,"logs/effective_random_seeds.csv"))
      writeLines(paste(RNGkind(),collapse=" / "),file.path(cfg$out,"logs/RNGkind.txt"))
      writeLines(capture.output(sessionInfo()),file.path(cfg$out,"logs/sessionInfo.txt"))
      packages<-c("ggplot2","dplyr","tidyr","purrr","adegenet","ape","yaml","openxlsx","digest")
      list(environment=data.frame(Package=packages,Version=vapply(packages,function(p)as.character(packageVersion(p)),"")),mode=cfg$mode)
    },
    `01_import_validate`={
      qe<-new.env(parent=globalenv());sys.source(cfg$qc_config,qe)
      if(cfg$mode=="private")assert(!any(grepl("synthetic|DEMO|SYN[0-9]",unlist(qe$qc_config$confirmed_patterns),ignore.case=TRUE)),"Private QC confirmation contains demonstration values")
      result<-run_initial_qc(read_qc_csv(cfg$input),qe$qc_config)
      out<-file.path(cfg$out,"qc/initial");write_initial_qc(result,out,cfg$input)
      assert(result$passed,"Initial QC blocked release. Review qc/initial")
      list(qc=result,config=qe$qc_config,checks=result$checks)
    },
    `02_build_datasets`={
      q<-s$`01_import_validate`$qc;c<-s$`01_import_validate`$config
      assert(!dir.exists(cfg$clean_out),"Clean run directory exists")
      dir.create(cfg$clean_out,recursive=TRUE)
      file.copy(list.files(file.path(cfg$out,"qc/initial/04_clean_data"),full.names=TRUE),cfg$clean_out,recursive=TRUE)
      info<-file.path(cfg$out,"qc/project_info");dir.create(info)
      file.copy(list.files(file.path(cfg$out,"qc/initial/00_project_info"),full.names=TRUE),info)
      populations<-unique(as.character(q$full$atoll));regions<-unique(as.character(q$full$region))
      mapping<-setNames(vapply(populations,function(p)unique(q$full$region[q$full$atoll==p]),""),populations)
      displays<-c(setNames(toupper(c$single),c$single),setNames(c("DYS385","DYF387S1"),names(c$multicopy)))
      displays["ygatah4"]<-"Y-GATA-H4"
      study<-list(population_label=if(cfg$mode=="synthetic")"Synthetic study" else "Study",single=c$single,
        multicopy=setNames(c("DYS385","DYF387S1"),names(c$multicopy)),displays=displays,
        rapidly_mutating=c("DYS449","DYS518","DYS570","DYS576","DYS627","DYF387S1"),populations=populations,regions=regions,population_regions=mapping,
        allele_note="Typed denominators exclude no-call observations. Multicopy states describe observed peaks; no missing copy is inferred. DC in this allele workbook is 1-MP, called PD in the locus summary.")
      list(full=q$full,complete=q$complete,deletion=q$deletion,study=study,slots=c$multicopy,provenance=q$provenance)
    },
    `03_sampling_statistics`={
      meta<-full[!duplicated(full$atoll),c("atoll","region","latitude_centroid","longitude_centroid","latitude_weighted","longitude_weighted","islands_sampled","total_samples")]
      pop<-s$`01_import_validate`$qc$counts_by_population;reg<-s$`01_import_validate`$qc$counts_by_region
      list(metadata=meta,population=pop,region=reg)
    },
    `04_allele_frequencies`={
      list(full=frequency_analysis(full,study),complete=frequency_analysis(complete,study),locus=locus_summary(full,study,s$`02_build_datasets`$slots))
    },
    `05_forensic_statistics`={
      analytic<-list(SINGLE_COPY_LOCI=study$single,MULTICOPY_STATE_COLUMNS=study$multicopy,LOCUS_DISPLAY_NAMES=study$displays[study$single],FORENSIC_MARKER_DISPLAY_NAMES=study$displays)
      e<-recipe_environment("analysis_forensic",analytic,list(),study)
      list(full=e$calculate_forensic_all_levels(full,markers),complete=e$calculate_forensic_all_levels(complete,markers),
        full_summary=scoped_metrics(full,markers),complete_summary=scoped_metrics(complete,markers))
    },
    `06_sensitivity_analysis`={
      a<-random_removal_sensitivity(full,complete,markers,cfg$iterations$sensitivity,stage_seed(cfg,"sensitivity"),TRUE)
      b<-random_removal_sensitivity(full,complete,markers,cfg$iterations$sensitivity,stage_seed(cfg,"sensitivity"),FALSE)
      c<-paired_stratified_bootstrap(full,complete,markers,cfg$iterations$paired_bootstrap,stage_seed(cfg,"paired_bootstrap"))
      list(stratified=a,unstratified=b,paired_bootstrap=c,stratified_summary=a$summary,unstratified_summary=b$summary,bootstrap_summary=c$summary)
    },
    `07_marker_panels`={
      p<-read.csv(cfg$panel_file);assert(all(p$marker %in% markers)&&!anyDuplicated(p[c("panel","marker")]),"Invalid marker panel definitions")
      panels<-lapply(split(p,p$panel),function(z)z$marker[order(z$order)])
      a<-panel_and_contribution(full,panels,markers);b<-panel_and_contribution(complete,panels,markers)
      list(panels=panels,definitions=p,full=a$performance,complete=b$performance,contribution=a$contribution,complete_contribution=b$contribution)
    },
    `08_haplotype_accumulation`=accumulation_analysis(haplotype_keys(full,markers),cfg$iterations$accumulation,stage_seed(cfg,"accumulation"),sort(unique(c(1L,nrow(full),seq_len(nrow(full)%/%10L)*10L)))),
    `09a_prepare_structure`=prepare_structure(complete,study$single),
    `09b_export_arlequin`={
      out<-file.path("external/arlequin",cfg$run_id);assert(!dir.exists(out),"Arlequin export directory exists")
      maps<-export_arlequin(s$`09a_prepare_structure`,out)
      writeLines(c("Run the three .arp projects in external Arlequin using microsatellite squared size difference (RST).",
        paste("Requested permutations:",cfg$iterations$structure),"Keep original XML/text outputs; do not edit reported labels or values.",
        "Set arlequin_results in the private config before starting the run. On resume, supply those files without changing the frozen config.",
        "The synthetic mode uses a separately labelled R verification, not simulated Arlequin output."),file.path(out,"RUN_INSTRUCTIONS.txt"))
      list(maps=maps,directory=out)
    },
    `09c_import_arlequin`={
      st<-s$`09a_prepare_structure`;maps<-s$`09b_export_arlequin`$maps
      if(cfg$mode=="private"){
        assert(!is.null(cfg$arlequin_results)&&all(file.exists(unlist(cfg$arlequin_results))),"AWAITING_EXTERNAL_ARLEQUIN: export is ready; provide the three configured result files and resume")
        result<-import_arlequin(cfg$arlequin_results,maps,ncol(st$integer),cfg$iterations$structure,nrow(st$integer))
      }else{
        x<-as.matrix(st$integer);meta<-st$arlequin_metadata
        point<-amova_point(x,meta$atoll,meta$region)
        a<-data.frame(Source=c("Among regions","Among populations within regions","Within populations"),Variance_Component=point[c("Va","Vb","Vc")],Percentage_of_Variation=point[c("percent_region","percent_population","percent_within")],P_Raw=NA_real_)
        result<-list(amova=list(AMOVA=a),external=FALSE,software="Synthetic R method-of-moments verification; not Arlequin")
        for(g in c("atoll","region")){
          pair<-pairwise_rst_verification(x,meta[[g]],cfg$iterations$structure,stage_seed(cfg,"structure"))
          labs<-sort(unique(meta[[g]]));r<-p<-matrix(0,length(labs),length(labs),dimnames=list(labs,labs))
          for(i in seq_len(nrow(pair))){j<-match(pair$population1[i],labs);k<-match(pair$population2[i],labs);r[j,k]<-r[k,j]<-pair$RST[i];p[j,k]<-p[k,j]<-pair$p_greater[i]}
          result[[g]]<-list(rst=r,p=p,long=pairwise_long(r,p),settings=data.frame(Method=result$software))
        }
      }
      result$ordination<-ordination_analysis(result$atoll$rst);result
    },
    `10_spatial_statistics`=spatial_analysis(s$`09c_import_arlequin`$atoll$rst,s$`03_sampling_statistics`$metadata,cfg$iterations$mantel,stage_seed(cfg,"mantel"),cfg$iterations$correlogram,stage_seed(cfg,"correlogram","mantel")),
    `11_geographic_dapc`=geographic_dapc(s$`09a_prepare_structure`,cfg$dapc,stage_seed(cfg,"geographic_dapc")),
    `12_haplogroup_analysis`={
      settings<-cfg$haplogroup;settings$dapc<-cfg$dapc
      settings$dapc$pc_grid<-cfg$haplogroup$pc_grid
      settings$dapc$group_order<-"frequency"
      settings$dapc$rank_margin_candidate<-TRUE
      settings$dapc$minimum_discriminant_pcs<-TRUE
      haplogroup_analysis(full,read_qc_csv(cfg$prediction_leaves),read_qc_csv(cfg$parent_mapping),settings,stage_seed(cfg,"haplogroup_dapc"),
        unlist(settings$unaffected_loci),read_qc_csv(cfg$reported_parents))
    },
    `13_rarefaction_sharing`=rarefaction_sharing_analysis(full,markers,cfg$iterations$rarefaction,stage_seed(cfg,"rarefaction")),
    `14_bootstrap_uncertainty`=bootstrap_uncertainty(full,markers,cfg$iterations$bootstrap,stage_seed(cfg,"bootstrap")),
    `15_generate_tables`=generate_analysis_tables(s,cfg),
    `16_generate_figures`=generate_analysis_figures(s,cfg),
    `17_final_qc`=final_pipeline_qc(s,cfg),
    `18_publication_package`=publication_package(s,cfg),
    stop("Unknown stage ",id))
}
