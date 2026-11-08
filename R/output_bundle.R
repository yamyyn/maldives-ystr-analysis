make_analysis_bundle <- function(s,cfg) {
  d<-s$`02_build_datasets`;full<-d$full;study<-d$study
  pop<-study$populations;regions<-study$regions;map<-study$population_regions;n<-nrow(full)
  markers<-c(study$single,names(study$multicopy))
  f<-s$`05_forensic_statistics`$full;overall<-f$Overall_Forensic_Parameters
  admin<-f$Atoll_Forensic_Parameters
  colors<-setNames(rep(c("#CC79A7","#E69F00","#0072B2","#009E73","#D55E00","#56B4E9"),length.out=length(regions)),regions)
  grays<-setNames(gray.colors(length(regions),start=.25,end=.7),regions);shapes<-setNames(rep(c(8,21,22,24,23,25),length.out=length(regions)),regions)
  hg<-s$`12_haplogroup_analysis`;observed_hg<-hg$frequency$Parent[hg$frequency$N>0]
  hcolors<-setNames(scales::hue_pal()(length(observed_hg)),observed_hg)
  breaks<-head(cumsum(as.integer(table(factor(map[pop],levels=regions)))),-1)
  ax<-function(maximum)list(limits=c(0,maximum),breaks=pretty(c(0,maximum),4))
  spatial<-s$`10_spatial_statistics`
  style<-list(font=cfg$graphics$font,dpi=cfg$graphics$dpi,region_grays=grays,haplogroup_colors=hcolors,panel_xlim=c(.45,6.4),
    panel_axes=list(dc=ax(1.05),fraction_unique=ax(1.05),hd=ax(1.05),hmp=ax(max(.001,s$`07_marker_panels`$full$HMP)*1.15)),
    accumulation_x=ax(n),accumulation_y=ax(n),ibd_x=list(limits=range(c(0,spatial$pairs$geographic_km)),breaks=pretty(spatial$pairs$geographic_km)),
    ibd_y=list(limits=range(c(0,spatial$pairs$rst))*1.1,breaks=pretty(spatial$pairs$rst)),
    correlogram_y=list(limits=c(-1,1),breaks=seq(-1,1,.5)),nocall_breaks=pretty(c(0,nrow(d$deletion))),sharing_vlines=breaks+.5,sharing_hlines=length(pop)-breaks+.5)
  common<-list(REGION_COLORS=colors,REGION_COLOURS=colors,REGION_GRAYS=grays,REGION_SHAPES=shapes,
    REGION_ORDER=regions,REGION_LEVELS=regions,REGION_GROUPS=regions,ADMIN_ORDER=pop,HEAT_LOW="#FFFFCC",HEAT_MID="#9E9AC8",HEAT_HIGH="#54278F")
  lt<-s$`04_allele_frequencies`$locus
  plot_df<-data.frame(locus=lt$Locus,diversity=lt$PD,system_type=ifelse(gsub("\\*$","",lt$Locus) %in% unname(study$multicopy),"Multicopy pattern state","Single-copy locus"))
  plot_df<-plot_df[order(plot_df$diversity),];plot_df$locus<-factor(plot_df$locus,levels=plot_df$locus)
  dat<-data.frame(Administrative_division=factor(admin$Administrative_Division,levels=pop),Region=factor(admin$Region,levels=regions),N=admin$N,DC=admin$DC,Fraction_unique_haplotypes=admin$Unique_Profile_Fraction,HMP=admin$HMP)
  panels<-s$`07_marker_panels`$panels
  pn<-c("Minimal","Yfiler","Full","Non_RM","RM")
  # Select explicitly by both keys; never rely on row grouping offsets.
  p<-do.call(rbind,lapply(pn,function(z)s$`07_marker_panels`$full[s$`07_marker_panels`$full$Panel==z&s$`07_marker_panels`$full$Scope=="Overall",]))
  pstats<-data.frame(`Marker set`=pn,N=p$N,`Analytical systems`=p$Systems,`Distinct haplotypes`=p$Distinct_Haplotypes,`Singleton haplotypes`=p$Unique_Haplotypes,HMP=p$HMP,HD=p$HD,DC=p$DC,check.names=FALSE)
  panel_df<-data.frame(dc=p$DC,fraction_unique=p$Unique_Profile_Fraction,hd=p$HD,hmp=p$HMP,x_pos=c(1,2,3,5,6),group=c(rep("Reference / assay progression",3),rep("Yfiler Plus subsets",2)))
  ord<-s$`09c_import_arlequin`$ordination
  coords<-data.frame(Admin=rownames(ord$scores),Region=factor(map[rownames(ord$scores)],levels=regions),ord$scores)
  da_context<-function(a,region=TRUE){
    score<-a$scores;score$Group<-factor(score$Group,levels=a$groups$Group)
    score$Region<-factor(full$region[match(score$sample_id,full$sample_id)],levels=regions)
    cent<-aggregate(cbind(DA1,DA2)~Group,score,mean)
    if(!region)cent$Region<-factor(map[as.character(cent$Group)],levels=regions)
    list(score=score,centroids=cent)
  }
  dr<-da_context(s$`11_geographic_dapc`$region);dp<-da_context(s$`11_geographic_dapc`$atoll,FALSE);dh<-da_context(hg$dapc)
  rare<-s$`13_rarefaction_sharing`
  rare_plot<-function(g){z<-rare[[g]]$summary;out<-data.frame(Group=factor(z$Group,levels=if(g=="atoll")pop else regions),Mean=z$Expected_Rarefied_Distinct,Lower=z$Distinct_LCL_2.5,Upper=z$Distinct_UCL_97.5);if(g=="atoll")out$Region<-factor(map[z$Group],levels=regions);out}
  se<-recipe_environment("sharing",common,style,study)
  sharing_plot<-function(g){order<-if(g=="atoll")pop else regions;z<-se$matrix_to_long(rare[[g]]$sharing$shared_distinct,order);z<-z[z$RowIndex>z$ColIndex,];z$Row<-factor(z$Row,levels=rev(order));z$Col<-factor(z$Col,levels=order);z$Label<-as.character(z$Shared);z$LabelColor<-ifelse(z$Shared>3,"white","black");z}
  pe<-recipe_environment("pairwise",common,style,study);str<-s$`09c_import_arlequin`
  pair_plot<-function(g){order<-if(g=="atoll")pop else regions;z<-pe$make_lower_long(str[[g]]$rst,order)
    sig<-str[[g]]$long;z$bonf_significant<-vapply(seq_len(nrow(z)),function(i){row<-sig[(sig$Population_1==z$row_group[i]&sig$Population_2==z$col_group[i])|(sig$Population_2==z$row_group[i]&sig$Population_1==z$col_group[i]),];nrow(row)==1L&&row$P_Bonferroni<.05},TRUE)
    z$label<-sprintf("%.4f",z$rst);z$text_colour<-ifelse(z$rst>mean(range(z$rst)),"white","black");z}
  hfreq<-hg$frequency[hg$frequency$N>0,];hfreq<-hfreq[order(-hfreq$N,hfreq$Parent),];horder<-hfreq$Parent
  hgover<-data.frame(Final_Haplogroup=factor(hfreq$Parent,levels=horder),Percent=hfreq$Percent,Label=sprintf("%d (%.1f%%)",hfreq$N,hfreq$Percent))
  hgreg<-hg$region_frequency;hgreg<-hgreg[hgreg$Parent %in% horder,];hgreg<-data.frame(Final_Haplogroup=factor(hgreg$Parent,levels=horder),Region=factor(hgreg$Group,levels=regions),Percent=hgreg$Percent)
  prov<-s$`01_import_validate`$qc$provenance;patterns<-unique(prov$No_call_pattern[nzchar(prov$No_call_pattern)])
  pattern<-if(length(patterns))do.call(rbind,lapply(seq_along(patterns),function(i)data.frame(locus=unname(study$displays[markers]),pattern=paste("Observed pattern",i),locus_call_state=ifelse(markers %in% strsplit(patterns[i],"|",fixed=TRUE)[[1]],"No allele call","Allele observed")))) else data.frame(locus=character(),pattern=character(),locus_call_state=character())
  pattern$locus<-factor(pattern$locus,levels=unname(study$displays[markers]))
  calls<-hg$calls;no<-nzchar(prov$No_call_pattern)
  bar<-if(any(no))aggregate(list(n=rep(1L,sum(no))),list(atoll_region=paste(full$atoll[no],full$region[no]),final_haplogroup=calls$parent[no]),sum) else data.frame(atoll_region=character(),final_haplogroup=character(),n=integer())
  tip<-ord$display_tree$tip.label;tr<-map[tip]
  contexts<-list(locus_diversity=list(plot_df=plot_df),pcoa=list(plot_data=coords,LOCKED_AXIS1_PERCENT=ord$axis_percent[1],LOCKED_AXIS2_PERCENT=ord$axis_percent[2]),
    population_forensic=list(dat=dat),panel_performance=list(panel_df=panel_df,GROUP_SEPARATOR_X=4,PANEL_BREAKS=c(1,2,3,5,6),PANEL_LABELS=pn),
    accumulation=list(plot_data=s$`08_haplotype_accumulation`$summary,line_types=c(`Distinct haplotypes`="solid",`Singleton haplotypes`="22"),ribbon_fills=c(`Distinct haplotypes`="grey78",`Singleton haplotypes`="grey58")),
    rarefaction=list(admin_df=rare_plot("atoll"),region_df=rare_plot("region"),ADMIN_SUBSAMPLE_N=rare$atoll$m,REGION_SUBSAMPLE_N=rare$region$m),
    sharing=list(admin_long=sharing_plot("atoll"),region_long=sharing_plot("region")),pairwise=list(admin_df=pair_plot("atoll"),region_df=pair_plot("region"),ADMIN_BLOCK_ENDS=breaks),
    predicted_distribution=list(overall=hgover,region_freq=hgreg,FINAL_HG_ORDER=horder,BAR_COLOR="#2C7FB8",BAR_GRAY="grey55"),no_call=list(shared_pattern=pattern,bar_dat=bar),
    spatial=list(ibd_df=spatial$pairs,annotation_label=sprintf("Mantel r = %.3f\nP = %.4f",spatial$summary$R,spatial$summary$P_greater),corr_df=spatial$correlogram),
    dapc_region=list(region_df=dr$score,region_centroids=dr$centroids,REGION_DA1=s$`11_geographic_dapc`$region$axis_percent[1],REGION_DA2=s$`11_geographic_dapc`$region$axis_percent[2]),
    dapc_population=list(admin_df=dp$score,admin_centroids=dp$centroids,ADMIN_DA1=s$`11_geographic_dapc`$atoll$axis_percent[1],ADMIN_DA2=s$`11_geographic_dapc`$atoll$axis_percent[2]),
    dapc_predicted=list(hg_df=dh$score,hg_centroids=dh$centroids,HG_COLOR_PALETTE=hcolors,HG_DA1=hg$dapc$axis_percent[1],HG_DA2=hg$dapc$axis_percent[2]),
    nj=list(display_tree=ord$display_tree,raw_tree=ord$raw_tree,tip_meta=data.frame(Color=unname(colors[tr]),Gray=unname(grays[tr]),Shape=unname(shapes[tr])),COLOR_PALETTE=colors,GRAYSCALE_PALETTE=grays))
  contexts<-lapply(contexts,function(x)modifyList(common,x))
  sample<-s$`03_sampling_statistics`$population;meta<-s$`03_sampling_statistics`$metadata
  sample_table<-data.frame(`Administrative division code`=sample$Group,`Administrative division`=sample$Group,Region=map[sample$Group],
    `N, full forensic dataset`=sample$Full_n,`N, complete-profile categorical dataset`=sample$Complete_n,`Excluded structural-anomaly profiles`=sample$Full_n-sample$Complete_n,
    `Weighted latitude`=as.numeric(meta$latitude_weighted[match(sample$Group,meta$atoll)]),`Weighted longitude`=as.numeric(meta$longitude_weighted[match(sample$Group,meta$atoll)]),check.names=FALSE)
  total<-sample_table[1,,drop=FALSE];total[1,1:3]<-c("Total","Total","All");total[1,4:6]<-c(n,nrow(d$complete),n-nrow(d$complete));total[1,7:8]<-NA;sample_table<-rbind(sample_table,total)
  eligible_meta<-s$`09a_prepare_structure`$arlequin_metadata
  sample_table$`N, allele-size-eligible dataset`<-as.integer(table(factor(eligible_meta$atoll,levels=sample_table[[1]])))
  sample_table$`N, allele-size-eligible dataset`[sample_table[[1]]=="Total"]<-nrow(eligible_meta)
  table2<-data.frame(Population=study$population_label,N=overall$N,`Analytical systems`=length(markers),`Distinct haplotypes`=overall$Distinct_Haplotypes,`Singleton haplotypes`=overall$Unique_Haplotypes,HMP=overall$HMP,HD=overall$HD,DC=overall$DC,check.names=FALSE)
  mainhg<-data.frame(`Final Haplogroup category`=c(as.character(hg$frequency$Parent),"Total"),Count=c(hg$frequency$N,n),`Frequency (%)`=c(hg$frequency$Percent,100),check.names=FALSE)
  sensitivity<-s$`06_sensitivity_analysis`$stratified_summary;sensitivity<-sensitivity[sensitivity$Scope=="Overall",]
  sensitivity_table<-data.frame(Parameter=sensitivity$Metric,`Full cohort`=sensitivity$Full,`Complete cohort`=sensitivity$Complete,`Absolute change`=sensitivity$Delta,`Relative change (%)`=100*sensitivity$Delta/sensitivity$Full,Interpretation=paste("Empirical P",format(sensitivity$P_two_sided,digits=4)),check.names=FALSE)
  h<-haplotype_keys(full,markers);u<-unique(h)
  catalogue<-do.call(rbind,lapply(seq_along(u),function(i){ii<-which(h==u[i]);data.frame(Haplotype_ID=sprintf("H%04d",i),Haplotype_25System=u[i],Count=length(ii),Empirical_frequency=length(ii)/n,Parent_Predictions=paste(sort(unique(calls$parent[ii])),collapse=";"),Parent_Prediction_Conflict=length(unique(calls$parent[ii]))>1,Regions=paste(sort(unique(full$region[ii])),collapse=";"),Atolls=paste(sort(unique(full$atoll[ii])),collapse=";"))}))
  pairmatrix<-function(g,value){m<-matrix(0,nrow(str[[g]]$rst),ncol(str[[g]]$rst),dimnames=dimnames(str[[g]]$rst));z<-str[[g]]$long;for(i in seq_len(nrow(z))){m[z$Population_1[i],z$Population_2[i]]<-m[z$Population_2[i],z$Population_1[i]]<-z[[value]][i]};m}
  pwsummary<-do.call(rbind,lapply(c("atoll","region"),function(g){z<-str[[g]]$long;data.frame(Level=g,Comparisons=nrow(z),Raw_significant=sum(z$P_Raw<.05),Bonferroni_significant=sum(z$P_Bonferroni<.05),BH_significant=sum(z$P_BH<.05),Notes=str$software)}))
  amova<-str$amova$AMOVA;amova<-amova[seq_len(3),]
  amovatab<-data.frame(Component=amova$Source,`Variance component`=amova$Variance_Component,`Percentage of variation`=amova$Percentage_of_Variation,Statistic=if("Fixation_Index" %in% names(amova))amova$Fixation_Index else NA_real_,`P value`=amova$P_Raw,check.names=FALSE)
  pairbindings<-list(amova_table=str$amova$AMOVA,amova_manuscript=amovatab,atoll_pairs=str$atoll$long,atoll_rst_matrix=str$atoll$rst,atoll_p_raw_matrix=str$atoll$p,atoll_p_bonf_matrix=pairmatrix("atoll","P_Bonferroni"),atoll_p_bh_matrix=pairmatrix("atoll","P_BH"),
    region_pairs=str$region$long,region_rst_matrix=str$region$rst,region_p_raw_matrix=str$region$p,region_p_bonf_matrix=pairmatrix("region","P_Bonferroni"),region_p_bh_matrix=pairmatrix("region","P_BH"),pairwise_significance_summary=pwsummary,
    negative_rst_audit=ord$negative_rst,all_settings=data.frame(Software=str$software),permutation_documentation=data.frame(Requested=cfg$iterations$structure),source_file_index=data.frame(Source=if(str$external)unlist(cfg$arlequin_results) else "Synthetic R verification"),figure_captions=data.frame(Caption="Raw RST retained; multiplicity corrections stated separately"))
  list(style=style,study=study,contexts=contexts,skip_figures=if(!any(no))"no_call" else character(),raw=full,complete=full$sample_id %in% d$complete$sample_id,frequency=s$`04_allele_frequencies`$full,forensic=f,
    original_observations=if(is.null(cfg$original_observations))NULL else read_qc_csv(cfg$original_observations),
    supplemental=list(recoding=list(audit_summary=data.frame(Check="Posterior aggregation and leading constituent fitness",Status="COMPUTED"),mapping=setNames(hg$mapping,c("Predictor_Haplogroup","Final_Haplogroup")),final_summary=mainhg,best=list(source="Configured posterior and parent mapping")),pairwise=pairbindings),
    tables=list(table_sampling=sample_table,table_forensic=table2,table_locus=lt,table_panels=pstats,table_predicted=mainhg,population_forensic_table=dat,sensitivity_table=sensitivity_table,haplotype_catalogue=catalogue),
    panel_definitions=data.frame(`Marker set`=names(panels),`Analytical systems`=lengths(panels),`Systems included`=vapply(panels,paste,"",collapse=";"),check.names=FALSE),
    structure=list(amova=amovatab,pairwise=pwsummary,mantel=data.frame(`Genetic distance`="Raw RST",`Geographic distance`="Weighted-coordinate Haversine",`Mantel r`=spatial$summary$R,Permutations=cfg$iterations$mantel,`Directional P`=spatial$summary$P_greater,`Two-sided P`=spatial$summary$P_two_sided,check.names=FALSE)))
}
