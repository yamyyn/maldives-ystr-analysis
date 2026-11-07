aggregate_parent_predictions <- function(leaves,mapping,ids,probability_threshold=95,fitness_threshold=40,tolerance=.21,reported_parents=NULL) {
  assert(all(c("sample_id","leaf","probability_pct","fitness") %in% names(leaves)),"Prediction leaf schema invalid")
  assert(all(c("leaf","parent") %in% names(mapping))&&!anyDuplicated(mapping$leaf),"Invalid parent mapping")
  assert(!anyDuplicated(leaves[c("sample_id","leaf")])&&setequal(unique(leaves$sample_id),ids),"Prediction IDs/leaf rows invalid")
  leaves$probability_pct<-as.numeric(leaves$probability_pct);leaves$fitness<-as.numeric(leaves$fitness)
  assert(all(is.finite(leaves$probability_pct)&leaves$probability_pct>=0&leaves$probability_pct<=100)&all(is.finite(leaves$fitness)),"Invalid prediction values")
  assert(all(leaves$leaf %in% mapping$leaf),"Unmapped predictor leaf")
  sums<-tapply(leaves$probability_pct,leaves$sample_id,sum)
  assert(all(abs(sums-100)<=tolerance),"Leaf posteriors do not sum to 100 within configured rounding tolerance")
  assert(all(table(leaves$sample_id)==nrow(mapping)),"Every modeled leaf must be supplied per sample, including zero probabilities")
  leaves$parent<-mapping$parent[match(leaves$leaf,mapping$leaf)]
  allparents<-unique(mapping$parent)
  if(!is.null(reported_parents))assert(all(c("sample_id","parent") %in% names(reported_parents))&&!anyDuplicated(reported_parents$sample_id)&&setequal(reported_parents$sample_id,ids)&&all(reported_parents$parent %in% allparents),"Invalid independently reported parent labels")
  agg<-aggregate(leaves$probability_pct,leaves[c("sample_id","parent")],sum);names(agg)[3]<-"probability_pct"
  assert(nrow(agg)==length(ids)*length(allparents),"Missing modeled parent classes")
  calls<-do.call(rbind,lapply(ids,function(id){
    a<-agg[agg$sample_id==id,];a<-a[order(-a$probability_pct,a$parent),]
    parent<-if(is.null(reported_parents))a$parent[1] else reported_parents$parent[match(id,reported_parents$sample_id)]
    parent_probability<-a$probability_pct[match(parent,a$parent)]
    members<-leaves[leaves$sample_id==id&leaves$parent==parent,]
    members<-members[order(-members$probability_pct,members$leaf),]
    lead<-members[1,]
    data.frame(sample_id=id,parent=parent,parent_probability_pct=parent_probability,
      leading_leaf=lead$leaf,leading_leaf_probability_pct=lead$probability_pct,leading_leaf_fitness=lead$fitness,
      best_constituent_fitness=max(members$fitness),posterior_tie=sum(a$probability_pct==a$probability_pct[1])>1L,
      threshold_qualified=parent_probability>=probability_threshold&lead$fitness>=fitness_threshold)
  }))
  freq<-as.data.frame(table(factor(calls$parent,levels=allparents)),stringsAsFactors=FALSE);names(freq)<-c("Parent","N")
  freq$Percent<-100*freq$N/length(ids)
  list(leaves=leaves,all_parent_posteriors=agg,calls=calls,frequency=freq,mapping=mapping,
    terminology="Y-STR predicted parents; not SNP-confirmed haplogroups. Fitness is never summed or averaged.")
}
haplogroup_analysis <- function(full,leaves,mapping,settings,seed,unaffected_loci,reported_parents=NULL) {
  result<-aggregate_parent_predictions(leaves,mapping,full$sample_id,settings$probability_threshold,settings$fitness_threshold,settings$rounding_tolerance,reported_parents)
  call<-result$calls[match(full$sample_id,result$calls$sample_id),];call$atoll<-full$atoll;call$region<-full$region
  for(g in c("region","atoll")){
    tab<-as.data.frame(table(Group=full[[g]],Parent=factor(call$parent,levels=result$frequency$Parent)),stringsAsFactors=FALSE)
    tab$Percent<-100*tab$Freq/as.numeric(table(full[[g]])[tab$Group]);result[[paste0(g,"_frequency")]]<-tab
  }
  result$confidence<-transform(call,Probability_stratum=ifelse(parent_probability_pct>=settings$probability_threshold,"Above threshold","Below threshold"))
  raw<-full
  raw$dys389iiadj<-as.numeric(full$dys389ii)-as.numeric(full$dys389i)
  assert(all(unaffected_loci %in% names(raw)),"Unknown unaffected haplogroup DAPC locus")
  x<-encode_allele_states(raw[unaffected_loci]);result$matrix<-x
  # October 30: full exploratory fit; CV only groups with the configured minimum.
  result$dapc<-fit_dapc_scheme(x,call$parent,full$sample_id,settings$dapc,seed,settings$cv_minimum)
  counts<-table(call$parent)
  rare_min<-if(is.null(settings$rare_group_minimum))5L else settings$rare_group_minimum
  keep<-as.integer(counts[call$parent])>=rare_min
  result$rare_group_sensitivity<-fit_dapc_scheme(x[keep,,drop=FALSE],call$parent[keep],full$sample_id[keep],settings$dapc,seed,rare_min)
  result$calls<-call;result
}
