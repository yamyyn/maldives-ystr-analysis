rarefaction_sharing_analysis <- function(data,markers,B,seed) {
  h<-haplotype_keys(data,markers)
  be<-method_env("bootstrap_method")
  re<-method_env("rarefaction_method",list(SEED=seed,calc_haplotype_stats=be$calc_haplotype_stats))
  result<-list()
  for(group in c("atoll","region")){
    strata<-split(seq_along(h),data[[group]]);m<-min(lengths(strata))
    runs<-lapply(seq_along(strata),function(i)re$subsample_one_group(h[strata[[i]]],m,B,i+if(group=="region")100L else 0L))
    summary<-do.call(rbind,lapply(seq_along(runs),function(i)cbind(Group=names(strata)[i],runs[[i]]$summary)))
    se<-method_env("sharing_method");hd<-data.frame(Group=data[[group]],Haplotype_ID=h)
    sharing<-se$pairwise_sharing(hd,Group,names(strata))
    result[[group]]<-list(summary=summary,iterations=lapply(runs,`[[`,"iterations"),sharing=sharing,m=m)
  }
  result
}
bootstrap_uncertainty <- function(data,markers,B,seed) {
  h<-haplotype_keys(data,markers);e<-method_env("bootstrap_method")
  scopes<-c(list(Overall=seq_along(h)),setNames(split(seq_along(h),data$region),paste0("Region:",sort(unique(data$region)))),
    setNames(split(seq_along(h),data$atoll),paste0("Population:",sort(unique(data$atoll)))))
  runs<-lapply(seq_along(scopes),function(i)e$bootstrap_group(h[scopes[[i]]],B,seed+i))
  list(summary=do.call(rbind,lapply(seq_along(runs),function(i)cbind(Scope=names(scopes)[i],runs[[i]]$summary))),
    iterations=setNames(lapply(runs,`[[`,"iterations"),names(scopes)),use="Internal uncertainty diagnostics only")
}
