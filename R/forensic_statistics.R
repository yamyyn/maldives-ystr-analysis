method_env <- function(name,bindings=list())recipe_environment(name,bindings,list(font="Arial"),list(population_label="Study"))
haplotype_keys <- function(data,markers) {
  assert(all(markers %in% names(data))&&length(markers)>0,"Unknown or empty marker panel")
  assert(!anyNA(data[markers])&&!any(as.matrix(data[markers])==""),"Unresolved forensic states")
  apply(data[markers],1,function(x)paste(paste0(markers,"=",x),collapse="|"))
}
hap_metrics <- local({
  calculate<-NULL
  function(h){if(is.null(calculate))calculate<<-method_env("sensitivity")$compute_metrics_from_haps
    c(N=length(h),calculate(h))}
})
scoped_metrics <- function(data,markers) {
  h<-haplotype_keys(data,markers)
  scopes<-c(list(Overall=seq_len(nrow(data))),setNames(split(seq_len(nrow(data)),data$region),paste0("Region:",sort(unique(data$region)))),
    setNames(split(seq_len(nrow(data)),data$atoll),paste0("Population:",sort(unique(data$atoll)))))
  do.call(rbind,lapply(names(scopes),function(g)data.frame(Scope=g,as.list(hap_metrics(h[scopes[[g]]])),check.names=FALSE)))
}
random_removal_sensitivity <- function(full,complete,markers,B,seed,stratified=TRUE) {
  assert(B>=1&&all(complete$sample_id %in% full$sample_id),"Invalid resampling cohort")
  h<-haplotype_keys(full,markers);keep<-full$sample_id %in% complete$sample_id
  scope<-c(list(Overall=seq_len(nrow(full))),setNames(split(seq_len(nrow(full)),full$region),paste0("Region:",sort(unique(full$region)))),
    setNames(split(seq_len(nrow(full)),full$atoll),paste0("Population:",sort(unique(full$atoll)))))
  metric<-c("Distinct_Haplotypes","Unique_Haplotypes","HMP","HD","DC","Unique_Profile_Fraction")
  obs<-lapply(scope,function(ii){assert(sum(keep[ii])>1,"Sensitivity scope has fewer than two complete profiles");rbind(Full=hap_metrics(h[ii])[metric],Complete=hap_metrics(h[ii[keep[ii]]])[metric])})
  draws<-lapply(scope,function(ii)matrix(NA_real_,B,length(metric),dimnames=list(NULL,metric)))
  set.seed(seed)
  strata<-split(seq_len(nrow(full)),full$atoll)
  take<-function(ii,m)if(m==0L)integer() else ii[sample.int(length(ii),m,FALSE)]
  # One joint draw per iteration is reused for all scopes, retaining pairing.
  for(b in seq_len(B)){
    removed<-if(stratified)unlist(lapply(strata,function(ii){
      candidates<-ii[keep[ii]];n_remove<-sum(!keep[ii])
      assert(length(candidates)>=n_remove,"Not enough complete profiles for the division removal plan")
      take(candidates,n_remove)
    }),use.names=FALSE) else take(which(keep),sum(!keep))
    assert(all(keep[removed]),"Only complete profiles may be randomly removed")
    selected<-setdiff(seq_along(h),removed)
    for(g in names(scope)){
      ii<-intersect(scope[[g]],selected)
      if(length(ii)>1)draws[[g]][b,]<-hap_metrics(h[ii])[metric]
    }
  }
  summary<-do.call(rbind,lapply(names(scope),function(g){
    delta<-sweep(draws[[g]],2,obs[[g]]["Full",],"-");actual<-obs[[g]]["Complete",]-obs[[g]]["Full",]
    valid<-colSums(is.finite(delta));assert(all(valid==B),"A sensitivity scope has undefined draws; revise scope design explicitly")
    lo<-(1+colSums(sweep(delta,2,actual,"<=") ))/(B+1)
    hi<-(1+colSums(sweep(delta,2,actual,">=") ))/(B+1)
    data.frame(Scope=g,Metric=metric,Full=obs[[g]]["Full",],Complete=obs[[g]]["Complete",],Delta=actual,
      Mean_null_delta=colMeans(delta),Lower=apply(delta,2,quantile,.025),Upper=apply(delta,2,quantile,.975),
      P_lower=lo,P_upper=hi,P_two_sided=pmin(1,2*pmin(lo,hi)),Iterations=B,Stratified=stratified)
  }))
  summary$P_BH_primary<-summary$P_Bonferroni_primary<-NA_real_
  primary<-summary$Scope=="Overall"&summary$Metric %in% metric
  summary$P_BH_primary[primary]<-p.adjust(summary$P_two_sided[primary],"BH")
  summary$P_Bonferroni_primary[primary]<-p.adjust(summary$P_two_sided[primary],"bonferroni")
  list(summary=summary,draws=draws,observed=obs,seed=seed)
}
paired_stratified_bootstrap <- function(full,complete,markers,B,seed) {
  h<-haplotype_keys(full,markers);keep<-full$sample_id %in% complete$sample_id
  # QC-only paired bootstrap of complete and confirmed-excluded strata.
  strata<-split(seq_along(h),interaction(full$atoll,keep,drop=TRUE));metric<-c("HMP","HD","DC")
  set.seed(seed)
  draws<-t(replicate(B,{
    ii<-unlist(lapply(strata,function(z)z[sample.int(length(z),length(z),TRUE)]),use.names=FALSE)
    hap_metrics(h[ii[keep[ii]]])[metric]-hap_metrics(h[ii])[metric]
  }))
  list(iterations=as.data.frame(draws),summary=data.frame(Metric=metric,Lower=apply(draws,2,quantile,.025),Upper=apply(draws,2,quantile,.975),Use="QC-only paired bootstrap; not primary P values"))
}
panel_and_contribution <- function(data,panels,markers) {
  performance<-do.call(rbind,lapply(names(panels),function(p){
    z<-scoped_metrics(data,panels[[p]]);cbind(Panel=p,Systems=length(panels[[p]]),z)
  }))
  contribution<-do.call(rbind,lapply(markers,function(m){
    base<-scoped_metrics(data,markers);one<-scoped_metrics(data,m);drop<-scoped_metrics(data,setdiff(markers,m))
    data.frame(Marker=m,Scope=base$Scope,Single_PD=1-one$HMP,Single_HD=one$HD,
      Delta_HD=base$HD-drop$HD,Delta_DC=base$DC-drop$DC,Delta_HMP=drop$HMP-base$HMP)
  }))
  list(performance=performance,contribution=contribution)
}
accumulation_analysis <- function(h,B,seed,grid) {
  e<-method_env("accumulation_method");set.seed(seed)
  grid<-sort(unique(c(1L,grid,length(h))));assert(all(grid>=1&grid<=length(h)),"Invalid accumulation grid")
  raw<-do.call(rbind,lapply(seq_len(B),function(b)e$evaluate_one_permutation(h,grid,b)))
  cols<-c(Distinct_Haplotypes="Distinct haplotypes",Singleton_Haplotypes="Singleton haplotypes")
  summary<-do.call(rbind,lapply(names(cols),function(m)do.call(rbind,lapply(grid,function(n){
    v<-raw[[m]][raw$Sample_Size==n];data.frame(Sample_Size=n,Series=cols[[m]],Mean=mean(v),Lower=unname(quantile(v,.025)),Upper=unname(quantile(v,.975)))
  }))))
  list(iterations=raw,summary=summary,seed=seed)
}
