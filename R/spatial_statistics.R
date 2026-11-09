spatial_analysis <- function(rst,metadata,B,seed,correlogram_B=B,correlogram_seed=seed) {
  validate_pairwise(rst)
  cols<-c("latitude_weighted","longitude_weighted")
  assert(all(cols %in% names(metadata)),"Weighted coordinates required for spatial analysis")
  meta<-metadata[match(rownames(rst),metadata$atoll),,drop=FALSE]
  assert(!anyNA(meta$atoll)&&!anyDuplicated(meta$atoll),"Spatial labels do not match pairwise matrix")
  lat<-as.numeric(meta$latitude_weighted);lon<-as.numeric(meta$longitude_weighted)
  assert(all(is.finite(lat)&abs(lat)<=90&is.finite(lon)&abs(lon)<=180),"Invalid spatial coordinates")
  he<-method_env("metadata_qc_helpers");distance<-outer(seq_along(lat),seq_along(lat),Vectorize(function(i,j)he$haversine_km(lat[i],lon[i],lat[j],lon[j])))
  dimnames(distance)<-dimnames(rst);assert(all(distance[lower.tri(distance)]>0),"Duplicate weighted coordinates require reviewed resolution")
  me<-method_env("analysis_mantel");observed<-me$mantel_observed(rst,distance);draws<-me$mantel_permute(rst,distance,B,seed)
  summary<-data.frame(R=observed,P_greater=me$empirical_p_greater(observed,draws),P_two_sided=me$empirical_p_two_sided(observed,draws),Iterations=B,Seed=seed)
  ce<-method_env("analysis_correlogram");labels<-c("0-<100 km","100-<200 km","200-<300 km","300-<400 km",">=400 km")
  weights<-lapply(labels,function(l)ce$build_class_matrix(distance,l));vectors<-lapply(weights,function(x)x[lower.tri(x)])
  counts<-vapply(vectors,sum,0.0);defined<-counts>0&counts<choose(nrow(rst),2)
  coef<-rep(NA_real_,length(labels));names(coef)<-labels
  if(any(defined))coef[defined]<-ce$compute_correlogram_coefficients(rst[lower.tri(rst)],do.call(cbind,vectors[defined]))
  set.seed(correlogram_seed);permuted<-matrix(NA_real_,correlogram_B,length(labels))
  for(b in seq_len(correlogram_B)){
    ix<-sample.int(nrow(rst));z<-rst[ix,ix];if(any(defined))permuted[b,defined]<-ce$compute_correlogram_coefficients(z[lower.tri(z)],do.call(cbind,vectors[defined]))
  }
  p<-vapply(seq_along(labels),function(i)if(defined[i])ce$directional_p(coef[i],permuted[,i]) else NA_real_,0.0)
  adjusted<-rep(NA_real_,length(p));if(any(defined))adjusted[defined]<-ce$progressive_holm(p[defined])
  corr<-data.frame(distance_class=labels,mantel_r=coef,Pairs=counts,P_directional=p,P_two_sided=vapply(seq_along(labels),function(i)if(defined[i])ce$two_sided_p(coef[i],permuted[,i]) else NA_real_,0.0),
    P_progressive_holm=adjusted,P_simultaneous_holm=p.adjust(p,"holm"),Iterations=correlogram_B,Seed=correlogram_seed,
    Status=ifelse(defined,"ESTIMATED","UNDEFINED_EMPTY_OR_CONSTANT_CLASS"))
  ix<-which(lower.tri(rst),arr.ind=TRUE)
  list(summary=summary,distance=distance,permutations=draws,correlogram=corr,correlogram_permutations=permuted,
    pairs=data.frame(Population_1=rownames(rst)[ix[,1]],Population_2=colnames(rst)[ix[,2]],geographic_km=distance[ix],rst=rst[ix]))
}
