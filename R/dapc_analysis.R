encode_allele_states <- function(data) {
  blocks<-lapply(names(data),function(f){
    values<-as.numeric(as.character(data[[f]]));assert(all(is.finite(values)),"DAPC requires complete numeric alleles")
    levels<-sort(unique(values));block<-vapply(levels,function(v)as.numeric(values==v),numeric(length(values)))
    if(is.null(dim(block)))block<-matrix(block,ncol=1)
    colnames(block)<-paste0(f,"__",format(levels,trim=TRUE,scientific=FALSE));block
  })
  x<-do.call(cbind,blocks);assert(all(rowSums(x)==ncol(data)),"Invalid one-hot allele encoding");x
}
fit_dapc_scheme <- function(x,groups,ids,settings,seed,cv_minimum=2L) {
  assert(nrow(x)==length(groups)&&length(ids)==nrow(x)&&!anyDuplicated(ids),"Misaligned DAPC metadata")
  order_mode<-if(is.null(settings$group_order))"alphabetical" else settings$group_order
  assert(order_mode %in% c("alphabetical","frequency"),"Unknown DAPC group order")
  counts<-table(groups);levels<-if(order_mode=="frequency")names(counts)[order(-as.integer(counts),names(counts))] else sort(unique(as.character(groups)))
  group<-factor(groups,levels=levels)
  sizes<-table(group);eligible<-group %in% names(sizes[sizes>=cv_minimum])
  cv_group<-droplevels(group[eligible]);assert(nlevels(cv_group)>=2L,"Fewer than two CV-eligible groups")
  rows <- split(seq_len(nrow(x)), group, drop = TRUE)
  train_n <- ifelse(lengths(rows)<cv_minimum,lengths(rows),pmax(1L,pmin(lengths(rows)-1L,round(settings$training*lengths(rows)))))
  set.seed(seed)
  # Keep group-block order and sample singleton groups too, matching the study's RNG stream.
  training <- replicate(settings$repetitions,
    unlist(Map(function(i,n)i[sample.int(length(i),n)],rows,train_n),use.names=FALSE),simplify=FALSE)
  full_rank <- qr(scale(x,TRUE,FALSE),tol=1e-10)$rank
  requested_grid <- sort(unique(as.integer(unlist(settings$pc_grid))))
  grid <- requested_grid[requested_grid>=1L & requested_grid<full_rank]
  if(isTRUE(settings$rank_margin_candidate)||isTRUE(settings$october30_parent_grid))grid<-sort(unique(c(grid,full_rank-1L)))
  if(isTRUE(settings$minimum_discriminant_pcs)||isTRUE(settings$october30_parent_grid))grid<-grid[grid>=nlevels(group)-1L]
  assert(length(grid)>0L,"No feasible PC candidates for DAPC CV")
  train_matrices <- lapply(training,function(i)as.data.frame(x[i,,drop=FALSE]))
  pcas <- lapply(train_matrices,function(z)ade4::dudi.pca(z,center=TRUE,scale=FALSE,scannf=FALSE,nf=max(grid)))
  grid <- grid[grid<=min(vapply(pcas,function(pca)ncol(pca$li),integer(1)))]
  assert(length(grid)>0L,"No PC candidates supported by every training-subset PCA")
  raw <- do.call(rbind,lapply(seq_along(training),function(rep){
    i<-training[[rep]];j<-setdiff(seq_len(nrow(x)),i)
    do.call(rbind,lapply(grid,function(pc){
      fit<-adegenet::dapc(train_matrices[[rep]],group[i],n.pca=pc,n.da=nlevels(group)-1L,dudi=pcas[[rep]],var.contrib=FALSE)
      assigned<-stats::predict(fit,newdata=as.data.frame(x[j,,drop=FALSE]))$assign
      success<-mean(tapply(as.character(assigned)==as.character(group[j]),droplevels(group[j]),mean))
      data.frame(n.pca=pc,Replicate=rep,Group_Assignment_Success=success)
    }))
  }))
  summary <- do.call(rbind, lapply(split(raw, raw$n.pca), function(z) data.frame(
    n.pca=z$n.pca[1L], Mean_Group_Assignment_Success=mean(z$Group_Assignment_Success),
    RMSE=sqrt(mean((1-z$Group_Assignment_Success)^2)))))
  summary$Mean_Group_Assignment_Percent <- 100 * summary$Mean_Group_Assignment_Success
  best <- max(summary$n.pca[summary$RMSE == min(summary$RMSE)])
  parsed <- list(Summary=summary, Raw=raw)
  cv <- list("Training Sample IDs"=lapply(training, function(i) ids[i]),
             "Cross-Validation Results"=raw)
  train <- length(training[[1L]]) / nrow(x)
  fit<-adegenet::dapc(x,group,n.pca=best,n.da=nlevels(group)-1L,center=TRUE,scale=FALSE,var.contrib=TRUE)
  coord<-as.matrix(fit$ind.coord);if(ncol(coord)<2)coord<-cbind(coord,DA2=0)
  scores<-data.frame(sample_id=ids,Group=as.character(group),DA1=coord[,1],DA2=coord[,2])
  pred<-as.character(fit$assign);truth<-as.character(group)
  confusion<-as.data.frame(table(Truth=truth,Assigned=pred),stringsAsFactors=FALSE)
  accuracy<-data.frame(Group=levels(group),N=as.integer(table(group)),Apparent_accuracy=vapply(levels(group),function(g)mean(pred[truth==g]==g),0.0))
  eig<-as.numeric(fit$eig);axis<-100*eig/sum(eig)
  list(fit=fit,cv=cv,cv_summary=parsed$Summary,cv_raw=parsed$Raw,scores=scores,axis_percent=c(axis,0)[1:2],
    selected_pc=best,requested_pc_grid=requested_grid,effective_pc_grid=grid,group_order=order_mode,seed=seed,groups=data.frame(Group=names(sizes),N=as.integer(sizes)),cv_eligible=data.frame(sample_id=ids,Eligible=eligible),
    metadata=data.frame(Full_n=nrow(x),Full_groups=nlevels(group),CV_n=sum(eligible),CV_groups=nlevels(cv_group),
      Requested_training=settings$training,Feasible_training=train,Repetitions=settings$repetitions,Selected_PC=best,
      CV_mean_success=parsed$Summary$Mean_Group_Assignment_Success[match(best,parsed$Summary$n.pca)],
      Independent_validation=FALSE,Method="Repeated stratified holdout; not literal k-fold CV"),confusion=confusion,accuracy=accuracy)
}
geographic_dapc <- function(structure,settings,seed) {
  x<-encode_allele_states(structure$decimal);meta<-structure$metadata
  list(matrix=x,region=fit_dapc_scheme(x,meta$region,meta$sample_id,settings,seed),
    atoll=fit_dapc_scheme(x,meta$atoll,meta$sample_id,settings,seed))
}
