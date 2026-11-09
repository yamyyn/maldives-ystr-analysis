# Methodological settings and invented-data regression checks, never manuscript targets.
demo_cfg<-load_pipeline_config("config/analysis_settings.yml")
study_cfg<-yaml::read_yaml("config/private_settings.template.yml")
study_seeds<-read.csv(study_cfg$seed_file)
check(!identical(demo_cfg$seed_file,study_cfg$seed_file))
check(study_cfg$iterations$mantel==10000L&&study_cfg$iterations$correlogram==10000L)
check(stage_seed(list(seeds=study_seeds),"sensitivity")==20261002L)
check(stage_seed(list(seeds=study_seeds),"mantel")==20261005L)
check(stage_seed(list(seeds=study_seeds),"correlogram")==20261006L)
check(identical(unlist(study_cfg$haplogroup$pc_grid),c(17L,20L,25L,30L,35L,40L,50L,60L,70L,80L,90L,100L,110L,120L,130L,140L)))
check(yaml::read_yaml("CITATION.cff")$`cff-version`=="1.2.0")

# Independent grouped split generation includes singleton RNG consumption and no row sorting.
set.seed(913)
toy_x<-matrix(rnorm(17*6),17,6)
toy_groups<-c(rep("Beta",6),rep("Alpha",6),rep("Gamma",4),"Singleton")
toy_ids<-paste0("TOY_",seq_len(nrow(toy_x)))
settings<-list(training=.75,repetitions=4,pc_grid=c(1,2,3,5,9),group_order="frequency",
               rank_margin_candidate=TRUE,minimum_discriminant_pcs=TRUE)
fit<-suppressWarnings(fit_dapc_scheme(toy_x,toy_groups,toy_ids,settings,127))
levels<-names(table(toy_groups));levels<-levels[order(-as.integer(table(toy_groups)),levels)]
group<-factor(toy_groups,levels=levels);rows<-split(seq_along(toy_groups),group)
counts<-lengths(rows);train_n<-pmax(1L,pmin(counts-1L,round(.75*counts)))
set.seed(127)
independent_training<-replicate(4,unlist(Map(function(ii,n)ii[sample.int(length(ii),n)],rows,train_n),use.names=FALSE),simplify=FALSE)
check(identical(fit$cv[["Training Sample IDs"]],lapply(independent_training,function(ii)toy_ids[ii])))
check(all(vapply(fit$cv[["Training Sample IDs"]],function(ii)"TOY_17" %in% ii,TRUE)))
check(!any(fit$cv_eligible$Eligible[fit$cv_eligible$sample_id=="TOY_17"]))
check(all(fit$effective_pc_grid<qr(scale(toy_x,TRUE,FALSE),tol=1e-10)$rank))
check(all(fit$effective_pc_grid>=nlevels(group)-1L))
check(fit$selected_pc==max(fit$cv_summary$n.pca[fit$cv_summary$RMSE==min(fit$cv_summary$RMSE)]))
train<-independent_training[[1]];test<-setdiff(seq_len(nrow(toy_x)),train)
pca<-ade4::dudi.pca(as.data.frame(toy_x[train,]),center=TRUE,scale=FALSE,scannf=FALSE,nf=max(fit$effective_pc_grid))
pc<-fit$effective_pc_grid[1]
model<-adegenet::dapc(as.data.frame(toy_x[train,]),group[train],n.pca=pc,n.da=nlevels(group)-1,dudi=pca,var.contrib=FALSE)
pred<-predict(model,newdata=as.data.frame(toy_x[test,]))$assign
expected<-mean(tapply(pred==group[test],droplevels(group[test]),mean))
observed<-subset(fit$cv_raw,Replicate==1&n.pca==pc)$Group_Assignment_Success
check(abs(expected-observed)<1e-12)

# Mantel and correlogram use separately recorded seeds and iteration counts.
labels<-c("A","B","C","D")
rst<-matrix(c(0,.02,.07,-.01,.02,0,.04,.03,.07,.04,0,.08,-.01,.03,.08,0),4,dimnames=list(labels,labels))
meta<-data.frame(atoll=labels,latitude_weighted=c(0,.4,1.6,4.5),longitude_weighted=c(73,73.1,73.2,73))
sp<-spatial_analysis(rst,meta,11,123,17,456)
check(length(sp$permutations)==11L&&nrow(sp$correlogram_permutations)==17L)
check(all(sp$summary$Seed==123)&&all(sp$correlogram$Seed==456))
me<-method_env("analysis_mantel")
check(identical(sp$permutations,me$mantel_permute(rst,sp$distance,11,123)))
sp_other<-spatial_analysis(rst,meta,11,123,17,457)
check(identical(sp$permutations,sp_other$permutations))
check(!identical(sp$correlogram_permutations,sp_other$correlogram_permutations))
ce<-method_env("analysis_correlogram")
defined<-sp$correlogram$Status=="ESTIMATED"
check(isTRUE(all.equal(sp$correlogram$P_simultaneous_holm,p.adjust(sp$correlogram$P_directional,"holm"))))
check(all(is.finite(sp$correlogram$P_two_sided[defined])))
rm(toy_x,toy_groups,toy_ids,fit,sp,sp_other)
