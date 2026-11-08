# The categorical and allele-size branches have explicitly separate row sets.
allele_size_eligibility <- function(data, single, no_call_tokens=c("", "DEL", "?", "NULL")) {
  assert(length(single)>0L && !anyDuplicated(single) && all(single %in% names(data)), "Invalid included allele-size loci")
  raw <- as.data.frame(lapply(data[single], as.character), stringsAsFactors=FALSE)
  missing <- as.data.frame(lapply(raw,function(v)is.na(v) | v %in% no_call_tokens))
  numeric <- as.data.frame(lapply(raw,function(v)suppressWarnings(as.numeric(v))))
  invalid <- !as.matrix(missing) & (!is.finite(as.matrix(numeric)) | as.matrix(numeric)<=0)
  assert(!any(invalid), "Invalid allele-size call: resolve QC rather than silently recoding it")
  micro <- !as.matrix(missing) & is.finite(as.matrix(numeric)) &
           abs(as.matrix(numeric)-round(as.matrix(numeric)))>1e-12
  micro[is.na(micro)] <- FALSE
  eligible <- rowSums(as.matrix(missing))==0L & rowSums(micro)==0L
  list(eligible=eligible, numeric=numeric,
       audit=data.frame(Row=seq_len(nrow(data)),sample_id=data$sample_id,Eligible=eligible,
         Missing_loci=apply(as.matrix(missing),1,function(v)paste(single[v],collapse=";")),
         Microvariant_loci=apply(micro,1,function(v)paste(single[v],collapse=";"))))
}
prepare_structure <- function(complete,single) {
  eligibility <- allele_size_eligibility(complete,single)
  x <- eligibility$numeric
  assert(all(is.finite(as.matrix(x))),"Categorical structure source must contain complete numeric states")
  assert(all(c("dys389i","dys389ii") %in% names(x)),"Configure DYS389 adjustment explicitly for this assay")
  # Eligibility is assessed on typed calls BEFORE subtraction, so partial repeats cannot cancel.
  x$dys389iiadj<-x$dys389ii-x$dys389i;x$dys389ii<-NULL
  assert(all(x$dys389iiadj>0),"DYS389II adjusted values must be positive")
  integer <- x[eligibility$eligible,,drop=FALSE]
  assert(nrow(integer)>0L && all(as.matrix(integer)==round(as.matrix(integer))),"No eligible integer allele-size cohort")
  meta <- complete[c("sample_id","atoll","region")]
  list(decimal=x,integer=integer,metadata=meta,
    arlequin_metadata=meta[eligibility$eligible,,drop=FALSE],eligibility=eligibility$audit,
    cohort_summary=data.frame(Categorical_N=nrow(x),Allele_size_N=nrow(integer),Excluded_N=sum(!eligibility$eligible)),
    recoding=data.frame(Row=rep(which(eligibility$eligible),ncol(x)),Locus=rep(names(x),each=nrow(integer)),
      Observed=unlist(integer,use.names=FALSE),Exported=unlist(integer,use.names=FALSE),Action="UNCHANGED_INTEGER_CALL"))
}
validate_pairwise <- function(x,labels=NULL,p_values=FALSE) {
  assert(is.matrix(x)&&nrow(x)==ncol(x)&&nrow(x)>=3L,"Pairwise matrix must be square with >=3 groups")
  assert(identical(rownames(x),colnames(x))&&!anyDuplicated(rownames(x))&&!is.null(rownames(x)),"Invalid matrix labels")
  if(!is.null(labels))assert(setequal(rownames(x),labels),"Pairwise labels differ from cohort")
  assert(all(is.finite(x))&&max(abs(x-t(x)))<1e-10&&max(abs(diag(x)))<1e-10,"Invalid matrix symmetry/diagonal/values")
  if(p_values)assert(all(x>=0&x<=1),"P values outside [0,1]")
  invisible(TRUE)
}
pairwise_long <- function(rst,p) {
  validate_pairwise(rst);validate_pairwise(p,rownames(rst),TRUE)
  ix<-which(lower.tri(rst),arr.ind=TRUE)
  z<-data.frame(Population_1=colnames(rst)[ix[,2]],Population_2=rownames(rst)[ix[,1]],RST=rst[ix],P_Raw=p[ix])
  z$P_Bonferroni<-p.adjust(z$P_Raw,"bonferroni");z$P_BH<-p.adjust(z$P_Raw,"BH");z
}
export_arlequin <- function(structure,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  meta<-structure$arlequin_metadata;d<-cbind(meta,structure$integer)
  loci<-names(structure$integer)
  e<-method_env("arlequin_method",list(STRUCTURE_LOCI=loci,STRUCTURE_LOCI_DISPLAY=setNames(toupper(loci),loci)))
  # Original text writer, with a count-neutral header replacing case labels.
  header<-function(title,k)c("[Profile]",paste0('Title="',title,'"'),paste0("NbSamples=",k),"DataType=MICROSAT","GenotypicData=0","GameticPhase=0","RecessiveData=0","LocusSeparator=WHITESPACE","MissingData='?'","Frequency=ABS","","[Data]","[[Samples]]","")
  maps<-list()
  for(g in c("atoll","region")){
    labels<-sort(unique(meta[[g]]));tokens<-e$make_sample_name(g,labels)
    assert(!anyDuplicated(tokens),"Arlequin token collision")
    safe_ids<-e$make_safe_token(meta$sample_id);assert(!anyDuplicated(safe_ids),"Arlequin sample ID token collision")
    maps[[g]]<-data.frame(Label=labels,Arlequin_Name=tokens)
    blocks<-unlist(lapply(seq_along(labels),function(i)e$build_sample_block(d[d[[g]]==labels[i],,drop=FALSE],tokens[i])),use.names=FALSE)
    path<-file.path(out,paste0("pairwise_",g,".arp"));writeLines(c(header(paste("Pairwise",g),length(labels)),blocks),path)
    e$validate_common_arp(path,length(labels),nrow(d))
  }
  lines<-readLines(file.path(out,"pairwise_atoll.arp"))
  regions<-sort(unique(meta$region));hier<-c("[[Structure]]",'StructureName="Populations nested within regions"',paste0("NbGroups=",length(regions)))
  for(g in regions){
    pops<-unique(meta$atoll[meta$region==g]);tokens<-maps$atoll$Arlequin_Name[match(pops,maps$atoll$Label)]
    hier<-c(hier,"Group={",paste0('"',tokens,'"'),"}")
  }
  writeLines(c(lines,hier),file.path(out,"hierarchical_amova.arp"))
  e$validate_common_arp(file.path(out,"hierarchical_amova.arp"),nrow(maps$atoll),nrow(d))
  saveRDS(maps,file.path(out,"label_maps.rds"));write_csv(structure$recoding,file.path(out,"integer_call_audit.csv"))
  write_csv(structure$eligibility,file.path(out,"allele_size_eligibility_audit.csv"))
  maps
}
import_arlequin <- function(files,maps,n_loci,permutations,n_profiles=NULL) {
  assert(all(c("amova","atoll","region") %in% names(files)),"Supply all three Arlequin results")
  paths<-lapply(files,safe_path)
  e<-method_env("arlequin_parser",list(EXPECTED_LOCI=n_loci,REQUESTED_PERMUTATIONS=permutations,EXPECTED_ATOLLS=nrow(maps$atoll),ALPHA=.05))
  a<-e$parse_hierarchical_amova(paths$amova)
  if(!is.null(n_profiles)){
    expected_df<-c(nrow(maps$region)-1L,nrow(maps$atoll)-nrow(maps$region),n_profiles-nrow(maps$atoll),n_profiles-1L)
    assert(identical(as.integer(a$AMOVA$df),as.integer(expected_df)),"AMOVA degrees of freedom differ from frozen cohort/nesting")
  }
  assert(a$Reported_Permutations>=permutations,"AMOVA reported permutations below configured request")
  result<-list(amova=a,software="External Arlequin",external=TRUE)
  for(g in c("atoll","region")){
    z<-e$parse_pairwise_arlequin(paths[[g]],nrow(maps[[g]]))
    pass<-unlist(z$Settings[grepl("_Pass$",names(z$Settings))]);assert(all(pass),paste("Arlequin settings failed",g))
    assert(setequal(z$Labels$Arlequin_Name,maps[[g]]$Arlequin_Name),"External Arlequin labels do not match frozen export")
    labels<-maps[[g]]$Label[match(rownames(z$RST_Matrix),maps[[g]]$Arlequin_Name)]
    # P-value diagonal is undefined in Arlequin ('*'); zero is display-only.
    diag(z$P_Raw_Matrix)<-0
    dimnames(z$RST_Matrix)<-dimnames(z$P_Raw_Matrix)<-list(labels,labels)
    validate_pairwise(z$RST_Matrix,maps[[g]]$Label);validate_pairwise(z$P_Raw_Matrix,maps[[g]]$Label,TRUE)
    assert(length(z$Reported_Permutations)>0L&&all(z$Reported_Permutations>=permutations),"Pairwise reported permutations below configured request")
    result[[g]]<-list(rst=z$RST_Matrix,p=z$P_Raw_Matrix,long=pairwise_long(z$RST_Matrix,z$P_Raw_Matrix),settings=z$Settings,reported_permutations=z$Reported_Permutations)
  }
  assert(all(unlist(a$Settings[grepl("_Pass$",names(a$Settings))])),"Hierarchical Arlequin settings failed")
  result
}
ordination_analysis <- function(rst) {
  validate_pairwise(rst);visual<-pmax(rst,0);diag(visual)<-0
  assert(any(visual>0),"All pairwise values are nonpositive; no informative ordination")
  fit<-cmdscale(as.dist(visual),k=2,eig=TRUE,add=TRUE,x.ret=TRUE)
  e<-method_env("ordination_method");scores<-apply(fit$points,2,e$orient_axis)
  colnames(scores)<-c("PCoA1","PCoA2")
  rawtree<-ape::nj(as.dist(visual));display<-rawtree;display$edge.length<-pmax(rawtree$edge.length,0)
  eig<-fit$eig;percent<-100*pmax(eig,0)/sum(pmax(eig,0))
  list(scores=scores,axis_percent=percent[1:2],eigenvalues=eig,cailliez=fit$ac,raw_rst=rst,visual_rst=visual,
    negative_rst=data.frame(Row=which(rst<0&lower.tri(rst),arr.ind=TRUE)[,1],Value=rst[rst<0&lower.tri(rst)]),raw_tree=rawtree,display_tree=display)
}
