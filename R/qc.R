qc_row <- function(section,check,expected,observed,pass,critical=TRUE) {
  data.frame(Section=section,Check=check,Expected=as.character(expected),Observed=as.character(observed),
    Status=if(is.na(pass))"NOT_CONFIGURED" else if(pass)"PASS" else "FAIL",Critical=critical)
}

check_private_targets <- function(observed,targets_path=NULL) {
  if(is.null(targets_path))return(qc_row("Private baseline","Private project targets supplied","configured private specification",
    "not configured",NA,TRUE))
  targets<-read.csv(safe_path(targets_path),colClasses="character",check.names=FALSE,na.strings="")
  assert(identical(names(targets),c("key","expected","tolerance","critical"))&&!anyDuplicated(targets$key),"Invalid private target specification")
  do.call(rbind,lapply(seq_len(nrow(targets)),function(i){
    target<-targets[i,];value<-observed[[target$key]]
    configured<-!is.na(target$expected)&&!is.null(value)
    pass<-if(!configured)NA else if(!is.na(target$tolerance))
      is.finite(as.numeric(value))&&abs(as.numeric(value)-as.numeric(target$expected))<=as.numeric(target$tolerance)
      else identical(as.character(value),target$expected)
    qc_row("Private baseline",target$key,target$expected,if(is.null(value))"not observed" else value,pass,target$critical=="TRUE")
  }))
}

check_locus_table <- function(d) {
  assert(all(c("Locus","N","Na/Np","MCAF","GD","MP","PD","Ne") %in% names(d)),"Invalid original locus-table contract")
  rbind(qc_row("Locus summary","PD defined as 1-MP","difference <= 1e-12",max(abs(d$PD-(1-d$MP))),all(abs(d$PD-(1-d$MP))<=1e-12)),
    qc_row("Locus summary","GD uses unbiased correction","difference <= 1e-12",max(abs(d$GD-d$N/(d$N-1)*(1-d$MP))),all(abs(d$GD-d$N/(d$N-1)*(1-d$MP))<=1e-12)),
    qc_row("Locus summary","Effective number defined as 1/MP","difference <= 1e-12",max(abs(d$Ne-1/d$MP)),all(abs(d$Ne-1/d$MP)<=1e-12)),
    qc_row("Locus summary","Typed denominator valid",">1",min(d$N),all(d$N>1)),
    qc_row("Locus summary","Probability bounds","0 to 1",paste(range(d$PD),collapse=" to "),all(d$PD>=0&d$PD<=1)))
}

check_accumulation <- function(d) {
  do.call(rbind,lapply(unique(d$Series),function(series){
    z<-d[d$Series==series,];last<-z[which.max(z$Sample_Size),]
    rbind(qc_row("Accumulation",paste(series,"valid percentile bounds"),"0 <= lower <= upper <= sample size",
      "checked",all(z$Lower>=0&z$Lower<=z$Upper&z$Upper<=z$Sample_Size)),
      qc_row("Accumulation",paste(series,"fixed full-sample endpoint"),"lower = mean = upper",
        paste(last$Lower,last$Mean,last$Upper,sep="; "),abs(last$Lower-last$Mean)<1e-12&&abs(last$Upper-last$Mean)<1e-12))
  }))
}

check_nj <- function(raw,display) {
  rbind(qc_row("PCoA / NJ","Display topology preserved","TRUE",isTRUE(ape::all.equal.phylo(raw,display,use.edge.length=FALSE)),
      isTRUE(ape::all.equal.phylo(raw,display,use.edge.length=FALSE))),
    qc_row("PCoA / NJ","Original tip order preserved","TRUE",identical(raw$tip.label,display$tip.label),identical(raw$tip.label,display$tip.label)),
    qc_row("PCoA / NJ","Original edge order preserved","TRUE",identical(raw$edge,display$edge),identical(raw$edge,display$edge)),
    qc_row("PCoA / NJ","Display has no negative edges",0,sum(display$edge.length<0),all(display$edge.length>=0)))
}

write_qc_workbook <- function(path,qc,plot_data=NULL,source_audit=NULL) {
  wb<-openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb,"Figure_QC");openxlsx::writeData(wb,"Figure_QC",qc)
  if(!is.null(plot_data)) {
    if(!is.list(plot_data)||is.data.frame(plot_data))plot_data<-list(Plot_Data=plot_data)
    for(name in names(plot_data)) {
      sheet<-substr(gsub("[^A-Za-z0-9_]","_",name),1,31)
      openxlsx::addWorksheet(wb,sheet);openxlsx::writeData(wb,sheet,as.data.frame(plot_data[[name]]))
    }
  }
  if(!is.null(source_audit)){openxlsx::addWorksheet(wb,"Source_Audit");openxlsx::writeData(wb,"Source_Audit",source_audit)}
  for(sheet in names(wb))openxlsx::setColWidths(wb,sheet,cols=1:10,widths="auto")
  openxlsx::saveWorkbook(wb,path,overwrite=FALSE)
}
