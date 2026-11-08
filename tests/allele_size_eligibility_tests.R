# Invented inputs test filtering, branch alignment, panel-specific eligibility and cancellation.
toy<-data.frame(sample_id=c("T1","T2","T3","T4"),atoll=c("A","A","B","B"),region=c("N","N","S","S"),
  dys389i=c("12","12.2","13","13"),dys389ii=c("28","28.2","29","29"),dys518=c("36.2","37","38","39"))
all_loci<-c("dys389i","dys389ii","dys518")
z<-prepare_structure(toy,all_loci)
check(identical(z$metadata$sample_id,toy$sample_id))
check(identical(z$arlequin_metadata$sample_id,c("T3","T4")))
check(nrow(z$decimal)==nrow(toy)&&nrow(z$integer)==nrow(z$arlequin_metadata))
check(z$decimal$dys389iiadj[2]==16 && !z$eligibility$Eligible[2])
reduced<-prepare_structure(toy,setdiff(all_loci,"dys518"))
check(identical(reduced$arlequin_metadata$sample_id,c("T1","T3","T4")))
check(all(z$recoding$Observed==z$recoding$Exported))
check(identical(toy$dys518,c("36.2","37","38","39")))
null<-toy;null$dys518[1]<-"DEL"
check(!allele_size_eligibility(null,all_loci)$eligible[1])
bad<-toy;bad$dys518[1]<-"unexpected"
check(reject(allele_size_eligibility(bad,all_loci)))
