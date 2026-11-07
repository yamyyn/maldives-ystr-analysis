# Run from the extracted repository root. This checks source integrity only.
if(!requireNamespace("digest",quietly=TRUE))stop("Restore the documented R dependencies first")
manifest<-read.csv("SHA256_MANIFEST.csv",stringsAsFactors=FALSE)
stopifnot(all(c("File","SHA256") %in% names(manifest)),!anyDuplicated(manifest$File))
root<-normalizePath(".",winslash="/",mustWork=TRUE)
paths<-normalizePath(manifest$File,winslash="/",mustWork=TRUE)
stopifnot(all(startsWith(paths,paste0(root,"/"))))
actual<-vapply(paths,function(p){
  # Read bytes directly: some Windows hidden files report an incorrect
  # file.access() status despite being readable through a connection.
  bytes<-readBin(p,what="raw",n=file.info(p)$size)
  digest::digest(bytes,algo="sha256",serialize=FALSE)
},character(1))
if(any(actual!=manifest$SHA256))stop("Source checksum mismatch: ",paste(manifest$File[actual!=manifest$SHA256],collapse=", "))
message(length(paths)," source-file checksums verified. This is not private-data reproduction verification.")
