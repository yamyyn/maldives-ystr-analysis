require_packages <- function() {
  packages <- c("ggplot2","dplyr","tidyr","stringr","purrr","tibble","scales",
                "patchwork","ggrepel","openxlsx","ape")
  for(p in packages) {
    if(!requireNamespace(p,quietly=TRUE))stop("Required package missing: ",p,call.=FALSE)
    suppressPackageStartupMessages(library(p,character.only=TRUE))
  }
  suppressPackageStartupMessages(library(grid))
  if(.Platform$OS.type=="windows")grDevices::windowsFonts(Arial=grDevices::windowsFont("Arial"))
}

assert <- function(ok,message)if(!isTRUE(ok))stop(message,call.=FALSE)

safe_path <- function(path,exists=TRUE) {
  assert(is.character(path)&&length(path)==1L&&!is.na(path)&&nzchar(path),"A path is required")
  path<-gsub("\\\\","/",path)
  assert(!grepl("^(/|[A-Za-z]:|~)",path)&&!any(strsplit(path,"/",fixed=TRUE)[[1]] %in% c("..",".")),"Use repository-relative paths without dot segments")
  ancestor<-path
  while(!file.exists(ancestor)&&!dir.exists(ancestor))ancestor<-dirname(ancestor)
  root<-normalizePath(".",winslash="/",mustWork=TRUE)
  actual<-normalizePath(ancestor,winslash="/",mustWork=TRUE)
  assert(actual==root||startsWith(actual,paste0(root,"/")),"Path escapes repository")
  if(exists)assert(file.exists(path),paste("Missing path",path))
  path
}

read_config <- function(path) {
  env<-new.env(parent=baseenv());sys.source(safe_path(path),env)
  cfg<-env$config
  assert(is.list(cfg)&&is.logical(cfg$synthetic)&&length(cfg$synthetic)==1L,"Invalid config")
  safe_path(cfg$input);safe_path(cfg$output,FALSE)
  assert(startsWith(cfg$output,"outputs/"),"Output must be beneath outputs/")
  assert(all(cfg$formats %in% c("pdf","svg","png","tiff")),"Unsupported export format")
  if(!is.null(cfg$private_targets))safe_path(cfg$private_targets)
  cfg
}

recipe_environment <- function(family,bindings,style,study,table_meta=NULL) {
  assert(grepl("^[a-z_]+$",family),"Unsafe recipe name")
  env<-new.env(parent=globalenv())
  list2env(bindings,envir=env)
  env$style<-style;env$study<-study;env$table_meta<-table_meta
  sys.source(safe_path(file.path("R","recipes",paste0(family,".R"))),env)
  env
}

write_csv <- function(x,path)write.csv(x,path,row.names=FALSE,na="")

export_source_figure <- function(item,bindings,style,study,cfg,out) {
  files<-list()
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(mode in cfg$palettes) {
    family<-item$family
    if(family=="pcoa")family<-if(mode=="color")"pcoa_color" else "pcoa_grayscale"
    e<-recipe_environment(family,bindings,style,study)
    for(format in cfg$formats) {
      path<-file.path(out,paste0(item$id,"_",mode,".",format))
      assert(!file.exists(path),"Refusing to overwrite a figure")
      if(family=="nj") {
        e$make_nj_plot(path,format,mode)
      }else{
        if(format=="pdf")grDevices::cairo_pdf(path,width=item$width_mm/25.4,height=item$height_mm/25.4,family=style$font,bg="white")
        else if(format=="svg")grDevices::svg(path,width=item$width_mm/25.4,height=item$height_mm/25.4,family=style$font,bg="white")
        else if(format=="png")grDevices::png(path,width=item$width_mm,height=item$height_mm,units="mm",res=cfg$dpi,type="cairo",bg="white")
        else grDevices::tiff(path,width=item$width_mm,height=item$height_mm,units="mm",res=cfg$dpi,type="cairo",compression="lzw",bg="white")
        tryCatch({
          object<-if(item$builder=="combined")
            ((e$build_admin_plot(mode)+ggplot2::labs(tag="(a)")+ggplot2::theme(plot.tag=ggplot2::element_text(family=style$font,size=9.5,face="plain"),plot.tag.position=c(.01,.99))) /
             (e$build_region_plot(mode)+ggplot2::labs(tag="(b)")+ggplot2::theme(plot.tag=ggplot2::element_text(family=style$font,size=9.5,face="plain"),plot.tag.position=c(.01,.99))))+
              patchwork::plot_layout(heights=c(1.45,.72))
          else if(item$family=="nocall")e$make_combined(mode=="color")
          else if(item$builder=="build_plot")e$build_plot()
          else get(item$builder,e)(mode)
          if(inherits(object,"grob")){grid::grid.newpage();grid::grid.draw(object)}else print(object)
        },finally=grDevices::dev.off())
      }
      files[[length(files)+1L]]<-data.frame(id=item$id,designation=item$designation,
        locked_label=item$locked_label,revision_label=item$revision_label,
        palette=mode,format=format,file=basename(path),width_mm=item$width_mm,height_mm=item$height_mm,
        dpi=if(format %in% c("pdf","svg"))NA_integer_ else cfg$dpi,
        md5=unname(tools::md5sum(path)),bytes=file.info(path)$size)
    }
  }
  do.call(rbind,files)
}
