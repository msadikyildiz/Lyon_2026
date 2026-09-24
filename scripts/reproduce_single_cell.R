# Execute the supplied count-matrix analysis, retaining technical sample identities.
args <- commandArgs(trailingOnly=TRUE)
full_dge <- '--dge' %in% args
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
root <- dirname(dirname(normalizePath(script, winslash='/')))
base <- file.path(root,'data/single-cell')
out <- file.path(base,'generated');dir.create(out,recursive=TRUE,showWarnings=FALSE)
run <- file.path(root,'runs/single-cell');dir.create(run,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(dplyr);library(ggplot2)})
stopifnot(packageVersion('Seurat') >= '5.0.0')
# Parse the supplied Rmd instead of maintaining a second copy of its calculations.
lines <- readLines(file.path(base,'original-code/UTSW_Adam_Figures_Filtered_V3.Rmd'))
starts <- grep('^```\\{r',lines);chunks <- list()
for (i in seq_along(starts)) {
 start <- starts[i];finish <- which(seq_along(lines)>start & lines=='```')[1]
 chunks[[i]] <- paste(lines[(start+1):(finish-1)],collapse='\n')
}
run_chunk <- function(i,code=chunks[[i+1]]) {
 cat('Executing original chunk',i,'\n');flush.console()
 eval(parse(text=code),envir=.GlobalEnv)
}
run_chunk(2,sub('./data/Filtered_Features/',paste0(base,'/counts/Filtered_Features/'),chunks[[3]],fixed=TRUE))
run_chunk(3)
run_chunk(4)
run_chunk(5)
expected <- c('3a'=7827L,'3b'=10774L,'4a'=4764L,'4b'=2764L,'7a'=9619L,'7b'=13135L)
stopifnot(all(as.integer(after_filtering[names(expected)])==expected),ncol(UTSW_filtered)==48883L)
UTSW_filtered$technical_sample <- UTSW_filtered$sample
write.csv(df,file.path(out,'filtering_counts.csv'),row.names=FALSE)
run_chunk(8);run_chunk(9);run_chunk(10)
# Original RunPCA, FindClusters and RunUMAP defaults use fixed seeds (42, 0, 42).
code <- strsplit(chunks[[12]],'UTSW_filtered_UMAP <-',fixed=TRUE)[[1]][1]
run_chunk(11,code)
UTSW_filtered$cluster_published <- as.integer(as.character(UTSW_filtered$seurat_clusters))+1L
cells <- cbind(data.frame(barcode=colnames(UTSW_filtered)),UTSW_filtered@meta.data,Embeddings(UTSW_filtered,'umap'))
write.csv(cells,file.path(out,'cells_and_embeddings.csv'),row.names=FALSE)
write.csv(max_probe_df,file.path(out,'retained_probes.csv'),row.names=FALSE)
counts <- as.data.frame(table(cluster=UTSW_filtered$cluster_published,sample=UTSW_filtered$sample))
write.csv(counts,file.path(out,'cluster_sample_counts.csv'),row.names=FALSE)
probes <- max_probe_df$Feature[max_probe_df$Gene %in% c('hipA','rplJ')]
expression <- as.data.frame(t(as.matrix(GetAssayData(UTSW_filtered,assay='ProbeCollapseAssay',layer='data')[probes,,drop=FALSE])))
expression$barcode <- rownames(expression)
write.csv(expression,file.path(out,'figure6_expression.csv'),row.names=FALSE)
saveRDS(UTSW_filtered,file.path(run,'seurat.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'R-session.txt'))
cat('Saved 48,883 cells, technical sample metadata, cluster counts and UMAP coordinates.\n');flush.console()
if (full_dge) {
 DefaultAssay(UTSW_filtered) <- 'ProbeCollapseAssay'
 for (group in c('cluster_published','sample')) {
  Idents(UTSW_filtered) <- group
  markers <- FindAllMarkers(UTSW_filtered,assay='ProbeCollapseAssay',logfc.threshold=.1,min.pct=.01,return.thresh=.01)
  write.csv(markers,file.path(out,paste0(group,'_DGE_recomputed.csv')),row.names=FALSE)
 }
 for (culture in c('4','7')) {
  markers <- FindMarkers(UTSW_filtered,ident.1=culture,ident.2='parent',group.by='sample',assay='ProbeCollapseAssay',logfc.threshold=.1,min.pct=.01)
  markers$gene <- rownames(markers)
  write.csv(markers,file.path(out,paste0('sample',culture,'_vs_parent_recomputed.csv')),row.names=FALSE)
  targeted <- FindMarkers(UTSW_filtered,ident.1=culture,ident.2='parent',group.by='sample',assay='ProbeCollapseAssay',features=c('hipA-3','rplJ-2'),min.pct=0,logfc.threshold=0)
  targeted$gene <- rownames(targeted);targeted$culture <- culture
  write.csv(targeted,file.path(out,paste0('marker_claim_check_',culture,'.csv')),row.names=FALSE)
 }
}
