setClass(
  "gGSEAResult",
  slots = list(
    result   = "list",
    plot     = "list",
    grp_name = "character"
  )
)

gGSEA <- function(DEGlist, grplist) {
  result <- list()
  plot_ls <- list()
  
  for (i in seq_along(DEGlist)) {
    df <- data.frame(gene = row.names(DEGlist[[i]]),
                     logFC = DEGlist[[i]]$avg_log2FC)
    
    gene <- as.vector(df$logFC)
    names(gene) <- df$gene
    gene <- sort(gene, decreasing = T)
    
    edo2 <- tryCatch(
      GSEA(gene,
           exponent = 0,
           minGSSize = 1,
           maxGSSize = 1000,
           pvalueCutoff  = 0.05,
           pAdjustMethod = "BH",
           TERM2GENE = grplist,
           nPermSimple = 10000,
           by = "fgsea",
           eps = 0), error = function(e) NULL)
    
    if (!is.null(edo2) && nrow(edo2@result) > 0) {
      p <- dotplot(edo2, showCategory=50, font.size = 10, label_format = 60) + ggtitle(paste0("S", i-1))
      result[[i]] <- edo2@result
      plot_ls[[i]] <- p}
    p <- NULL
    edo2 <- NULL
  }
  
  result <- lapply(result, function(x) {
    if (is.null(x)) {
      return(data.frame())
    } else {
      return(x)
    }
  })
  grp_name <- deparse(substitute(grplist))
  write.xlsx(result, file = paste0("MK_Intergrated_data/Merged_", grp_name, ".xlsx", asTable = T, rowNames = T))
  new(
    "gGSEAResult",
    result   = result,
    plot     = plot_ls,
    grp_name = grp_name
  )
}
