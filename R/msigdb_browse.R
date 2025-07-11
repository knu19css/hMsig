#' msigdb_browser
#'
#' You can load MsigDB gene sets list on R from vector.
#' @importFrom shiny validate
#' @importFrom jsonlite validate
#' @importFrom DT dataTableOutput renderDataTable
#' @import DBI
#' @import RSQLite
#' @import dplyr
#' @param species A specie, Hu or Mm.
#' @param name The gene sets' name to load.
#' @return The gene set's that you selected.
#' @export

msigdb_browse <- function(species, name) {
  if(species == "Hu") {
    geneset_path <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName="
  } else if (species == "Mm") {
    geneset_path <- "https://www.gsea-msigdb.org/gsea/msigdb/mouse/download_geneset.jsp?geneSetName="
  } else {
    stop("It's a specie that doesn't support.")
  }

  temp.list <- c()

  temp.list <- lapply(name, function(n) {
    geneset_url <- paste0(geneset_path, n, "&fileType=json")
    temp <- fromJSON(geneset_url)
    temp_gene <- c(list("gene"), temp[[1]]$geneSymbols)
    geneset <- data.frame(gene = unlist(temp_gene))
    return(geneset)
  })
}
