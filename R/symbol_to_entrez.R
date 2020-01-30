#' Obtain entrezgene IDs for gene symbols
#'
#' @param df A data.frame containing a column named either "hgnc_symbol" or
#' "mgi_symbol".
#' @param dir_annotation The directory to store annotation reference.
#' @param verbose Set to FALSE to avoid messages.
#' @return The data.frame df with an added column containing entrezgene IDs.
#' @export
symbol_to_entrez <- function(df,
                             dir_annotation = getwd(),
                             verbose = TRUE) {
  if("hgnc_symbol" %in% names(df)) {
    species <- 1 # Human species.
  } else if("mgi_symbol" %in% names(df)) {
    species <- 2 # Mouse species.
  } else {
    stop("The data frame must contain a column nammed 'hgnc_symbol' or 'mgi_symbol'.")
  }

  if(species == 1) {
    # ! TODO: use instead -> as.list(org.Hs.egSYMBOL)
    # Human
    load_file <- file.path(dir_annotation, "hgnc_to_entrez.rds")
    if(file.exists(load_file)) {
      if(verbose)
        cat("\t- loading gene info from", load_file, "\n")
      gene_info <- readRDS(load_file)
    } else {
      if(verbose)
        cat("Obtaining gene info from biomart.\n")

      if(!("mart_human" %in% ls())) init_mart_human()
      gene_info <- biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                                  mart = mart_human)
      gene_info %>%
        dplyr::group_by(hgnc_symbol) %>%
        dplyr::summarise(entrezgene_id = first(entrezgene_id)) ->
        gene_info

      gene_info$hgnc_symbol <- gsub("-", "_", gene_info$hgnc_symbol)

      save_file <- load_file
      if(verbose)
        cat("\t- saving gene info to", save_file, "\n")
      saveRDS(gene_info, save_file)
    }

    return(dplyr::left_join(df, gene_info, "hgnc_symbol"))

  } else {
    # Mouse
    load_file <- file.path(dir_annotation, "hgnc_to_entrez.rds")
    if(file.exists(load_file)) {
      if(verbose)
        cat("\t- loading gene info from", load_file, "\n")
      gene_info <- readRDS(load_file)
    } else {
      if(verbose)
        cat("Obtaining gene info from biomart.\n")

      if(!("mart_mouse" %in% ls())) init_mart_mouse()
      gene_info <- biomaRt::getBM(attributes = c("mgi_symbol", "entrezgene_id"),
                                  mart = mart_mouse)
      gene_info %>%
        dplyr::group_by(mgi_symbol) %>%
        dplyr::summarise(entrezgene_id = first(entrezgene_id)) ->
        gene_info


      gene_info$mgi_symbol <- gsub("-", "_", gene_info$mgi_symbol)

      save_file <- load_file
      if(verbose)
        cat("\t- saving gene info to", save_file, "\n")
      saveRDS(gene_info, save_file)
    }

    return(dplyr::left_join(df, gene_info, "mgi_symbol"))
  }
}
