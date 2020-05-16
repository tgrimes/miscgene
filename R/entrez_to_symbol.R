#' Obtain gene symbols for entrezgene IDs
#'
#' @param df A data.frame containing a column named "entrezgene_id".
#' @param species One of "mouse" or "human".
#' @param dir_annotation The directory to store annotation reference.
#' @param verbose Set to FALSE to avoid messages.
#' @return The data.frame df with an added column containing gene symbols.
#' @export
entrez_to_symbol <- function(df, species,
                             dir_annotation = getwd(),
                             verbose = TRUE) {
  if(!any(names(df) %in% c("entrezgene_id")))
    stop("df must contain a column nammed 'entrezgene_id'.")
  species <- pmatch(tolower(species), c("human", "mouse"))
  if(is.na(species))
    stop("species must be either 'human' or 'mouse'.")

  if(species == 1) {
    # ! TODO: use instead -> as.list(org.Hs.egSYMBOL)
    # Human
    load_file <- file.path(dir_annotation, "entrez_to_hgnc.rds")
    if(file.exists(load_file)) {
      if(verbose)
        cat("\t- loading gene info from", load_file, "\n")
      gene_info <- readRDS(load_file)
    } else {
      if(!("mart_human" %in% ls())) init_mart_human()
      gene_info <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                                  mart = mart_human)
      gene_info %>%
        dplyr::group_by(entrezgene_id) %>%
        dplyr::summarise(hgnc_symbol = first(hgnc_symbol)) ->
        gene_info

      gene_info$hgnc_symbol <- gsub("-", "_", gene_info$hgnc_symbol)

      save_file <- load_file
      if(verbose)
        cat("\t- saving gene info to", save_file, "\n")
      if(!dir.exists(dir_annotation)) dir.create(dir_annotation, recursive = TRUE)
      saveRDS(gene_info, save_file)
    }

  } else {
    # Mouse
    load_file <- file.path(dir_annotation, "entrez_to_mgi.rds")
    if(file.exists(load_file)) {
      if(verbose)
        cat("\t- loading gene info from", load_file, "\n")
      gene_info <- readRDS(load_file)
    } else {
      if(verbose)
        cat("Obtaining gene info from biomart.\n")

      if(!("mart_mouse" %in% ls())) init_mart_mouse()
      gene_info <- biomaRt::getBM(attributes = c("entrezgene_id", "mgi_symbol"),
                                  mart = mart_mouse)
      gene_info %>%
        dplyr::group_by(entrezgene_id) %>%
        dplyr::summarise(mgi_symbol = first(mgi_symbol)) ->
        gene_info

      gene_info$mgi_symbol <- gsub("-", "_", gene_info$mgi_symbol)

      save_file <- load_file
      if(verbose)
        cat("\t- saving gene info to", save_file, "\n")
      if(!dir.exists(dir_annotation)) dir.create(dir_annotation, recursive = TRUE)
      saveRDS(gene_info, save_file)
    }
  }

  if(!is.numeric(df$entrezgene_id)) {
    if(is.factor(df$entrezgene_id)) {
      df$entrezgene_id <- as.numeric(as.character(df$entrezgene_id))
    } else {
      df$entrezgene_id <- as.numeric(df$entrezgene_id)
    }
  }

  return(dplyr::left_join(df, gene_info, by = "entrezgene_id"))
}
