#' Obtain gene annotations
#'
#' @param df A data.frame containing a column named "entrezgene_id".
#' @param species One of "mouse" or "human".
#' @param attribute A gene attribute. Examples include: "gene_biotype",
#' "description", "chromosome_name", etc. See listAttributes() for
#' more options.
#' @param dir_annotation The directory to store annotation reference.
#' @param verbose Set to FALSE to avoid messages.
#' @return The data.frame df with an added column containing gene symbols.
#' @export
entrez_to_attribute <- function(df, species, attribute,
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
    file_name <- paste0("entrez_to_", attribute, "_hgnc.rds")
    load_file <- file.path(dir_annotation, file_name)
    if(!file.exists(load_file) && !("mart_human" %in% ls()))
      init_mart_human()
    mart <- mart_human
  } else {
    # Mouse
    file_name <- paste0("entrez_to_", attribute, "_hgnc.rds")
    load_file <- file.path(dir_annotation, file_name)
    if(!file.exists(load_file) && !("mart_human" %in% ls()))
      init_mart_mouse()
    mart <- mart_mouse
  }

  if(file.exists(load_file)) {
    if(verbose)
      cat("\t- loading gene info from", load_file, "\n")
    gene_info <- readRDS(load_file)
  } else {
    # Note: no filter is used. Information for all entrezgene IDs are collected
    # and saved for future use. This prevents the need to retrieve information
    # multiple times due to future calls.
    gene_info <- biomaRt::getBM(attributes = c("entrezgene_id", attribute),
                                mart = mart)

    # Temporarily change name for indexing.
    colnames(gene_info)[2] <- "attribute"
    gene_info %>%
      dplyr::filter(!is.na(entrezgene_id)) %>%
      dplyr::group_by(entrezgene_id) %>%
      dplyr::summarise(attribute = attribute[1]) %>%
      dplyr::mutate(attribute = gsub("( .Source.*)", "", attribute)) %>%
      dplyr::ungroup() ->
      gene_info
    colnames(gene_info)[2] <- attribute # Undo the name change.

    save_file <- load_file
    if(verbose)
      cat("\t- saving gene info to", save_file, "\n")
    if(!dir.exists(dir_annotation)) dir.create(dir_annotation, recursive = TRUE)
    saveRDS(gene_info, save_file)
  }

  if(!is.numeric(df$entrezgene_id))
    df$entrezgene_id <- as.numeric(df$entrezgene_id)

  return(dplyr::left_join(df, gene_info, by = "entrezgene_id"))
}
