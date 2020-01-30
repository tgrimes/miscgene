#' Obtain gene ontology (GO) annotations
#'
#' @param df A data.frame containing a column named "entrezgene_id".
#' @param species One of "mouse" or "human".
#' @return A data.frame containing GO annotations for the genes in df.
#' @export
entrez_to_GO <- function(df, species) {
  if(!any(names(df) %in% c("entrezgene_id")))
    stop("df must contain a column nammed 'entrezgene_id'.")
  species <- pmatch(tolower(species), c("human", "mouse"))
  if(is.na(species))
    stop("species must be either 'human' or 'mouse'.")
  if(species == 1) {
    # ! TODO: use instead -> as.list(org.Hs.egSYMBOL)
    # Human
    if(!file.exists(load_file) && !("mart_human" %in% ls()))
      init_mart_human()
    mart <- mart_human
  } else {
    # Mouse
    if(!file.exists(load_file) && !("mart_human" %in% ls()))
      init_mart_mouse()
    mart <- mart_mouse
  }

  # The arguments "filters" and "values" are used to avoid loading the
  # GO features for all genes (which would be a lot of data).
  # Note that these features are not saved for future reference like
  # in symbol_to_entrez(), because not all genes are loaded.
  gene_info <- biomaRt::getBM(attributes = c("entrezgene_id", "go_id", "name_1006",
                                             "definition_1006", "go_linkage_type",
                                             "namespace_1003"),
                              filters = "entrezgene_id", values = df$entrezgene_id,
                              mart = mart)

  gene_info <- gene_info %>%
    dplyr::filter(!is.na(entrezgene_id),
           go_linkage_type %in% c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")) %>%
    dplyr::group_by(entrezgene_id, go_id) %>%
    dplyr::slice(1)

  # library(GO.db)
  # as.list(GOMFPARENTS)
  # as.list(GOMFCHILDREN)
  # as.list(GOMFOFFSPRING)["GO:0004519"]

  return(gene_info)
}
