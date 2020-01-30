#' Initialize biomaRt for mmusculus
#'
#' @return None. Creates a global variable called mart_mouse.
#' @export
init_mart_mouse <- function() {
  # Obtain MGI gene symbols and other information for each Entrez ID.
  # listMarts(); listAttributes(mart) # Useful functions to obtain more info.
  # listDatasets(useMart('ensembl'));

  # Add global variable.
  mart_mouse <<- biomaRt::useMart(biomart = "ensembl",
                                  dataset = "mmusculus_gene_ensembl")
}

#' Initialize biomaRt for hsapiens
#'
#' @return None. Creates a global variable called mart_human.
#' @export
init_mart_human <- function() {
  # Obtain HGNC gene symbols and other information for each Entrez ID.
  # listMarts(); listAttributes(mart) # Useful functions to obtain more info.
  # listDatasets(useMart('ensembl'));

  # Add global variable.
  mart_human <<- biomaRt::useMart(biomart = "ensembl",
                                  dataset = "hsapiens_gene_ensembl")
}


