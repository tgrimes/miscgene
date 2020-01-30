#' Obtain reactome pathways
#'
#' @param species One of "mouse" or "human".
#' @param file_name The file name for pathway data. If null, defaults to
#' "reactome_[species].rds".
#' @param dir_reactome The directory to store the pathway data. Defaults to
#' the current working directory.
#' @return A named list of vectors. Each vector corresponds to a Reactome pathway
#' and contains the entrezgene IDs of the genes in that pathway.
#' @export
get_reactome_pathways <- function(species = c("human", "mouse"),
                                  file_name = NULL,
                                  dir_reactome = getwd()) {
  species <- tolower(species[1])
  if(!(species %in% c("human", "mouse")))
    stop("species must be either 'human' or 'mouse'.")

  if(is.null(file_name)) {
    file_name <- file.path(dir_reactome,
                           paste0("reactome", "_", species, ".rds"))
  } else {
    # Check that ".rds" is appended to the provided file name.
    if(!grepl(".rds", file_name)) file_name <- paste0(file_name, ".rds")
    file_name <- file.path(dir_reactome, file_name)
  }

  # Convert species name; used later with reactome.db.
  species <- ifelse(species == "human", "Homo_sapiens", "Mus_musculus")

  if(file.exists(file_name)) {
    cat("\t- loading reactome pathway data from", file_name, "\n")
    reactome_to_entrez <- readRDS(file_name)
  } else {
    cat("Getting reactome pathway information for species:", species, "\n")
    # Get list of reactome pathway names
    # Map names to ID
    # Map ID to entrez ID
    # Save results as p by g matrix. (Optional; not performed.)

    # library(reactome.db)
    # Map reactome NAME to reactome ID. Subset on pathways for this species.
    reactome_to_id <- as.list(reactome.db::reactomePATHNAME2ID)
    index <- which(grepl(gsub("_", " ", species),
                         names(as.list(reactome.db::reactomePATHNAME2ID))))
    reactome_to_id <- reactome_to_id[index]

    # Map reactome ID to entrez ID.
    id_to_entrez <- as.list(reactome.db::reactomePATHID2EXTID)
    # Subset on reactome IDs obtained in first step.
    index <- which(names(id_to_entrez) %in% sapply(reactome_to_id, dplyr::first))
    id_to_entrez <- id_to_entrez[index]
    # Remove any empty pathways.
    sizes <- sapply(id_to_entrez, function(x) length(unique(x)))
    id_to_entrez <- id_to_entrez[sizes >= 1]

    # Map reactome ID back to reactome NAME.
    id_to_reactome <- as.list(reactome.db::reactomePATHID2NAME)
    # Subset on reactome IDs obtained in previous step.
    index <- which(names(id_to_reactome) %in% names(id_to_entrez))
    id_to_reactome <- id_to_reactome[index]

    # Finally, map reactome NAME to entrez ID
    reactome_to_entrez <- id_to_entrez
    names(reactome_to_entrez) <- sapply(names(id_to_entrez), function(x) {
      dplyr::first(id_to_reactome[[which(names(id_to_reactome) == x)[1]]])
    })
    names(reactome_to_entrez) <- gsub("^[ A-Za-z]*: ", "", names(reactome_to_entrez))

    # # Vector of all entrez ID involved among pathways.
    # entrez <- sort(unique(unlist(sapply(reactome_to_entrez, c))))
    # if(!is.null(subset_on_entrez)) {
    #   entrez <- entrez[entrez %in% subset_on_entrez]
    # }
    # G <- matrix(0, nrow = length(entrez), ncol = length(reactome_to_entrez))
    # colnames(G) <- names(reactome_to_entrez)
    #
    # # Each row of G corresponds to an entrez id. For each reactome pathway,
    # # determine which entrez ids are in the pathway, then set the corresponding
    # # rows in G to be 1.
    # for(i in 1:ncol(G)) {
    #   index_entrez_in_pathway <- which(entrez %in% reactome_to_entrez[[i]])
    #   if(length(index_entrez_in_pathway) > 0) {
    #     G[index_entrez_in_pathway, i] <- 1
    #   }
    # }
    # rownames(G) <- entrez
    #
    # # Remove any empty pathways.
    # G <- G[, !(apply(G, 2, sum) == 0)]

    cat("\t- saving reactome pathway data to", file_name, "\n")
    if(!dir.exists(dir_reactome)) dir.create(dir_reactome, recursive = TRUE)
    saveRDS(reactome_to_entrez, file_name)
  }

  pathway_list <- reactome_to_entrez
  return(pathway_list)
}

#' Modify a pathway list to combine overlapping pathways.
#'
#' @param pathway_list A list of pathways obtained from
#' \code{\link{get_reactome_pathways}}.
#' @param threshold A percentage between 0 and 1. If two pathways
#' overlap by more than this amount, they are combined into one pathway.
#' @return A modified list with overlapping pathways combined together.
#' @export
remove_redundant_pathways <- function(pathway_list, threshold = 0.9) {
  similar_to <- sapply(pathway_list,
                       function(p) find_similar_pathways(pathway_list, p, threshold))
  backup <- similar_to
  similar_to <- backup
  for(i in 1:length(similar_to)) {
    if(length(similar_to[[i]]) > 1) {
      index <- similar_to[[i]]
      names(pathway_list)[i] <- paste0(names(pathway_list)[i], " (See also: ",
                                      paste0(names(pathway_list)[index],
                                             collapse = ";; "), ")")

      for(j in index) {
        similar_to[[j]] <- i
      }
    }
  }
  index <- unique(unlist(similar_to))
  return(pathway_list[index])
}

#' Find pathways that overlap with a reference pathway.
#'
#' @param pathway_list A list of pathways obtained from
#' \code{\link{get_reactome_pathways}}.
#' @param pathway A reference pathway.
#' @param threshold A percentage between 0 and 1. If two pathways
#' overlap by more than this amount, then they are considered similar.
#' @return The indices to pathway_list that contain pathways similar
#' to the reference pathway.
find_similar_pathways <- function(pathway_list, pathway, threshold = 0.9) {
  similar <- sapply(1:length(pathway_list), function(i) {
    # Jaccard index >= threshold implies pathways are similar.
    (length(intersect(pathway, pathway_list[[i]])) /
       length(union(pathway, pathway_list[[i]]))) > threshold
  })
  return(which(similar))
}
