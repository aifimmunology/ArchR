##########################################################################################
# S4Vectors/BiocGenerics Within Methods
##########################################################################################

#' Negated Value Matching
#'
#' This function is the reciprocal of %in%. See the match funciton in base R.
#'
#' @param x The value to search for in `table`.
#' @param table The set of values to serve as the base for the match function.
#' @export
"%ni%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)

#' Generic matching function for S4Vector objects
#'
#' This function provides a generic matching function for S4Vector objects primarily to avoid ambiguity.
#'
#' @param x An `S4Vector` object to search for in `table`.
#' @param table The set of `S4Vector` objects to serve as the base for the match function.
#' @export
'%bcin%' <- function(x, table) S4Vectors::match(x, table, nomatch = 0) > 0

#' Negated matching function for S4Vector objects
#'
#' This function provides the reciprocal of %bcin% for S4Vector objects primarily to avoid ambiguity.
#'
#' @param x An `S4Vector` object to search for in `table`.
#' @param table The set of `S4Vector` objects to serve as the base for the match function.
#' @export
'%bcni%' <- function(x, table) !(S4Vectors::match(x, table, nomatch = 0) > 0)

##########################################################################################
# Helper to try to reformat fragment files appropriately if a bug is found
##########################################################################################

#' Reformat Fragment Files to be Tabix and Chr Sorted
#'
#' This function provides help in reformatting Fragment Files for reading in createArrowFiles.
#' It will handle weird anomalies found that cause errors in reading tabix bgzip'd fragment files.
#'
#' @param fragmentFiles A character vector the paths to fragment files to be reformatted
#' @param checkChrPrefix A boolean value that determines whether seqnames should be checked to contain
#' "chr". IF set to `TRUE`, any seqnames that do not contain "chr" will be removed from the fragment files.
#' @export
reformatFragmentFiles <- function(
  fragmentFiles = NULL,
  checkChrPrefix = getArchRChrPrefix()
  ){

  .validInput(input = fragmentFiles, name = "fragmentFiles", valid = c("character"))
  .validInput(input = checkChrPrefix, name = "checkChrPrefix", valid = c("boolean"))

  options(scipen = 999)
  .requirePackage("data.table")
  .requirePackage("Rsamtools")
  for(i in seq_along(fragmentFiles)){
    message(i, " of ", length(fragmentFiles))
    dt <- data.table::fread(fragmentFiles[i])
    dt <- dt[order(dt$V1,dt$V2,dt$V3), ]
    if(checkChrPrefix){
      idxRemove1 <- which(substr(dt$V1,1,3) != "chr")
    }else{
      idxRemove1 <- c()
    }
    idxRemove2 <- which(dt$V2 != as.integer(dt$V2))
    idxRemove3 <- which(dt$V3 != as.integer(dt$V3))
    #get all
    idxRemove <- unique(c(idxRemove1, idxRemove2, idxRemove3))
    if(length(idxRemove) > 0){
      dt <- dt[-idxRemove,]
    }
    if(nrow(dt) == 0){
      if(checkChrPrefix){
        stop("No fragments found after checking for integers and chrPrefix!")
      }else{
        stop("No fragments found after checking for integers!")
      }
    }
    #Make sure no spaces or #
    dt$V4 <- gsub(" |#", ".", dt$V4)
    fileNew <- gsub(".tsv.bgz|.tsv.gz", "-Reformat.tsv", fragmentFiles[i])
    data.table::fwrite(dt, fileNew, sep = "\t", col.names = FALSE)
    Rsamtools::bgzip(fileNew)
    file.remove(fileNew)
    .fileRename(paste0(fileNew, ".bgz"), paste0(fileNew, ".gz"))
  }
}


##########################################################################################
# Helper For cleaning up coverage file paths
##########################################################################################

#' Iterate over existing coverage file paths and files 
#' to clean up and remove any that are not found, and add any new ones.
#' @param proj An ArchRProject object
#' @param removeOld A boolean value that determines whether old files should be removed if not found.
#' @export
correctGroupCoveragePaths <- function(proj, removeOld = TRUE) {
  # Get the main directory
  mainDir <- getOutputDirectory(proj)
  #File all coverage files
  allH5s <- list.files(file.path(mainDir, "GroupCoverages"), pattern = ".h5$", full.names = TRUE, recursive = TRUE)
  unMatchedFiles = list()
  for (i in seq_along(proj@projectMetadata$GroupCoverages)) {

    # Get the coverage metadata for that call of group coverages
    coverageMetadata <- proj@projectMetadata$GroupCoverages[[i]]$coverageMetadata

    # Get the folder name of the files
    coverageFolder = names(proj@projectMetadata$GroupCoverages)[i]

    if(!any(grepl(coverageFolder, allH5s)) & !removeOld){
      warning(paste("No files found for group coverage:", coverageFolder))
      next
    }else if(!any(grepl(coverageFolder, allH5s)) & removeOld){
      warning(paste("Removing paths for group coverage not found:", coverageFolder))
      proj@projectMetadata$GroupCoverages[[i]] <- NULL
      if (dir.exists(file.path(mainDir, "GroupCoverages", coverageFolder))) {
        unlink(file.path(mainDir, "GroupCoverages", coverageFolder), recursive = TRUE)
      }
      next
    }

    # Get the base names of the files
    newBaseNames <- basename(allH5s)
    # Generate the folder + file names for alignment
    newFileNames = paste(coverageFolder, newBaseNames, sep = "/")
    # Filter these filenames to only include ones that actually exist. This will enable matching of both cell type name and group coverage folders at the same time
    # By matching folder/groupCoverage.h5 together, we will control for the possibility that the group coverage folder name matches with a cell type coverage name. 
    newFilePaths= unlist(lapply(newFileNames, function(x) { grep(x, allH5s, value = TRUE) })) 
    # Generate the folder + filenames for old files
    oldBaseNames <- basename(coverageMetadata$File)
    oldFileNames <- paste(coverageFolder, oldBaseNames, sep = "/")
    matchedIndex <- match(oldFileNames, newFileNames)
    # Get the corrected paths
    correctedPaths <- newFilePaths[matchedIndex]
    # Get the unmatched files
    unMatchedFiles <- c(unMatchedFiles, newFilePaths[is.na(matchedIndex)])
    # Remove any old files that are not found

    
    if(any(oldFileNames %in% newFileNames) & removeOld){
      warning(paste("Removing paths for files not found:", baseName))
      coverageMetadata = coverageMetadata[oldFileNames %in% newFileNames,]
    }

    proj@projectMetadata$GroupCoverages[[i]]$coverageMetadata$File <- unlist(correctedPaths)
  }

  # Remove any unmatched files
  if(removeOld){
    for(i in seq_along(unMatchedFiles)){
      warning(paste("Removing unmatched file:", unMatchedFiles[i]))
      file.remove(unMatchedFiles[i])
    }
  }

  ## Check if any .h5s are not found in any paths. Then delete them.
  allPaths =unlist(lapply(proj@projectMetadata$GroupCoverages, function(x) x$coverageMetadata$File))
  if(any(!allH5s %in% allPaths) & removeOld){
    warning("Removing unmatched .h5 files")
    for(i in seq_along(allH5s)){
      if(!allH5s[i] %in% allPaths){
        file.remove(allH5s[i])
      }
    }
  } else if(any(!allH5s %in% allPaths)){
    warning("Some .h5 files not found in paths")
  }else{
    message("All .h5 files found in paths")
  }

  saveArchRProject(proj)
  return(proj)
}

##########################################################################################
# Helper For cluster identity
##########################################################################################

#' Create a Confusion Matrix based on two value vectors
#'
#' This function creates a confusion matrix based on two value vectors.
#'
#' @param i A character/numeric value vector to see concordance with j.
#' @param j A character/numeric value vector to see concordance with i.
#' @export
confusionMatrix <- function(
  i = NULL, 
  j = NULL
  ){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}


#' Re-map a character vector of labels from an old set of labels to a new set of labels
#'
#' This function takes a character vector of labels and uses a set of old and new labels
#' to re-map from the old label set to the new label set.
#'
#' @param labels A character vector containing lables to map.
#' @param newLabels A character vector (same length as oldLabels) to map labels to from oldLabels.
#' @param oldLabels A character vector (same length as newLabels) to map labels from to newLabels
#' @export
mapLabels <- function(labels = NULL, newLabels = NULL, oldLabels = names(newLabels)){

  .validInput(input = labels, name = "labels", valid = c("character"))
  .validInput(input = newLabels, name = "newLabels", valid = c("character"))
  .validInput(input = oldLabels, name = "oldLabels", valid = c("character"))

  if(length(newLabels) != length(oldLabels)){
    stop("newLabels and oldLabels must be equal length!")
  }

  if(!requireNamespace("plyr", quietly = TRUE)){
    labels <- paste0(labels)
    oldLabels <- paste0(oldLabels)
    newLabels <- paste0(newLabels)
    labelsNew <- labels
    for(i in seq_along(oldLabels)){
        labelsNew[labels == oldLabels[i]] <- newLabels[i]
    }
    paste0(labelsNew)
  }else{
    paste0(plyr::mapvalues(x = labels, from = oldLabels, to = newLabels))
  }

}







