
#' Remove unneeded matrices from Arrow File. 
#' 
#' This function will iterate over the arrows within an ArchRProject, only keeping the matrices selected. 
#' Use this function if you want a minimal ArchR Project, or to recover a corrupted ArchR Project. 
#' Sometimes, an ArchR project will become corrupted. This can occur when a multithreading application gets interrupted, 
#' leaving some arrows with matrices, that others do not have. Removing those problematic matrices from the ArchR Project will fix any subsetting issues.
#' The minimum matrix to keep is is the TileMatrix - keep must include the TileMatrix, or else DietArchRProject won't run.
#'
#' @param input An `ArchRProject` object or a character vector containing the paths to the ArrowFiles to be used.
#' @param keep A vector of matrices to keep. You can check on all available matrices via getAvailableMatrices(). If an error occurs during copying or subsetting, some matrices may only appear within a subset of arrow files. Dieting the ArchR project down to just the available matrices will remove all matrices that do not occur in all arrow files. 
#' @param verbose A boolean value that determines whether standard output is printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
##

## 
DietArchRProject <- function(ArchRProj = NULL, keep = "TileMatrix", verbose = FALSE, logFile = createLogFile("DietArchRProject")){	
	if(! any(keep %in% "TileMatrix")){
		stop("You must keep at least the TileMatrix within each arrow file.")
	}else{
		allMatrices <- getAvailableMatrices(ArchRProj)
		for(i in 1:length(allMatrices)){
			if(!allMatrices[i] %in% keep){
				.dropGroupsFromProject(ArchRProj, dropGroups = allMatrices[i], verbose = verbose)
			}
		}
	}
    
	return(ArchRProj)
}


#################################################################################

##### This function will remove all of a given matrix from an ArchRProject
.dropGroupsFromProject <- function(
	ArchRProj = NULL,
	dropGroups = NULL,
	level = 0,
	verbose = FALSE,
	logFile = NULL){

	arrowFiles <- getArrowFiles(ArchRProj)
	if(getArchRThreads() == 1){

		for(i in 1:length(arrowFiles)){
			.dropGroupsFromArrow(arrowFiles[i], dropGroups=dropGroups, level = level, verbose = verbose, logFile = logFile)
		}
	}else{

		mclapply(arrowFiles, function(x){
			.dropGroupsFromArrow(x, dropGroups=dropGroups, level = level, verbose = verbose, logFile = logFile)
		}, mc.cores = getArchRThreads())
	}
}