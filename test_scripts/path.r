#target bin directory
targetDir <- "~/lib/bin/"

#list of scripts to link
scripts <- c(
	"csv2json.R","joinCounts.R","runLibraryQC.R",
	"calcEnrichment.R","runSelectionQC.R","mavevisLocal.R",
    "condenseQC.R", "scaleScores.R"
)

#find scripts folder in the local library installation
scriptsFolder <- system.file("scripts/",
	package = "tileseqMave",
	mustWork = TRUE
)

#function to create symlinks
linkScript <- function(scriptName,targetDir="~/lib/bin/") {
	scriptPath <- paste0(scriptsFolder,scriptName)
	if (!file.exists(scriptPath)) {
		stop(scriptName, "not found!")
	}
	file.symlink(from=scriptPath,to=paste0(targetDir,scriptName))
}

#run the function
lapply(scripts,linkScript,targetDir=targetDir)
