library(ArchR)
set.seed(1)

## Setting default number of Parallel threads to 8
addArchRThreads(threads = 8)
addArchRGenome("mm10")

proj <- loadArchRProject(path = "GutProject")

Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

cellNames <- rownames(Cdata)[which(Cdata$superClusters%in%c("lung","stomach","esophagus","pharynx"))]

subsetArchRProject(
  ArchRProj = proj,
  cells = cellNames,
  outputDirectory = "GutProject_foregut",
  dropCells = TRUE,
  logFile = NULL
)

proj_fore <- loadArchRProject(path = "GutProject_foregut")

Cdata <- getCellColData(ArchRProj = proj_fore, select = NULL, drop = FALSE)

#find macs2 path
pathToMacs2 <- findMacs2()

proj_fore <- addGroupCoverages(ArchRProj = proj_fore, groupBy = "superClusters",force=T)

proj_fore <- addReproduciblePeakSet(
    ArchRProj = proj_fore, 
    groupBy = "superClusters", 
    pathToMacs2 = pathToMacs2,
    additionalParams = "--nomodel --nolambda -q 0.01 --extsize 100",
    force=TRUE
)

proj_fore <- saveArchRProject(ArchRProj = proj_fore)

forePeaks <- getPeakSet(proj_fore)

names(forePeaks) <- 1:154731
forePeaks_df <- as.data.frame(forePeaks)
forePeaks_df <- forePeaks_df[,1:3]
forePeaks_df[,2] <- forePeaks_df[,2]-1
write.table(forePeaks_df,"GutProject_foregut/GutProject_foregut_peakset.bed",row.names=F,col.names=F,sep="\t",quote=F)


