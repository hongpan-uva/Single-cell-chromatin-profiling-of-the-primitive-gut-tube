library(ArchR)
set.seed(1)

## Setting default number of Parallel threads to 8
addArchRThreads(threads = 8)
addArchRGenome("mm10")

proj <- loadArchRProject(path = "GutProject")

cellNames <- getCellNames(proj)
cellNames_r1 <- cellNames[which(substr(cellNames,1,7)=="E9.5_r1")]
cellNames_r2 <- cellNames[which(substr(cellNames,1,7)=="E9.5_r2")]

subsetArchRProject(
  ArchRProj = proj,
  cells = cellNames_r1,
  outputDirectory = "GutProject_r1",
  dropCells = TRUE,
  logFile = NULL
)

subsetArchRProject(
  ArchRProj = proj,
  cells = cellNames_r2,
  outputDirectory = "GutProject_r2",
  dropCells = TRUE,
  logFile = NULL,
)


proj_r1 <- loadArchRProject(path = "GutProject_r1")
#proj_r1 <- addClusters(input = proj_r1, reducedDims = "IterativeLSI", name = "Clusters_singleSample")
#proj_r1 <- saveArchRProject(ArchRProj = proj_r1)
Cdata_r1 <- getCellColData(ArchRProj = proj_r1, select = NULL, drop = FALSE)

proj_r2 <- loadArchRProject(path = "GutProject_r2")
#proj_r2 <- addClusters(input = proj_r2, reducedDims = "IterativeLSI", name = "Clusters_singleSample")
#proj_r2 <- saveArchRProject(ArchRProj = proj_r2)
Cdata_r2 <- getCellColData(ArchRProj = proj_r2, select = NULL, drop = FALSE)

proj <- loadArchRProject(path = "GutProject")
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

vec1 <- Cdata_r1$Clusters_singleSample
vec1 <- paste0("r1_",vec1)
names(vec1) <- rownames(Cdata_r1)

vec3 <- Cdata_r2$Clusters_singleSample
vec3 <- paste0("r2_",vec3)
names(vec3) <- rownames(Cdata_r2)

vec2 <- Cdata$newClusters
names(vec2) <- rownames(Cdata)

vec2_organ <- Cdata$superClusters
names(vec2_organ) <- rownames(Cdata)

save.image("individually_clustering.Rdata")

#Sankey diagram plotly
load("individually_clustering.Rdata")

linkmtx <- matrix(0,nrow=35,ncol=35)
organ_labels <- names(table(vec2_organ))
rownames(linkmtx) <- c(paste0("r1_C",1:13),organ_labels,paste0("r2_C",1:13))
colnames(linkmtx) <- c(paste0("r1_C",1:13),organ_labels,paste0("r2_C",1:13))

#source: single sample clustering
#target: whole clustering
for(i in paste0("r1_C",1:13)){
  for(j in organ_labels){
    cellnames1 <- names(vec1[which(vec1==i)])
    cellnames2 <- names(vec2_organ[which(vec2_organ==j)])
    linkmtx[i,j] <- length(intersect(cellnames1,cellnames2))
  }
}

for(i in paste0("r2_C",1:13)){
  for(j in organ_labels){
    cellnames1 <- names(vec3[which(vec3==i)])
    cellnames2 <- names(vec2_organ[which(vec2_organ==j)])
    linkmtx[i,j] <- length(intersect(cellnames1,cellnames2))
  }
}

nodelist = list(
  label = c(paste0("r1_C",1:13),organ_labels,paste0("r2_C",1:13)),
  color = c(rep("navy",13),rep("black",9),rep("#cd2626",13)),
  x = c(rep(0.1,13),rep(0.5,9),rep(0.9,13)),
  y = c(rep(0.1,13),rep(0.1,9),rep(0.1,13)),
  pad = 15,
  thickness = 20,
  line = list(
    color = "black",
    width = 0.5
  )
)

nodeindex = 0:34
names(nodeindex) = c(paste0("r1_C",1:13),organ_labels,paste0("r2_C",1:13))

#get linklist
sourceVec=c()
targetVec=c()
valueVec=c()

for(i in  0:34){
  for(j in 0:34){
    if(linkmtx[i+1,j+1]>0){
      sourceVec <- c(sourceVec, i)
      targetVec <- c(targetVec, j)
      valueVec <- c(valueVec, linkmtx[i+1,j+1])
    }
  }
}

linklist = list(
  source = sourceVec,
  target = targetVec,
  value = valueVec
)


#plot sankey diagram
library(plotly)

nodelist$label <- rep("",35)

fig <- plot_ly(
  type = "sankey",
  arrangement = "snap",
  orientation = "h",
  node = nodelist,
  link = linklist,
)
  
fig <- fig %>% layout(
    title = "Basic Sankey Diagram",
    font = list(
      size = 10
    )
)

fig
