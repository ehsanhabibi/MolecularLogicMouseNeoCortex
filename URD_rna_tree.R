source("seuratV3ToURD.R")

#Converte seurat to urd object
object.merged.urd <- seuratV3ToURD(seurat.object = scDevSC.seu)

#Calculate diffusion map
object.merged.urd <- calcDM(object.merged.urd, knn=200)
  
# Do the flood
# For our data set (~95K cells) this step was run on the cluster using the script "URD_Flood_server.R". In total, we performed 200 simulations. 
root.cells <- root_and_tips$root

#Next, we loaded floods from the cluster
floods <- lapply(list.files(path="floods/", pattern="flood-", full.names=T), readRDS)

#Then, we processed the simulations.
object.merged.urd <- floodPseudotimeProcess(object.merged.urd, floods, floods.name="pseudotime", max.frac.NA=0.4, pseudotime.fun=mean, stability.div=20)

pseudotimePlotStabilityOverall(object.merged.urd)

plotDim(object.merged.urd, "pseudotime", plot.title = "Pseudotime", point.size = 0.2, legend.title = "Pseudotime" )

gg.data <- cbind(object.merged.urd@pseudotime, object.merged.urd@meta[rownames(object.merged.urd@pseudotime),])

# Plot
ggplot(gg.data, aes(x=pseudotime, color=orig.ident, fill=orig.ident)) + geom_density(alpha=0.4) + theme_bw()


#Biased Random Walks 
diffusion.logistic <- pseudotimeDetermineLogistic(object.merged.urd, "pseudotime", optimal.cells.forward=40, max.cells.back=80, pseudotime.direction="<", do.plot=T, print.values=T)


#Run walks on cluster
biased.tm <- pseudotimeWeightTransitionMatrix(object.merged.urd, pseudotime = "pseudotime", logistic.params = diffusion.logistic, pseudotime.direction = "<")

########################################################################

#Tips is a list, names are terminal cell types like Tips[['glia']], Tips[[SCPN]], ... and each is a vector of terminal cells, the number might vary for each tip
Tips <- root_and_tips$tips

res.n <- length(names(Tips))

cluster.assignments <- data.frame(
  cluster=1:res.n,
  name=rep(NA, res.n),
  tip=rep(NA, res.n),
  row.names=1:res.n
)

tip.names <- names(Tips)

for (i in 1:length(tip.names)){
  
  cluster.assignments[i,"name"] <- tip.names[i]
  
}

# Remove any clusters that weren't assigned an identity
cluster.assignments <- cluster.assignments[!is.na(cluster.assignments$name),]

# Renumber clusters
cluster.assignments$cluster.new <- 1:nrow(cluster.assignments)

object.merged.urd@group.ids$`fs-Cluster` <- NA
object.merged.urd@group.ids$`fs-Cluster-Num` <- NA

# Copy cell identities over for each cluster
for (i in 1:nrow(cluster.assignments)) {
  cells <- Tips[[i]]
  object.merged.urd@group.ids[cells,"fs-Cluster"] <- cluster.assignments[i,"name"]
  object.merged.urd@group.ids[cells,"fs-Cluster-Num"] <- as.character(cluster.assignments[i,"cluster.new"])
}


# Define the root cells
root.cells <- root_and_tips$root

# Define the tip cells
tip.to.walk <- setdiff(unique(object.merged.urd@group.ids[,"fs-Cluster-Num"]), NA)


this.tip <- tips[as.numeric(tip.to.walk)] # tip.to.walk was passed by the cluster job array.

tip.cells <- rownames(object.merged.urd@group.ids)[which(object.merged.urd@group.ids[,"fs-Cluster-Num"] == this.tip)]
# Do the random walks on the cluster. We performed 350K walks per tip
these_walks <- simulateRandomWalk(start.cells=tip.cells, transition.matrix=biased.tm, end.cells=root.cells, n=walks.to.do, end.visits=1, verbose.freq=round(walks.to.do/20), max.steps=5000)



tip.walk.files <- list.files(path = path, pattern = paste0("walks-"), full.names = T)

walks <- list()

for (tip in seq){
  print(tip)
  tip.files <- list.files(path = path, pattern = paste0("walks-",tip,"\\."), full.names = T)
  walks[[tip]] <- unlist(lapply(tip.files, readRDS), recursive = F)
}

save(walks,file=paste0(path,"/walks_",walks.to.do/1000,"K.obj"))

names(walks) <- names(tip.cells)
  
object.merged.urd_after_adding_walks <- processRandomWalksFromTips(object.merged.urd, walks, verbose=T)
  
library(rgl)
  
object.merged.urd_after_adding_walks <- loadTipCells(object.merged.urd_after_adding_walks, tips="fs-Cluster-Num")
  

path2 = paste0(path,"tree/")

dir.create(path = path2)

divergence.method <- c("preference")
visit.threshold <- c(0.7)
minimum.visits <- c(2)
cells.per.pseudotime.bin <- c(80)
bins.per.pseudotime.window <- c(8)
p.threshs = c(0.01)


name = paste0(path2,"tree_v",visit.threshold,"_m",minimum.visits,"_b",bins.per.pseudotime.window,"_c",cells.per.pseudotime.bin,"_d",divergence.method,"_p",p.thresh)
print(name)
tree <- buildTree(object.merged.urd_after_adding_walks,
                              pseudotime = "pseudotime",
                              tips.use=names(Tips),
                              divergence.method = divergence.method,
                              cells.per.pseudotime.bin = cells.per.pseudotime.bin,
                              bins.per.pseudotime.window = bins.per.pseudotime.window,
                              save.all.breakpoint.info = T,
                              p.thresh=p.thresh,
                              minimum.visits=minimum.visits,
                              visit.threshold=visit.threshold,
                              min.cells.per.segment=10,
                              min.pseudotime.per.segment=.01,
                              verbose = T,
                              dendro.node.size=100,
                              dendro.cell.jitter=0.15,
                              dendro.cell.dist.to.tree=0.05,
                              use.only.original.tips=T,
                              weighted.fusion=T,
                              save.breakpoint.plots=name
            		)
            
tip.names <- unique(tree@group.ids[,c("fs-Cluster", "fs-Cluster-Num")])
tip.names <- tip.names[complete.cases(tip.names),]
tip.names$`fs-Cluster-Num` <- as.character(tip.names$`fs-Cluster-Num`)
            
tree <- nameSegments(tree,
                     segments= tip.names$`fs-Cluster-Num`,
                     segment.names = tip.names$`fs-Cluster`
	          )
