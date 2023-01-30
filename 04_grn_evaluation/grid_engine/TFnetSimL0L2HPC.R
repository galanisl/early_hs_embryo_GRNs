###preliminaries
args<-commandArgs(trailingOnly=F)
taskID<-as.numeric(gsub("-","",args[length(args)],fixed=T))
source("/home/ucbptba/R/RMPIsetup.R")


###parameters
numCores<-mpi.comm.size()-1
nRep<-numCores
nBinner<-nRep
adj<-1e-15


###folders
dataDir<-"/home/ucbptba/data"
saveDir<-"/home/ucbptba/Scratch/output"
libDir<-"/home/ucbptba/Scratch/localRpackages"
archiveDir<-"/home/ucbptba/Scratch/data"

###load functions and packages
mpi.bcast.Robj2slave(libDir,comm=1,all=F)
require("igraph")
mpi.remote.exec(require("igraph"))
require("Hmisc")
mpi.remote.exec(require("Hmisc"))
require("dagitty",lib.loc=libDir)
mpi.remote.exec(require("dagitty",lib.loc=libDir))
require("L0Learn",lib.loc=libDir)
mpi.remote.exec(require("L0Learn",lib.loc=libDir))


### process over cluster	
workerFunction<-function(iRep){

	set.seed(iRep)
	
	##load data
	load(paste0(archiveDir,"/TFnetSim.rd"))
	
	GTedgeList<-list()
	edgeList<-list()
	aucValidSym<-list()
	
	for (j in 1:length(pAll)){

	   GTedgeList[[j]]<-list()
	   edgeList[[j]]<-list()
	   aucValidSym[[j]]<-list()
	  
		for (i in 1:length(nAll)){
	
			print(paste0("p=",pAll[j],"; n=",nAll[i]))
			tempGTedges<-as.matrix(dagitty::edges(syntheticData[[iRep]][[j]][[i]][["GT"]])[,c(1:2)])	
	  		GTedgeList[[j]][[i]]<-tempGTedges 		
	  		tempData<-syntheticData[[iRep]][[j]][[i]][["data"]]
	  		tempEdges<-do.call(rbind,lapply(rownames(tempData),function(tempFeature){
	  			tempParents<-setdiff(rownames(tempData),tempFeature)
	  			tempY<-tempData[tempFeature,]
	  			tempY<-tempY/sd(tempY)
	  			tempX<-tempData[tempParents,]
	  			tempX<-t(tempX/pmax(apply(tempX,1,sd),1e-6))
	  			suppressWarnings(tempFit<-tryCatch(L0Learn.cvfit(tempX,tempY,penalty="L0L2",nFolds=3,algorithm="CDPSI",maxSuppSize=length(tempParents)),error=function(e){return(NA)}))
				if(!all(is.na(tempFit))){
					bestGammaLoc<-which.min(lapply(tempFit[["cvMeans"]],min))
					bestLambdaLoc<-which.min(tempFit$cvMeans[[bestGammaLoc]])
					bestLambda<-tempFit$fit$lambda[[bestGammaLoc]][bestLambdaLoc]
					bestGamma<-tempFit$fit$gamma[bestGammaLoc]
					tempCoeffs<-as.matrix(coef(tempFit,lambda=bestLambda,gamma=bestGamma))[,1]
					tempCoeffs<-tempCoeffs[names(tempCoeffs)!="Intercept"]
				}else{
					tempCoeffs<-rep(0,length(tempTFs))
				}
				names(tempCoeffs)<-tempParents
				tempCoeffs<-data.frame(FROM=names(tempCoeffs),TO=rep(tempFeature,length(tempCoeffs)),STAT=abs(as.numeric(tempCoeffs)),stringsAsFactors=F)
				return(tempCoeffs)		
	  		}))
	  		edgeList[[j]][[i]]<-tempEdges  		
	  		symEdgeVals<-matrix(0,nrow=nrow(tempData),ncol=nrow(tempData),dimnames=list(rownames(tempData),rownames(tempData)))
	  		rm(list="tempData")
	  		for(iRow in 1:nrow(tempEdges)){
	  			symEdgeVals[tempEdges[iRow,"FROM"],tempEdges[iRow,"TO"]]<-tempEdges[iRow,"STAT"]
	  		}
	  		rm(list="tempEdges")
			symEdgeVals<-pmax(symEdgeVals,t(symEdgeVals))
			symGT<-symEdgeVals*0
			tempGraph<-as.matrix(as_adjacency_matrix(graph_from_edgelist(tempGTedges,directed=T),type="both",sparse=F))
			rm(list="tempGTedges")
			symGT[rownames(tempGraph),colnames(tempGraph)]<-tempGraph
			rm(list="tempGraph")
			symGT<-pmax(symGT,t(symGT))	
			symEdgeVals<-symEdgeVals[upper.tri(symEdgeVals)]
			symGT<-symGT[upper.tri(symGT)]
			aucValidSym[[j]][[i]]<-(somers2(symEdgeVals,symGT)[["Dxy"]]+1)/2
			rm(list=c("symEdgeVals","symGT"))
		}
		
	}

	aucValidSym<-do.call(cbind,lapply(aucValidSym,unlist,recursive=F))
	rownames(aucValidSym)<-paste(nAll)
	colnames(aucValidSym)<-paste(pAll)
	return(aucValidSym)

}


mpi.bcast.Robj2slave(comm=1,all=T)
mpi.bcast.Rfun2slave(comm=1)

procTime<-system.time({	
	aucValidSymAll<-mpi.apply(1:nBinner,workerFunction)
})[[3]]
print(paste("Time: ",signif(procTime/3600,digits=2)," hours",sep=""))
	
load(paste0(archiveDir,"/TFnetSim.rd"))
		
### save results
save(list=c("aucValidSymAll","nAll","pAll","r"),file=paste0(saveDir,"/TFnetSimL0L2res.rd"))

### ending

.Last()
