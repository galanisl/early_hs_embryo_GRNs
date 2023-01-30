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
require("BiocVersion",lib.loc=libDir)
mpi.remote.exec(require("BiocVersion",lib.loc=libDir))
require("BiocManager",lib.loc=libDir)
mpi.remote.exec(require("BiocManager",lib.loc=libDir))
require("GENIE3",lib.loc=libDir)
mpi.remote.exec(require("GENIE3",lib.loc=libDir))


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
		
			print(paste0("rep ",iRep,"/",nRep,"; p=",pAll[j],"; n=",nAll[i]))		
	  		tempGTedges<-as.matrix(dagitty::edges(syntheticData[[iRep]][[j]][[i]][["GT"]])[,c(1:2)])
	  		GTedgeList[[j]][[i]]<-tempGTedges
	  		tempData<-syntheticData[[iRep]][[j]][[i]][["data"]]
			GENIE3regs<-GENIE3(tempData)
			rm(list="tempData")
			GENIE3regs[is.na(GENIE3regs)]<-0
			GENIE3regs<-abs(GENIE3regs)
			tempGraph<-graph_from_adjacency_matrix(as.matrix(GENIE3regs),mode="directed",weighted=T)
			edgeList[[j]][[i]]<-data.frame(as_edgelist(tempGraph),STAT=E(tempGraph)$"weight",stringsAsFactors=F)
			symEdgeVals<-pmax(GENIE3regs,t(GENIE3regs))
			rm(list="GENIE3regs")
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
save(list=c("aucValidSymAll","nAll","pAll","r"),file=paste0(saveDir,"/TFnetSimGENIE3res.rd"))

### ending

.Last()
