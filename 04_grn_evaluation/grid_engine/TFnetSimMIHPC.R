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
softDir<-"/home/ucbptba/software"
saveDir<-"/home/ucbptba/Scratch/output"
libDir<-"/home/ucbptba/Scratch/localRpackages"
archiveDir<-"/home/ucbptba/Scratch/data"

###load functions and packages
mpi.bcast.Robj2slave(libDir,comm=1,all=F)
require("igraph")
mpi.remote.exec(require("igraph"))
require("Hmisc")
mpi.remote.exec(require("Hmisc"))
require("splines")
mpi.remote.exec(require("splines"))
require("dagitty",lib.loc=libDir)
mpi.remote.exec(require("dagitty",lib.loc=libDir))



### process over cluster	
workerFunction<-function(iRep){

	set.seed(iRep)
	
	##load data
	load(paste0(archiveDir,"/TFnetSim.rd"))
	source(paste0(softDir,"/mutual_information-master/mutual_information2.R"))
	
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
			MIregs<-do.call(cbind,lapply(rownames(tempData),function(tempFeature){
				return(sapply(rownames(tempData),function(tempParent){
					return(tryCatch(mutual.information2(x1=tempData[tempParent,],x2=tempData[tempFeature,],bins=max(1,min(floor(ncol(tempData)/3),10)),spline.order=3),error=function(e){return(0)}))
				}))
			}))
			colnames(MIregs)<-rownames(tempData)
			rm(list="tempData")
			MIregs[is.na(MIregs)]<-0
			MIregs<-abs(MIregs)
			Zscores<-sapply(colnames(MIregs),function(tempName){
			  tempRegs<-MIregs[rownames(MIregs)!=tempName,tempName]
			  if(max(abs(tempRegs))<1e-12){
			    tempRes<-rep(0,length(tempRegs))
			  }else{
		      tempRes<-qnorm(0.5+1e-12+(0.5-2e-12)*ecdf(tempRegs)(tempRegs))
			  }
			  tempOut<-rep(0,nrow(MIregs))
			  names(tempOut)<-rownames(MIregs)
			  tempOut[names(tempOut)!=tempName]<-tempRes
				return(tempOut)
			})
			rownames(Zscores)<-rownames(MIregs)
			rm(list="MIregs")
			Zscores<-Zscores^2
			Zscores<-Zscores+t(Zscores)
			Zscores<-sqrt(Zscores)		
			pVals<-2*apply(Zscores,2,pnorm,lower.tail=F)
			rownames(pVals)<-rownames(Zscores)
			rm(list="Zscores")
			regStats<-(-log10(pVals))
			rm(list="pVals")		
			tempGraph<-graph_from_adjacency_matrix(as.matrix(regStats),mode="directed",weighted=T)
			edgeList[[j]][[i]]<-data.frame(as_edgelist(tempGraph),STAT=E(tempGraph)$"weight",stringsAsFactors=F)
			symEdgeVals<-pmax(regStats,t(regStats))
			rm(list="regStats")
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
save(list=c("aucValidSymAll","nAll","pAll","r"),file=paste0(saveDir,"/TFnetSimMIres.rd"))

### ending

.Last()
