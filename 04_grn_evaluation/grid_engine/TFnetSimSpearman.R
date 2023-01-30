rm(list=ls()) 

adj<-1e-15 

source("localDirs.R")

require("igraph")
require("Hmisc")
require("parallel")
require("dagitty")

numCores<-max(1,detectCores()-1)

load(paste0(dataDir,"/TFnetSim.rd"))

nRep<-length(syntheticData)

aucValidSymAll<-mclapply(1:nRep,function(iRep){

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
			corRegs<-cor(t(tempData),method="spearman")
			rm(list="tempData")
			corRegs[is.na(corRegs)]<-0
			corRegs<-abs(corRegs)
			tempGraph<-graph_from_adjacency_matrix(as.matrix(corRegs),mode="directed",weighted=T)
			edgeList[[j]][[i]]<-data.frame(as_edgelist(tempGraph),STAT=E(tempGraph)$"weight",stringsAsFactors=F)
			symEdgeVals<-pmax(corRegs,t(corRegs))
			rm(list="corRegs")
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
	
},mc.cores=numCores)
	

save(list=c("aucValidSymAll","nAll","pAll","r"),file=paste0(fileDir,"/TFnetSimSPEARMANres.rd"))
