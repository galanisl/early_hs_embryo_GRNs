rm(list=ls())

source("localDirs.R")
adj<-1e-15
nAll<-c(10*(1:9),100*(1:10))
pAll<-10*(1:10)
nRep<-250
r<-0.07


###load functions and packages
CRANpackages<-c("dagitty")
for(iP in 1:length(CRANpackages)){
  tempPackage<-CRANpackages[iP]
  pckg = try(require(tempPackage,character.only=T))
  if(!pckg){
    print(paste0("Installing '",tempPackage,"' from CRAN"))
    install.packages(tempPackage,repos="http://cran.r-project.org")
    require(tempPackage,character.only=T)
  }
}

set.seed(1)

syntheticData<-lapply(1:nRep,function(iRep){
	print(paste0(iRep,"/",nRep))
	return(lapply(pAll,function(p){	
		return(lapply(nAll,function(n){
			tempGTedges<-data.frame(matrix(NA,nrow=0,ncol=2))
			while(nrow(tempGTedges)==0){
				tempG<-randomDAG(p,r)
				tempGTedges<-dagitty::edges(tempG)
			}
			tempMat<-t(as.matrix(simulateSEM(tempG,b.default=1,N=n,standardized=F)))		
			return(list(GT=tempG,data=tempMat))
		}))
	}))
})

save(list=c("syntheticData","nAll","pAll","r"),file=paste0(dataDir,"/TFnetSim.rd"))
