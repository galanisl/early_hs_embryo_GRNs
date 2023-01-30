rm(list=ls())

adj<-1e-15
source("localDirs.R")
topG<-500
topT<-25
nSamp<-100


CRANpackages<-c("vioplot","epitools","igraph","openxlsx")
for(iP in 1:length(CRANpackages)){
  tempPackage<-CRANpackages[iP]
  pckg = try(require(tempPackage,character.only=T))
  if(!pckg){
    print(paste0("Installing '",tempPackage,"' from CRAN"))
    install.packages(tempPackage,repos="http://cran.r-project.org")
    require(tempPackage,character.only=T)
  }
}

scale<-function(tempVar){
	return((tempVar-mean(tempVar))/(max(sd(tempVar),1e-6)))
}

GSEA<-function(sigFeatures,allFeatures,GSdefs){
  sigTemp<-allFeatures%in%sigFeatures
  system.time(GSEA<-do.call(rbind,lapply(names(GSdefs),function(GSname){#apply enrichment test to each gene set
    GS<-GSdefs[[GSname]];#genes in gene set
    inGS<-allFeatures%in%GS;#is each gene in the gene set?
    contTable<-array(c(sum(sigTemp&inGS),sum(sigTemp&!inGS),sum(!sigTemp&inGS),sum(!sigTemp&!inGS)),c(2,2));#contingency table for the Fisher tests
    FisherTest<-oddsratio.fisher(contTable);#Fisher's exact test
    FisherTestP<-fisher.test(contTable,alternative="greater")$p.value;#one-sided Fisher test pvalue
    sigGeneNames<-paste(allFeatures[sigTemp&inGS],collapse="; ");#symbols of significant genes
    res<-c(FisherTest$measure[2,1],FisherTest$measure[2,2:3],FisherTestP);
    return(res)
  })))[[3]]->tempTime
  colnames(GSEA)<-c("OR","OR lower 95% CI","OR upper 95% CI","p-val")
  rownames(GSEA)<-names(GSdefs)
  print(paste0("GSEA: t=",signif(tempTime,digits=2),"s"))
  GSEA<-GSEA[order(as.numeric(GSEA[,"p-val"])),];#sort by p-val
  return(GSEA)
}

invisible(lapply(c("tpm","bc","fpkm","lc"),function(normMethod){

#  invisible(lapply(c(F,T),function(doImpute){
doImpute<-F  
    if(doImpute){
      impSufx<-"Imp"
    }else{
      impSufx<-""
    }
    
    load(paste0(dataDir,"/TFnetPreImpEmbrData",impSufx,"_",normMethod,".rd"))
    load(paste0(fileDir,"/TFnetRrnd",impSufx,"_",normMethod,".rd"))
    
    invisible(sapply(c("Epiblast","PE","TE"),function(tempType){
      methodTabs<-sapply(c("TFnetL0L2","TFnet0L0L2","TFnetGENIE3","TFnet0GENIE3","TFnetMI","TFnet0MI","TFnetCor","TFnet0Cor"),function(tempMethod){
        
        load(paste0(saveDir,"/",tempMethod,"res",impSufx,"_",normMethod,"_",tempType,".rd"))
        if(tempMethod%in%c("TFnet0GENIE3","TFnet0MI","TFnet0Cor")){
          allGRNregsTab<-allGRNregsTab0
        }
        print(paste0(tempType,"; ",impSufx," ",normMethod,"; ",tempMethod,", ",nrow(allGRNregsTab)))
        allGRNregsTab<-allGRNregsTab[order(abs(allGRNregsTab[,"SCORE"]),decreasing=T),]
        print(paste0(tempType,"; ",impSufx," ",normMethod,"; ",tempMethod,", ",nrow(allGRNregsTab)))
        return(allGRNregsTab)
      },simplify=F,USE.NAMES=T)
      
      topRstats<-sapply(methodTabs,function(tempTabs){
        return(tempTabs[,"R"])
      },simplify=F,USE.NAMES=T)
      topSstats<-sapply(methodTabs,function(tempTabs){
        return(scale(abs(tempTabs[,"SCORE"])))
      },simplify=F,USE.NAMES=T)
      topRstats[["TFnetRND"]]<-allSampR[[tempType]][["sampR"]]
      topRstats[["TFnet0RND"]]<-allSampR[[tempType]][["samp0R"]]
      
      pdf(file=paste0(plotDir,"/TFnetRviols",impSufx,"_",normMethod,"_",tempType,".pdf",sep=""),width=12,height=9);{
        par(mar=c(12,5,3,1),xpd=NA,las=3)
        vioplot(topRstats,col=c("turquoise4","turquoise2","chartreuse4","chartreuse2","tomato4","tomato2","goldenrod4","goldenrod2","darkgrey","lightgrey"),lwd=1.5,colMed="white",colMed2="grey 75",cex.axis=1.9)
        mtext("R",side=2,line=3,cex=1.9)
        mtext(paste(tempType,normMethod,impSufx),side=3,line=0.5,cex=1.9,las=1)
      }
      dev.off()

      pdf(file=paste0(plotDir,"/TFnetSviols",impSufx,"_",normMethod,"_",tempType,".pdf",sep=""),width=10,height=9);{
        par(mar=c(12,5,3,1),xpd=NA,las=3)
        tempPlot<-tryCatch(vioplot(topSstats,col=c("turquoise4","turquoise2","chartreuse4","chartreuse2","tomato4","tomato2","goldenrod4","goldenrod2"),lwd=1.5,colMed="white",colMed2="grey 75",cex.axis=1.9),error=function(e){return(NA)})
        if(all(is.na(tempPlot))){
        	topSstats<-lapply(topSstats,function(temp){return(temp+rnorm(length(temp),mean=0,sd=1e-12))})
        		vioplot(topSstats,col=c("turquoise4","turquoise2","chartreuse4","chartreuse2","tomato4","tomato2","goldenrod4","goldenrod2"),lwd=1.5,colMed="white",colMed2="grey 75",cex.axis=1.9)
        }
        mtext("Standardised network score",side=2,line=3,cex=1.9)
        mtext(paste(tempType,normMethod,impSufx),side=3,line=0.5,cex=1.9,las=1)
      }
      dev.off()
      
      return(NULL)
      
    }))
    
    allTables<-sapply(c("Epiblast","PE","TE"),function(tempType){
      methodTabs<-sapply(c("TFnetL0L2","TFnet0L0L2","TFnetGENIE3","TFnet0GENIE3","TFnetMI","TFnet0MI","TFnetCor","TFnet0Cor"),function(tempMethod){
        
        load(paste0(saveDir,"/",tempMethod,"res",impSufx,"_",normMethod,"_",tempType,".rd"))
        if(tempMethod%in%c("TFnet0GENIE3","TFnet0MI","TFnet0Cor")){
          allGRNregsTab<-allGRNregsTab0
        }
        print(paste0(tempType,"; ",impSufx," ",normMethod,"; ",tempMethod,", ",nrow(allGRNregsTab)))
        allGRNregsTab<-allGRNregsTab[order(abs(allGRNregsTab[,"SCORE"]),decreasing=T),]
        print(paste0(tempType,"; ",impSufx," ",normMethod,"; ",tempMethod,", ",nrow(allGRNregsTab)))
        return(allGRNregsTab)
      },simplify=F,USE.NAMES=T)
      
      return(methodTabs)
      
    },simplify=F,USE.NAMES=T)

    workBook<-createWorkbook()
    
    invisible(sapply(c("TE","Epiblast"),function(tempType){

      if(tempType=="TE"){
        geneSets<-read.csv(paste0(dataDir,"/TE_gene_sets.csv"),as.is=T,stringsAsFactors=F)
      }
      if(tempType=="Epiblast"){
        geneSets<-read.csv(paste0(dataDir,"/EPI_gene_sets.csv"),as.is=T,stringsAsFactors=F)
      }
      
      GSnames<-geneSets[,"id"]
      GSdefs<-strsplit(geneSets[,"genes"],split=",",fixed=T)
      names(GSdefs)<-GSnames
      
      allGenes<-rownames(allExprChrom[[tempType]][["expr"]])
      GSdefs<-sapply(GSdefs,function(temp){return(temp[temp%in%allGenes])})
      GSnames<-names(GSdefs)
      
      topVstats<-sapply(c("TFnetL0L2","TFnet0L0L2","TFnetGENIE3","TFnet0GENIE3","TFnetMI","TFnet0MI","TFnetCor","TFnet0Cor"),function(tempMethod){     
        allGRNregsTab<-allTables[[tempType]][[tempMethod]]
        allGRNregsTab<-allGRNregsTab[order(abs(allGRNregsTab[,"SCORE"]),decreasing=T),]
        
        tempGraph<-graph_from_edgelist(as.matrix(allGRNregsTab[,c("FROM","TO")]),directed=T)
        E(tempGraph)$weight<-as.matrix(allGRNregsTab[,"SCORE"])	
        tempAdj<-as.matrix(as_adjacency_matrix(tempGraph,attr="weight")[unique(allGRNregsTab[,"FROM"]),unique(allGRNregsTab[,"TO"])])
        topTFs<-head(names(sort(rowSums(abs(tempAdj)>0),decreasing=T)),topT)
        tempV<-sapply(topTFs,function(tempTF){
          tempRegs<-tempAdj[tempTF,]
          tempRegs<-tempRegs[abs(tempRegs)>0]
          tempRegs<-tempRegs[order(abs(tempRegs),decreasing=T)]
          tempRegs<-head(tempRegs,topG)
          print(paste0(tempType,", ",tempMethod,impSufx,", ",normMethod,", TF ",match(tempTF,topTFs),"/",length(topTFs),": ",length(tempRegs)))
          tempGSEA<-GSEA(names(tempRegs),allGenes,GSdefs)
          return(tempGSEA[,c("OR","p-val")])
        },simplify=F,USE.NAMES=T)
        return(tempV)
      },simplify=F,USE.NAMES=T)

      topVstats<-sapply(c("TFnetL0L2","TFnet0L0L2","TFnetGENIE3","TFnet0GENIE3","TFnetMI","TFnet0MI","TFnetCor","TFnet0Cor"),function(tempMethod){
        tempTable<-do.call(cbind,lapply(names(topVstats[[tempMethod]]),function(tempName){
          temp<-topVstats[[tempMethod]][[tempName]]
          colnames(temp)<-paste(tempName,colnames(temp))
          return(temp)
        }))
        addWorksheet(workBook,paste0(tempType,"_",tempMethod))
        writeData(workBook,paste0(tempType,"_",tempMethod),tempTable,startRow=1,startCol=1)
        return(sapply(topVstats[[tempMethod]],function(temp){return(temp[,"p-val"])},simplify=F,USE.NAMES=T))
      },simplify=F,USE.NAMES=T)
      
      regGenes<-intersect(rownames(allExprChrom[[tempType]][["expr"]]),unique(allExprChrom[[tempType]][["chrom"]][,2]))
      
      topVstatsRND<-lapply(1:nSamp,function(tempI){
        print(paste0(tempI,"/",nSamp))
        topGenes<-sample(regGenes,topG)
        tempGSEA<-GSEA(topGenes,rownames(allExprChrom[[tempType]][["expr"]]),GSdefs)
        return(tempGSEA[,c("OR","p-val")])
      })
      
      tempTable<-do.call(cbind,lapply(1:length(topVstatsRND),function(tempI){
        temp<-topVstatsRND[[tempI]]
        colnames(temp)<-paste0("rndSamp",tempI," ",colnames(temp))
        return(temp)
      }))
      
      addWorksheet(workBook,paste0(tempType,"_random"))
      writeData(workBook,paste0(tempType,"_random"),tempTable,startRow=1,startCol=1)
      
      topVstatsRND<-lapply(topVstatsRND,function(temp){
        return(temp[,"p-val"])
      })

      pdf(file=paste0(plotDir,"/enrichPhists",impSufx,"_",normMethod,"_",tempType,".pdf"),height=25,width=10);{
        par(mfrow=c(5,2),mar=c(5,5,2,1))
        invisible(sapply(c("TFnetL0L2","TFnet0L0L2","TFnetGENIE3","TFnet0GENIE3","TFnetMI","TFnet0MI","TFnetCor","TFnet0Cor"),function(tempMethod){
          tempP<-unlist(topVstats[[tempMethod]])
          tempHist<-hist(tempP,breaks=10,plot=F)
          tempHist[["counts"]]<-tempHist[["counts"]]/length(topVstats[[tempMethod]])
          plot(tempHist,freq=T,cex.axis=1.4,main="",xlab="",ylab="",ylim=c(0,200),col="lightgray")
          mtext(bquote(italic(p)*"-val"),side=1,line=3.5,cex=1.5)
          mtext("Number of gene sets",side=2,line=3,cex=1.5)
          mtext(paste0(tempMethod,impSufx,", ",normMethod,": V=",tempHist[["counts"]][1]),side=3,line=-3,cex=1.5)
          return(NULL)
        }))
        tempP<-unlist(topVstatsRND)
        tempHist<-hist(tempP,breaks=10,plot=F)
        tempHist[["counts"]]<-tempHist[["counts"]]/nSamp
        plot(tempHist,freq=T,cex.axis=1.4,main="",xlab="",ylab="",ylim=c(0,200),col="lightgray")
        mtext(bquote(italic(p)*"-val"),side=1,line=3.5,cex=1.5)
        mtext("Number of gene sets",side=2,line=3,cex=1.5)
        
        mtext(paste0("Random samples of\n n=",topG," genes: V=",tempHist[["counts"]][1]),side=3,line=-3,cex=1.5)
        
      }
      dev.off()
      
    }))
      
    saveWorkbook(workBook,file=paste0(saveDir,"/TFnetV",normMethod,".xlsx"),overwrite=T)
      
#    return(NULL)
    
#  }))
  
  return(NULL)
  
}))
