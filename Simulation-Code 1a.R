#
# Carla de la Fuente Canto and Yves Vigouroux,  2022
# Evaluation of nine statistics to identify QTL in bulk segregant analysis using next generation sequencing approaches
#
#
# This script was developed to test different statistical methods for NGS-based Bulk Segregant Analysis (BSA). We used R version 3.6.3 [64bits]
# Here we compile the steps explained in the Material and Methods of the current paper grouped by the tasks defined in point 3 and 4:
#
# 1.  Code for simulation of segregating population and BSA population: 
#	    simulation of genotype, phenotype and bulk size for High and Low bulks for BSA.
#     Calculation of bulks allele frequencies based on allele depth
#
# 2.  Code for running statistical methods for BSA-QTL mapping.
#     A total of nine statistics are considered in the analysis of allele frequency differences between the High and Low bulks. 
#     Four statistics are single-marker based calculations (deltaSNP, G, EDm and LOD).
#     The other five statistics (TdeltaSNP, Gprime, ED100^4, SmLOD and AFDexp) consist on the smoothed version single-marker
#     based stats referred to a sliding window of consecutive SNPs across the genome (W). 
#
# 3.  Codes to evaluate the detection of QTLs with different statistical methods in three recombination rate model species
#     
#     Supplementary Code 1a. Compare the mapping of a QTL using different statistical methods across three recombination models.
#         In this case we simulate a single QTL in the middle locus of the chromosome and sequencing noise is simulated according to
#         a binomial distribution. 
#     Supplementary Code 1b runs the same simulation in 1a but considering real data from rice in the simulation of sequencing noise.
#     Supplementary Code 2a. Evaluate the accuracy of QTL detection using different statistics. In this case we use one thousand runs
#         of simulation. Each loop of simulation defines a random QTL position. The absolute difference between the initial position 
#         of the QTL and the position retrieved with each method is used to compare the accuracy of the different methods.This code uses binomial
#         distribution in the simulation of sequencing noise.
#     Supplementary Code 2b runs the same simulation in 2a but considering real data from rice in the simulation of sequencing noise.
# 
# 4. Supplementary Code 3: Define confidence intervals for each method combined with each recombination model. We use ten thousand runs
#         of simulations with no QTL effect. The 95% quantiles were used as significant threshold values to define 
#         the confidence intervals of each test statistic.
# 
#
#
########################################################################################################################
#####
##### Code 1a. Compare the mapping of a QTL using different statistical methods across three recombination models ----
#####
########################################################################################################################

# Preliminary steps: ----
## a - Load R libraries 
## b - Define code parameters 

  ## a - R libraries ----
  library(QTLseqr)
  library(DescTools)
  library(modeest)
  library(locfit)
  library(doParallel)
  library(plyr)
  library(dplyr)
  library(rootSolve)
  library(data.table) 
  library(foreach)

setwd("C:/...") # Define a working directory and copy there the rice vcf file "OgOb-all-merged.DPech.recode.vcf"
  
  ## b - Code parameters ----
  nbsim<-1   # Number of simulations 1
  nbchr<-1000   # Number of chromosomes in in a population of 500 individuals
  loci<-10000   # Number of loci or markers in the model chromosome
  qtleffect<-1  # QTL effect equivalent to one time the standard deviation of a normally distributed trait
  # qtleffect<-0.5 # QTL effect equivalent to half standard deviation of a normally distributed trait
  
# 1.Population design and case studies ----
##  1.1 and 1.2 Define GENOTYPE and PHENOTYPE considering three recombination models: Pearl millet, Rice and Foxtail millet.
##  1.3 Link genotype and phenotype defining a QTL and Bulks of contrasted lines

for (n in 1:nbsim){ 
  
  # 1.1 Genotype ---- 
  # Simulate the number of recombination events (q) in a set of 'nbchr' chromosomes from a F2 population derived
  # from bi-parental cross between diploid homozygous lines. We consider recombination frequency of an average chromosome of
  # 90 cM length in Pearl millet (Moumouni et al., 2015), 130 cM length in Rice (IRGSP, 2005) and 215 cM length in Foxtail 
  # millet (Masumoto et al., 2016). 
  # Recombination events follow a Poisson distribution λ = L (Haldane, 1919) where L is the length of the chromosome in Morgan. 
  # Crossovers location are uniformly distributed over the interval (0, L).
  # 0s and 1s were used to code the reference and alternate parental alleles respectively
  
  q<-rpois(nbchr,0.90) #Pearl millet 
  #q<-rpois(nbchr,1.30) #Rice 
  #q<-rpois(nbchr,2.15) #Foxtail millet
  
  ## Genotype matrix:
  # If the q value is >0 (i.e there is at least one crossover), define a random position for the crossovers (=vec)
  # Fill each chromosome with 0s or 1s in the intervals defined by 'vec'
  
  geno<-matrix(ncol = nbchr, 
               nrow = loci, 
               dimnames = list(c(1:loci), 
                               c(1:nbchr)))
  
  for (j in 1:nbchr){ 
    
    allele<-ifelse(runif(1) < 0.5, 0, 1)
    geno[,j]<-allele
    
    if (q[j]>0) {
      vec<-c(sort(sample(loci, q[j])),loci)
      
      for (k in 1 : (length(vec)-1)) {
        
        allele<-ifelse(allele==1,0,1)
        geno[vec[k]:vec[k+1],j]<-allele
        
      }
    }
    
    
  }
  
  # 1.2 Phenotype ----
  ## Simulate normally distributed phenotype for a group of segregant diploid individuals
  pheno<-rnorm(nbchr/2,0,1)
  phenoqtl<-pheno
  
  # 1.3 QTL and Bulks----
  qtlpos<-5000 # QTL position in locus 5000 for task 3.1

  # Define contrasted BULKS of high and low phenotype taking the 10% extreme individuals with extreme phenotype
  # QTL additive effect:
  for (i in 1: (nbchr/2)) {
    phenoqtl[i]<-pheno[i]+geno[qtlpos,2*i-1]*qtleffect +geno[qtlpos,2*i]*qtleffect
  }
  
  individuals<-c(1:(nbchr/2))
  bsa<-data.frame(individuals, phenoqtl)
  bsa<-bsa[order(bsa$phenoqtl),]
  
  lowpheno<-bsa$individuals[1:(0.10*length(bsa$individuals))]
  poslow<-sort(c(2*lowpheno,2*lowpheno-1))
  lowbulk<-geno[,poslow]
  
  highpheno<-bsa$individuals[((0.90*length(bsa$individuals))+1):length(bsa$individuals)]
  poshigh<-sort(c(2*highpheno,2*highpheno-1))
  highbulk<-geno[,poshigh]
  
  
  # Allele depths of the reference (REF) and alternate (ALT) genome in each bulk
  # Format the data according to QTLseqr input file format.
  # Allele depths and frequencies for each bulk in new data frame 'simdata'
  AD_ALT.LOW_i<-apply(lowbulk,1,function(x)sum(x == 1))
  AD_REF.LOW_i<-apply(lowbulk,1,function(x)sum(x == 0))
  DP.LOW<-AD_ALT.LOW_i+AD_REF.LOW_i
  
  AD_ALT.HIGH_i<-apply(highbulk,1,function(x)sum(x == 1))
  AD_REF.HIGH_i<-apply(highbulk,1,function(x)sum(x == 0))
  DP.HIGH<-AD_ALT.HIGH_i+AD_REF.HIGH_i
  
  simdata<-data.frame("chr",1:loci,AD_ALT.LOW_i,AD_REF.LOW_i,DP.LOW,AD_ALT.HIGH_i,AD_REF.HIGH_i,DP.HIGH)
  colnames(simdata)[1]<-"CHROM"
  colnames(simdata)[2]<-"POS" 
  
  #Add sequencing noise to the data according to a binomial distribution
  simdata$AD_ALT.LOW<-0
  simdata$AD_REF.LOW<-0
  simdata$AD_ALT.HIGH<-0
  simdata$AD_REF.HIGH<-0
  
  depth<-100 

  for (i in 1:nrow(simdata)){
    simdata$AD_ALT.LOW[i]<-rbinom(1,depth,(simdata$AD_ALT.LOW_i[i]/(simdata$AD_ALT.LOW_i[i]+simdata$AD_REF.LOW_i[i])))
    simdata$AD_REF.LOW[i]<-rbinom(1,depth,(simdata$AD_REF.LOW_i[i]/(simdata$AD_ALT.LOW_i[i]+simdata$AD_REF.LOW_i[i])))
    simdata$AD_ALT.HIGH[i]<-rbinom(1,depth,(simdata$AD_ALT.HIGH_i[i]/(simdata$AD_ALT.HIGH_i[i]+simdata$AD_REF.HIGH_i[i])))
    simdata$AD_REF.HIGH[i]<-rbinom(1,depth,(simdata$AD_REF.HIGH_i[i]/(simdata$AD_ALT.HIGH_i[i]+simdata$AD_REF.HIGH_i[i])))

  }
  
  simdata<-simdata[,-c(3:8)]
  
  write.csv(simdata,
            file = "simdata.csv",
            row.names = FALSE)
  
# 2. Statistical methods for BSA-QTL mapping ----
## Calculation of four marker-based statistics (deltaSNP, G, EDm and LOD) and five window-based statistics 
## for bandwith fixed window (TdeltaSNP, Gprime, SmLOD and AFDexp) and consecutive groups of markers (ED100^4)
## Calculations grouped by R package or R code published for the method used:
## 2.1 QTLseqr (Mansfeld and Grumet, 2018): deltaSNP, TdeltaSNP, G and Gprime
## 2.2 Euclidean Distance based statistics (Hill et al., 2003; Omboki et al., 2018; Zhang et al., 2019): EDm and ED100^4
## 2.3 Quantitative Gene Sequencing or QTG-Seq method (Zhang et al., 2019): LOD and Smooth LOD
## 2.4 Block Regression Mapping or BRM method (Huang et al., 2019): AFDexp 
  
  ################################################################################################
  ##				deltaSNP, TdeltaSNP, G and Gprime methods
  ################################################################################################  
  ## 2.1 QTLseqr (Mansfeld and Grumet, 2018): deltaSNP, TdeltaSNP, G and Gprime ----
  # deltaSNP and TdeltaSNP based on the method proposed in Takagi et al., 2013
  # G and Gprime according to the method proposed in Magwene et al., 2011
  # NOTE: simdata matrix is formatted according to QTLseqr input format
  
  simdata<-importFromTable(file = "simdata.csv",
                           highBulk = "HIGH",
                           lowBulk = "LOW",
                           chromList = c("chr"))
  
  simdata<-runQTLseqAnalysis(simdata,
                             windowSize = 300,
                             popStruc = "F2",
                             bulkSize = c(50,50),
                             replications = 1000,
                             depth = 50:100,
                             filter = 0.1,
                             intervals = c(90, 95, 99))
  
  simdata<- runGprimeAnalysis(simdata,
                              outlierFilter = "deltaSNP",
                              windowSize = 300,
                              filterThreshold = 0.1)

  
  ################################################################################################
  ##				EUCLIDEAN DISTANCE METHODS: EDm and ED100^4
  ################################################################################################
  ## 2.2 Euclidean Distance based statistics (Hill et al., 2003; Omboki et al., 2018; Zhang et al., 2019): EDm and ED100^4 ---- 
  # EDm for marker-based statistics
  # ED100^4 (=ED100_4) for fixed sliding windows of one hundred consecutive SNP markers 
  
  for (i in 1:nrow(simdata)){
    fqlow<-c((simdata$AD_ALT.LOW[[i]]/simdata$DP.LOW[[i]]),(simdata$AD_REF.LOW[[i]]/simdata$DP.LOW[[i]]))
    fqhigh<-c((simdata$AD_ALT.HIGH[[i]]/simdata$DP.HIGH[[i]]),(simdata$AD_REF.HIGH[[i]]/simdata$DP.HIGH[[i]]))
    simdata$EDm[i]<-dist(rbind(fqlow, fqhigh))
  }
  
  simdata$ED100<-0
  
  range<-100
  nsnp<-nrow(simdata)
  
  for (i in 1:(range/2)){ #First window or group of consecutive markers
    simdata$ED100[i]<-sum(simdata$EDm[1:(range)])
  }
  
  for (i in (range/2):((nsnp-(range/2)))) {
    simdata$ED100[i]<-sum(simdata$EDm[(i-(range/2)):(i+(range/2))])
  }
  
  for (i in (nsnp-(range/2)):(nsnp)){ #Last window or group of consecutive markers
    simdata$ED100[i]<-sum(simdata$EDm[(nsnp-(range/2)):nsnp])
  }
  
  simdata$ED100_4<-(simdata$ED100)^4 
}
  ################################################################################################
  ## 								QTGseq METHOD: LOD and Smooth LOD 
  ################################################################################################
  ## 2.3 QTG-Seq method (Zhang et al., 2019): LOD and Smooth LOD ----
  # LOD as single-marker based statistic
  # Smoothed LOD statistic (=SmLOD) for a fixed bandwidth equivalent to 3 Mbp  
  # Steps:
  ## a - Format and save the input data.
  ## b - Load and run the LOD function from Zhang et al. 2019 
  
  ## a - Format input data for QTG ----
  
  simdataQTG<-simdata[,c(1,2)]
  names(simdataQTG)[1] <- "Chromosome"
  names(simdataQTG)[2] <- "Pos"
  simdataQTG$Chromosome <- 1
  simdataQTG$"a in low pool"<-simdata$AD_ALT.LOW
  simdataQTG$"A in low pool"<-simdata$AD_REF.LOW
  simdataQTG$"a in high pool"<-simdata$AD_ALT.HIGH
  simdataQTG$"A in high pool"<-simdata$AD_REF.HIGH
  
  write.csv(simdataQTG, 
            file = "./simdataQTG.csv", 
            row.names = FALSE)

  ## b - Load and run the script LODFunction.R from Zhang et al. 2019 ----
  ## available at: https://github.com/caulilin/QTG_Seq/tree/master/SmoothLOD
  
  LOD<-function(dir= NULL,filegen = NULL,width = NULL,DrawPlot= NULL,chrom =  NULL,chromnum =  NULL,
                col= NULL,Plotformat1= NULL,Resolution= NULL){
    data1<-read.csv(paste(dir,"/",filegen,sep = ""),sep=",",header = TRUE,stringsAsFactors=FALSE)
    
    Calculate<-function(data,width){
      # columns are:
      # 1. chromosome
      # 2. chromosome coordinate 
      # 3. observed number of allele A in the pool with low phenotype
      # 4. observed number of allele a in the pool with low phenotype
      # 5. observed number of allele A in the pool with high phenotype
      # 6. observed number of allele a in the pool with high phenotype
      chr = data1[,1];cc<-unique(chr)
      ff<-numeric();tt<-numeric()
      for(iii in 1:length(cc)){
        data<-data1[which(data1[,1]==cc[iii]),]
        chrom = data[,1]; pos=data[,2];e = data[,3];f = data[,4];g = data[,5];h = data[,6]
        lfrq<-cbind(e,f);hfrq<-cbind(g,h)
        Al<-lfrq[,1];al<-lfrq[,2];Ah<-hfrq[,1];ah<-hfrq[,2];pos<-as.matrix(pos)
        P_QTGL<-Al/(Al+al);P_QTGH<-Ah/(Ah+ah)
        P_QTGL<-matrix(P_QTGL,ncol = 1);P_QTGH<-matrix(P_QTGH,ncol = 1)
        nn<-dim(P_QTGH)[1]
        nQTG<-cbind(Al,Ah)
        lod<-matrix(NA,nrow = nn,ncol = 1)
        for(i in 1:nn){
          c1<-Al+al;c2<-Ah+ah
          unexist_QTN<-(choose(c1[i],nQTG[i,1])*((1/2)^c1[i]))*(choose(c2[i],nQTG[i,2])*((1/2)^c2[i]))
          exist_QTN<-choose(c1[i],nQTG[i,1])*((P_QTGL[i])^nQTG[i,1])*((1-P_QTGL[i])^(c1[i]-nQTG[i,1]))*choose(c2[i],nQTG[i,2])*((P_QTGH[i])^nQTG[i,2])*((1-P_QTGH[i])^(c2[i]-nQTG[i,2]))
          lod[i]<-log10(exist_QTN/unexist_QTN)
        }
        result<-cbind(pos,lod)
        
        x=pos
        y=lod     
        w=width
        olen = length(x)
        #beginning intervals
        xfirst = x[1,]      #start site
        startband = x <= xfirst + w  
        xstart = xfirst - (x[startband] - xfirst)
        ystart = y[startband]
        nstart = length(xstart)-1
        #end of intervals
        xlast=x[length(x),]#last site
        endband = x >= xlast - w
        xend = xlast + (xlast - x[endband])
        yend = y[endband]
        nend = length(xend) - 1   
        a=rev(xstart)   #reversal
        a=a[-length(a)] #remove the last site of b
        b=rev(xend)     
        b=b[-1]         #remove the first site of b
        x<-c(a,x,b)
        c=rev(ystart)
        c=c[-length(c)] #remove the last site of b
        d=rev(yend)     
        d=d[-1]
        y<-c(c,y,d)
        lx = length(x);ly= length(y)
        
        cl<- makeCluster(6)
        registerDoParallel(cl)
        yw=foreach(i=1:lx)%dopar%
          {
            c = x[i]
            inband = ( x >= (c-w)&x <= (c+w))
            xfrac = (abs(x[inband] - c))/w
            xwt = (1.0 - abs(xfrac)^3)^3    #Weights 
            xwt[abs(xfrac) >= 1] = 0    
            ywin = sum(y[inband]*xwt)/sum(xwt)
          }
        stopCluster(cl)
        g2=cbind(y,yw)
        h1=length(a)+1
        h2=length(a)+olen
        g2=g2[h1:h2,]
        smooth_G=cbind(pos,g2)
        colnames(smooth_G)<-NULL
        G1<-as.matrix(smooth_G)
        G<-matrix(as.numeric(G1),nrow = nrow(G1))
        smoothG<-G[,3]*2*log(10)  
        smoothG<-as.matrix(smoothG,ncol=1)
        prob<-function(x){1-pchisq(x, 1)}
        p<-apply(smoothG,1,prob)
        p<-as.matrix(p,ncol=1)
        G<-cbind(chrom,G,p)
        colnames(G)<-c("Chrom","Pos","LOD","Smooth_LOD","P_Value")
        
        
        bb<-na.omit(G)
        bb1<-bb[,4]
        aa<-numeric()
        for (i in 1:(nrow(bb)-2)){
          if ((bb1[i+1]>bb1[i]) & (bb1[i+1]>bb1[i+2])){
            aa<-rbind(aa,bb1[i+1])
          }
        }
        bb2<-bb[match(aa,bb1),]
        loc<-bb2[order(bb2[,4],decreasing = TRUE),]  
        fresult<-loc[which(loc[,5]<=10^(-4)&loc[,3]>=3),]
        ff<-rbind(ff,G)
        tt<-rbind(tt,fresult)
      }
      if(dim(tt)[1]==0){
        warning("No Significant positions find")
        output<-list(ff=ff)
      }else{
        Significant<-seq(1:dim(tt)[1])
        fresult<-cbind(Significant,tt)
        colnames(fresult)<-c("QTG","Chrom","Pos","LOD","Smooth_LOD","P_Value")
        fresult1<-rbind(fresult,matrix("",nrow = dim(ff)[1]-dim(fresult)[1],ncol = 6))
        result<-cbind(ff,matrix("",nrow = dim(ff)[1],ncol = 1),fresult1)
        rm(G1,endband,g2,smooth_G,startband,yw,x,y,pos);gc()
        output<-list(result=result,ff=ff)
      }
      return(output)
    }
    calculateresult<-Calculate(data,width)
    if(length(names(calculateresult))==1){
      write.table(calculateresult$ff,paste(dir,"/","Final resultTESTSIM.csv",sep=""),sep=",",row.names=F,col.names = T) 
    }else{
      write.table(calculateresult$result,paste(dir,"/","Final resultTESTSIM.csv",sep=""),sep=",",row.names=F,col.names = T) 
    }
    tdf<-calculateresult$ff[,1:4]
    
    
    
    checkdata<-function(tdf){
      tdf<-as.data.frame(tdf)
      tdf<-na.omit(tdf)
      colnames(tdf)<-c("chr","bp","lod","smoothlod")
      x<-tdf
      if (!("chr" %in% names(x))) stop(paste("Column", chr, "not found!"))
      if (!("bp" %in% names(x))) stop(paste("Column", bp, "not found!"))
      if (!("lod" %in% names(x))) stop(paste("Column", p, "not found!"))
      ## warn if you don't have a snp column
      if (!("smoothlod" %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
      if (!is.numeric(x$chr)) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
      if (!is.numeric(x$bp)) stop(paste(bp, "column should be numeric."))
      if (!is.numeric(x$lod)) stop(paste(lod, "column should be numeric."))
      if (!is.numeric(x$smoothlod)) stop(paste(smoothlod, "column should be numeric."))
      return(x)
    }
    
    if(DrawPlot==TRUE){
      draw<-function(drawdata=NULL,col=NULL){
        x<-checkdata(drawdata)
        CHR=BP=P=index=NULL
        d=data.frame(CHR=x$chr, BP=x$bp, lod=x$lod,smoothlod=x$smoothlod)
        d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(lod) & is.numeric(smoothlod)))
        d <- d[order(d$CHR, d$BP), ]
        d$pos=NA
        d$index=NA
        ind = 0
        for (i in unique(d$CHR)){
          ind = ind + 1
          d[d$CHR==i,]$index = ind
        }
        
        nchr = length(unique(d$CHR))
        if (nchr==1) { ## For a single chromosome
          d$pos=d$BP/1000000
          xlabel = paste('Chromosome',unique(d$CHR),'(Mb)')
        } else { ## For multiple chromosomes
          lastbase=0
          ticks=NULL
          for (i in unique(d$index)) {
            if (i==1) {
              d[d$index==i, ]$pos=d[d$index==i, ]$BP
            } else {
              lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
              d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
            }
            ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
          }
          xlabel = 'Chromosome'
          labs <- unique(d$CHR)
        }
        
        
        # Initialize plot
        xmax = ceiling(max(d$pos) * 1.03)
        xmin = floor(max(d$pos) * -0.03)
        ymax<-(ceiling(max(d$lod)*1.25/10))*10
        margin_space<-0.5
        par(mar=c(3*margin_space,3*margin_space,margin_space,3*margin_space)+0.8*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
        suppressWarnings(def_args <- list(xaxt='n',yaxt="n",bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                                          xlim=c(xmin,xmax),ylim=c(0,ymax),
                                          xlab="", cex.lab=0.6,ylab=""))
        dotargs <- as.list(match.call())[-1L]
        dotargs <- list()
        ## And call the plot function passing NA, your arguments, and the default
        do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
        # Create a vector of alternatiting colors
        col=rep(col, max(d$CHR))
        
        
        
        # Add points to the plot
        if (nchr==1) {
          with(d,  points(pos, lod, col="LightGray", pch=20,cex=0.2))
        } else {
          # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
          icol=1
          for (i in unique(d$index)) {
            with(d[d$index==unique(d$index)[i], ], points(pos, lod, col="LightGray", pch=20,cex=0.1))
            icol=icol+1
          }
        }
        suppressWarnings(axis(4,ylim=c(0,ymax),cex.lab=0.5,mgp=c(3,-0.2,0),tcl=-0.2,cex.axis=0.4,col.axis = "DarkGray",col.ticks="DarkGray",col ="DarkGray"))
        mtext("LOD",side=4,line=0.5,cex=0.5,font = 1,col="DarkGray")
        
        
        
        par(new=TRUE,mar=c(3*margin_space,3*margin_space,margin_space,3*margin_space)+0.8*margin_space,mgp=c(0.7,0.2,0),cex.axis=0.4)
        suppressWarnings(def_args <- list(xaxt="n", bty='n', xaxs='i', yaxs='i', las=1, pch=20, yaxt="n" ,
                                          xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$smoothlod))),
                                          xlab=xlabel, ylab="",cex.lab=0.55,font = 1))#line=0.78,
        dotargs <- list()
        # And call the plot function passing NA, your arguments, and the default
        do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
        # Add an axis. 
        if (nchr==1) { #If single chromosome, ticks and labels automatic.
          suppressWarnings(axis(1,cex.lab=0.5,mgp=c(3,-0.15,0.3),tcl=-0.2,cex.axis=0.4))
        } else { # if multiple chrs, use the ticks and labels you created above.
          suppressWarnings(axis(1, at=ticks, labels=labs,cex.lab=0.5,mgp=c(3,-0.15,0.3),tcl=-0.2,cex.axis=0.4))
        }
        # Create a vector of alternatiting colors
        col=rep(col, max(d$CHR))
        # Add points to the plot
        if (nchr==1) {
          with(d, points(pos, smoothlod, pch=20, cex=0.1,col=col[1]))
        } else {
          # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
          icol=1
          for (i in unique(d$index)) {
            with(d[d$index==unique(d$index)[i], ], points(pos, smoothlod, col=col[icol],pch=20,cex=0.18))
            icol=icol+1
          }
        }
        axis(2,ylim=c(0,max(d$smoothlod)),cex.lab=0.5,mgp=c(3,0.2,0),tcl=-0.2,cex.axis=0.4,col.axis = col[1],col.ticks=col[1],col =col[1] )
        mtext("Smooth LOD",side=2,line=0.8,cex=0.5,font = 1,col = col[1])
      }
    }
    
    
    if(DrawPlot==TRUE){
      if(Resolution=="Low"){
        manwidth<-960;manhei<-600;manwordre<-20;manfigurere<-72
      }else if(Resolution=="High"){
        manwidth=15000; manhei=8000;units= "px";manwordre =30;manfigurere=1000
      }
      if(Plotformat1=="png"){
        png(paste(dir,"/","LOD.png",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat1=="tiff"){
        tiff(paste(dir,"/","LOD.tiff",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat1=="jpeg"){
        jpeg(paste(dir,"/","LOD.jpeg",sep=""),width=as.numeric(manwidth), height=as.numeric(manhei), units= "px", pointsize =as.numeric(manwordre),res=as.numeric(manfigurere))
      }else if(Plotformat1=="pdf"){
        
        pdf(paste(dir,"/","LOD.pdf",sep=""),width=15,fonts = "sans")
      }
      if("all"%in%chrom){
        draw(tdf,col)
        dev.off()
      }
      if("all"%in%chrom == FALSE){
        data1<-tdf
        data2<-NULL
        for (i in 1:length(unique(chrom))) {
          data2[[i]]<-data1[which(data1[,1]==chrom[i]),]
        }
        data2<-do.call(rbind,data2)
        draw(data2,col)
        dev.off()
      }
    }
  }
  
  LOD(dir="./",filegen = "simdataQTG.csv",width = 300,
      DrawPlot=FALSE,chrom ="1",col=c("blue", "red"),Plotformat1="tiff",Resolution="High")
  
  LODresult<-read.csv("./Final resultTESTSIM.csv", header=TRUE)
  
  ################################################################################################
  ## 								BRM method: AFDexp 
  ################################################################################################
  ## 2.4 BRM method (Huang et al., 2019): AFDexp ----
  # The method uses differences in allele frequency (ie. deltaSNP) as the marker-based statistic
  # Steps:
  ## a - Format and save the input data
  ## b - Calculate the uα/2 value for each recombination model. We use the code cal_ua_fk.R available at
  ##     https://github.com/huanglikun/BRM/tree/master/tools. The uα/2 values used were: 3.43 for Pearl
  ##     millet, 3.65 for rice and 3.71 for Foxtail millet 
  ## c - Create the links to files 'file_conf', 'file_chrlen' and 'file_bsa' 
  ## d - Then, we can run the code below adapted from Huang et al., 2019 BRM.R obtained from: https://github.com/huanglikun/BRM
  
  ## a - Format input data for BRM ----
  simdataBRM<-simdata[,c(1:2)]
  simdataBRM$CHROM<-"I"
  simdataBRM$AD_ALT.HIGH<-simdata$AD_ALT.HIGH
  simdataBRM$AD_REF.HIGH<-simdata$AD_REF.HIGH
  simdataBRM$AD_ALT.LOW<-simdata$AD_ALT.LOW
  simdataBRM$AD_REF.LOW<-simdata$AD_REF.LOW
  
  if((dir.exists(paste0("./","BRM"))==FALSE)){dir.create(paste0("./","BRM"))}
  base_dir <- getwd()
  write.table(simdataBRM, file = "./BRM/simdataBRM.bsa", row.names = FALSE , col.names = FALSE)
  
  #****************************
  ## b - Configuration file BRM_conf.txt where uα/2 value is created outside the program ----
  
  ## c - Create the links to files ----
  file_conf   <- "./BRM/BRM_conf.txt" # arguments configure file ###CHANGE THIS FILE FOR EACH MODEL CROP -ua###
  file_chrlen <- "./BRM/chr_length.tsv" # chromosome length file
  file_bsa    <- "./BRM/simdataBRM.bsa" # bsa file
  
  ## d - Script BRM.R from Huang et al., 2019 ----
  #########################################
  # functions
  #########################################
  #
  read_conf <- function(file_conf){
    conf <- list();
    #
    FH <- file(file_conf,"r");
    while(TRUE){
      rline <- readLines(FH,n=1);
      if (length(rline) == 0) break;
      if (length (grep ("^#", rline, perl=TRUE) ) > 0 ) next;
      if (length (grep ("^\\s", rline, perl=TRUE) ) > 0) next;
      if (nchar (rline) == 0) next;
      rline <- sub ("#.*", "", rline, perl=TRUE);
      rline <- sub ("\\s+$", "", rline, perl=TRUE);
      rline <- sub ("\\s*?=\\s*", "=", rline, perl=TRUE);
      arr <- strsplit (rline, "=")[[1]];
      conf[[ arr[1] ]] <- arr[2];
    }
    close(FH);
    return(conf);
  }
  
  #
  create_dir <- function(dir){
    if (file.exists(dir)){
      cat(dir," already exists.\n");
    }else{
      tryCatch({dir.create(dir,recursive=T)},warning=function(w){s <- as.character(w);stop(s)},error=function(e){s <- as.character(e);stop(s)});
    }
  }
  
  #
  chk_num_para <- function(x,name){
    if (is.null(x)){
      stop("Undefined ",name,".");
    }
    return(as.numeric(x));
  }
  
  #
  chk_dir <- function(path){
    dir <- dirname(path);
    if (dir != '.') create_dir(dir);
  }
  
  #
  chk_conf <- function(conf){
    if (is.null(conf$Design) || (conf$Design != "A" && conf$Design != "BH" && conf$Design != "BL")) stop("Incorrect Design.");
    if (is.null(conf$t) || (conf$t != 0 && conf$t != 1)) stop("Incorrect t.");
    conf$n1       <- chk_num_para(conf$n1,"n1");
    conf$n2       <- chk_num_para(conf$n2,"n2");
    conf$ua       <- chk_num_para(conf$ua,"ua");
    conf$UNIT     <- chk_num_para(conf$UNIT,"UNIT");
    conf$DEG      <- chk_num_para(conf$DEG,"DEG");
    conf$BLK      <- chk_num_para(conf$BLK,"BLK");
    conf$MIN      <- chk_num_para(conf$MIN,"MIN");
    conf$MINVALID <- chk_num_para(conf$MINVALID,"MINVALID");
    if (is.null(conf$Result1_File)) conf$Result1_File <- "result/result1.xls";
    if (is.null(conf$Result2_File)) conf$Result2_File <- "result/result2.xls";
    chk_dir(conf$Result1_File);
    chk_dir(conf$Result2_File);
    return(conf);
  }
  
  #####################
  # step 1
  #####################
  
  # block's middle position 
  def_block_pos <- function(chr, size){
    # block number
    n <- as.integer(2*chr/size)+2;
    if( n%%2 != 0 ) n <- n+1;
    n <- as.integer(n/2);
    # block index and the middle position of each block
    i <- c(1:n);
    pos <- (i-1)*size+floor(size/2); # middle position
    if(pos[n]>chr) pos[n]<-chr;
    return(pos);
  }
  
  # meta value of block
  cal_block_meta <- function(loc, val, chr, size, depth, MIN){
    # input: location vector, value vector
    # input: chr length, block size, location depth vector
    pos <- def_block_pos(chr, size);
    idx <- as.integer(0.5+loc/size)+1;
    #
    avg <- c();
    blockDepth <- c();
    for(i in 1:length(pos)){
      k  <- which(idx==i);
      no <- length(k);
      a <- NA;
      n <- 0;
      if (no > 0) {n <- sum(depth[k])};
      if( n >= MIN ) a <- sum(val[k])/n;
      avg <- c(avg, a);
    }
    return( list(pos=pos, avg=avg) );
  }
  
  check_variable <- function(A,B,C,D){
    idx <- c();
    # not missing
    idx <- c(idx, which(is.na(A)) );
    idx <- c(idx, which(is.na(B)) );
    idx <- c(idx, which(is.na(C)) );
    idx <- c(idx, which(is.na(D)) );
    #
    idx <- c(idx, which(!is.integer(A)) );
    idx <- c(idx, which(!is.integer(B)) );
    idx <- c(idx, which(!is.integer(C)) );
    idx <- c(idx, which(!is.integer(D)) );
    #
    idx <- c(idx, which(A+B==0) );
    idx <- c(idx, which(A+C==0) );
    idx <- c(idx, which(B+D==0) );
    idx <- c(idx, which(C+D==0) );
    #
    idx <- unique(idx);
    return(idx);
  }
  
  #
  run_step1 <- function(conf,file_chrlen,file_bsa,file_out,statistic){
    #
    # code from fANCOVO package
    # best smoothing parameter 
    opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
      as.crit <- function(x) {
        span   <- x$pars$span;
        traceL <- x$trace.hat;
        sigma2 <- sum(x$residuals^2)/(x$n - 1);
        aicc   <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2);
        gcv    <- x$n * sigma2/(x$n - traceL)^2;
        result <- list(span = span, aicc = aicc, gcv = gcv);
        return(result);
      }
      criterion <- match.arg(criterion);
      fn <- function(span) {
        mod <- update(model, span = span);
        as.crit(mod)[[criterion]];
      }
      result <- optimize(fn, span.range);
      return(list(span = result$minimum, criterion = result$objective));
    }
    
    # return
    fit_model <- list();
    # global options
    BLK <- as.numeric(conf$BLK) * as.numeric(conf$UNIT); # block size/step
    MIN <- as.numeric(conf$MIN); # least total depth in one block 
    DEG <- as.numeric(conf$DEG); # degree of polynomial
    MINVALID <- as.numeric(conf$MINVALID); # least valid blocks in one chromosome
    #
    if(MIN==0) MIN <- 1;
    if(MINVALID < 10) MINVALID <- 10;
    # read into memory
    tab <- read.table(file_bsa, as.is=TRUE);
    chr <- read.table(file_chrlen, as.is=TRUE);
    #
    dat_chr <- tab[,1];
    dat_pos <- tab[,2];
    A 		<- tab[,3];
    B 		<- tab[,4];
    C 		<- tab[,5];
    D 		<- tab[,6];
    ##
    idx <- check_variable(A, B, C, D);
    cat("Number of missing data:", length(idx), "\n");
    dat_all_chr <- dat_chr;
    dat_all_pos <- dat_pos;
    if(length(idx)>0){
      A <- A[-idx];
      B <- B[-idx];
      C <- C[-idx];
      D <- D[-idx];
      dat_all_chr <- dat_chr[-idx];
      dat_all_pos <- dat_pos[-idx];
    }
    ## 
    pool1Depth <- A + B;
    pool2Depth <- C + D;
    #
    cat("Number of imported chromosomes is", length(chr[,1]), ".\n");
    for(i in 1:length(chr[,1])){
      dat_chr <- dat_all_chr;
      dat_pos <- dat_all_pos;
      idx <- which(dat_chr == chr[i,1]);
      total <- length(idx);
      cat("Data size of", chr[i,1], "is", total, ".\n");
      if(total==0) next;
      #
      block_af1 <- cal_block_meta(dat_pos[idx], A[idx], chr[i, 2], BLK, pool1Depth[idx], MIN);
      block_af2 <- cal_block_meta(dat_pos[idx], C[idx], chr[i, 2], BLK, pool2Depth[idx], MIN);
      
      # initialize
      block <- list();
      #
      if(statistic=="AF1") {
        block <- block_af1;
      }else if (statistic=="AF2") {
        block <- block_af2;
      }else if (statistic=="AFD") {
        if (conf$Design == "BL"){
          block$avg <- block_af2$avg - block_af1$avg;
        }else{
          block$avg <- block_af1$avg - block_af2$avg;
        }
        block$pos <- block_af1$pos;
        block$num <- block_af1$num;
      }else if (statistic=="AAF"){
        block$avg <- (block_af1$avg + block_af2$avg) / 2;
        block$pos <- block_af1$pos;
        block$num <- block_af1$num;
      }else{
        stop("No statistic definded.");
      }
      #
      x <- as.numeric(block$pos);
      y <- as.numeric(block$avg);
      jdx <- which(!is.na(y));
      cat("Total blocks in ", chr[i,1], ": ", length(x), "\n", sep="");
      cat("Number of valid block:", length(jdx), "\n");
      # only consider those chromosomes that have at least MINVALID valid blocks
      if(length(jdx)<MINVALID) next;
      #
      fits <- list();
      fit0  <- loess(y[jdx]~x[jdx], degree=DEG);
      span1 <- opt.span(fit0, criterion="aicc")$span;
      fit1  <- loess(y[jdx]~x[jdx], degree=DEG, span=span1);
      fit_model[[ i ]] <- fit1;
      #
      plo <- predict(fit1, x, se=TRUE);
      value <- plo$fit;
      ###################################
      # manually correct fitted value
      if(statistic=="AFD")    value <- ifelse(value < -1,    -1, value);
      if(statistic=="AFD")    value <- ifelse(value >  1,     1, value);
      if(statistic=="AAF") value <- ifelse(value <  0,     0, value);
      if(statistic=="AAF") value <- ifelse(value >  1,     1, value);
      if(statistic=="AF1")     value <- ifelse(value <  0,     0, value);
      if(statistic=="AF1")     value <- ifelse(value >  1,     1, value);
      if(statistic=="AF2")     value <- ifelse(value <  0,     0, value);
      if(statistic=="AF2")     value <- ifelse(value >  1,     1, value);
      ####################################
      dat_chr <- rep(chr[i,1], length(block$pos));
      dat_pos <- block$pos;
      dat_avg <- round(block$avg,4);
      fit_avg <- round(value,4);
      out <- data.frame(dat_chr, dat_pos, dat_avg, fit_avg);
      if (i == 1){
        allout <- out;
      }else{
        allout <- rbind(allout,out);
      }
    }
    cat("\n");
    return(list(fits=fit_model,data=allout));
  }
  
  #
  run_step2 <- function(conf,data_aaf,data_af1,data_af2,file_out){
    #
    N1 <- as.numeric(conf$n1); # number of pool 1
    N2 <- as.numeric(conf$n2); # number of pool 2
    T  <- as.numeric(conf$t);  # level of population. For DH or RI etc., T=0; F2 or F3 etc., T=1
    ua <- as.numeric(conf$ua); # uα
    #
    dat_chr <- data_aaf[,1];
    dat_pos <- data_aaf[,2];
    dat_val <- data_aaf[,3];
    fit_val <- data_aaf[,4];
    fit_af1 <- data_af1[,4];
    fit_af2 <- data_af2[,4];
    # sample threshold
    dat_var <- (N1 + N2) / (2^T * N1 * N2) * fit_val * (1 - fit_val);
    dat_stdua <- ua * sqrt(dat_var);
    # theoretical threshold
    var_fix <- (N1 + N2) / (2^T * N1 * N2) * 0.25; # Fix AF = 0.5
    stdua_fix <- ua * sqrt(var_fix);
    # variance if this location is QTL
    var_qtl <- (fit_af1 * (1 - fit_af1) / (2^T * N1)) + (fit_af2 * (1 - fit_af2) / (2^T * N2));
    # 
    out <- data.frame(dat_chr, dat_pos, dat_val, fit_val, dat_stdua, dat_var, stdua_fix, var_qtl);
    return(list(theo=stdua_fix,data=out));
  }
  
  #
  return_cross_pos <- function(all_plots, x1, x2, y, tol){
    subregion  <- all_plots[which(all_plots[,1]>x1 & all_plots[,1]<x2),];
    y_distance <- abs(subregion[,2] - y);
    pos        <- subregion[which(y_distance==min(y_distance) & y_distance<tol),1];
  }
  
  #
  run_step3 <- function(conf,dat_afd,dat_var,file_out,fit_model_afd,fit_model_af1,fit_model_af2,threshold){
    # fixed global options
    tol <- 0.01 # tolerance threshold
    threshold <- as.numeric(threshold);
    #
    opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
      as.crit <- function(x) {
        span <- x$pars$span;
        traceL <- x$trace.hat;
        sigma2 <- sum(x$residuals^2)/(x$n - 1);
        aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2);
        gcv <- x$n * sigma2/(x$n - traceL)^2;
        result <- list(span = span, aicc = aicc, gcv = gcv);
        return(result);
      }
      criterion <- match.arg(criterion);
      fn <- function(span) {
        mod <- update(model, span = span);
        as.crit(mod)[[criterion]];
      }
      result <- optimize(fn, span.range);
      return(list(span = result$minimum, criterion = result$objective));
    }
    
    #
    locfit_by_loess <- function(x, y){
      # loess+AICc
      fit0  <- loess(y~x, degree = 2)
      span1 <- opt.span(fit0, criterion="aicc")$span
      fit1  <- loess(y~x, degree=2, span=span1)
    }
    
    #
    N1 <- as.numeric(conf$n1); # number of pool 1
    N2 <- as.numeric(conf$n2); # number of pool 2
    T  <- as.numeric(conf$t);  # level of population. For DH or RI etc., T=0; F2 or F3 etc., T=1
    # one-by-one for every chromosome
    chrom   <- c("#Chr.");
    peak_x  <- c("Pos.");
    peak_y  <- c("Val.");
    type    <- c("Peak Dir.");
    left_x  <- c("Start");
    right_x <- c("End");
    # title write to file
    out <- data.frame(chrom, peak_x, peak_y, type, left_x, right_x);
    write.table(out, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=FALSE);
    #
    dat_chr <- dat_afd[,1];
    dat_pos <- dat_afd[,2];
    dat_avg <- dat_afd[,3];
    #
    dat_std <- sqrt(dat_var[,6]);
    qtl_std <- sqrt(dat_var[,8]);
    #
    chr     <- data.frame(unique(dat_chr));
    #
    for(i in 1:length(chr[,1])){
      thischr <- chr[i,1];
      idx <- which(dat_chr == chr[i,1] & !is.na(dat_avg));
      total <- length(idx);
      cat("data size of", chr[i,1], "is", total, "\n");
      if(total==0) next;
      #
      x <- dat_pos[idx];
      y <- dat_avg[idx];
      std  <- dat_std[idx];
      qstd <- qtl_std[idx];
      fit1 <- locfit_by_loess(x, y);
      #
      newx <- seq(x[1],x[length(x)],1);
      newx_len <- length(newx);
      plo  <- data.frame(x=newx, y=predict(fit1, newx));
      peakxs   <- newx[which(diff(diff(plo$y)>0)!=0L)+1L];
      len  <- length(peakxs);
      # 
      if (len==0){
        k  <- ifelse(diff(c(plo$y[1],plo$y[length(plo)]))>0L,1L,-1L);
        ifelse(k>0L,k <- c(1L,-1L),k <- c(-1L,1L));
        pks <- data.frame(peak_x=c(plo$x[1],plo$x[length(plo$x)]),peak_y=c(plo$y[1],plo$y[length(plo$y)]),type=k);
        peak_std   <- c(qstd[1],qstd[length(qstd)]);
      }else{
        k    <- diff(diff(plo$y)>0L);
        k    <- k[which(k!=0L)];
        stat     <- k[1]>0;
        stat0    <- stat;
        
        peakxs   <- (peakxs[1:len-1]+peakxs[2:len])/2;
        
        pks  <- sapply(as.data.frame(rbind(c(newx[1],peakxs), c(peakxs,newx[length(newx)]))),
                       function(k){ 
                         stat <<- !stat;
                         optimize(f=function(x){
                           predict(fit1, newdata=data.frame(x))}, 
                           maximum=stat, interval=k)});
        pks  <- as.data.frame(t(pks));
        names(pks) <- c("peak_x","peak_y");
        pks$peak_x <- as.numeric(pks$peak_x);
        pks$peak_y <- as.numeric(pks$peak_y);
        pks$type   <- k;
        peak_af1   <- predict(fit_model_af1[[ i ]],pks$peak_x);
        peak_af2   <- predict(fit_model_af2[[ i ]],pks$peak_x);
        peak_std   <- as.numeric(sqrt(peak_af1 * (1 - peak_af1) / (2^T * N1) + peak_af2 * (1 - peak_af2) / (2^T * N2)));
        # add the first and the last position.
        firstPos   <- c(plo$x[1],plo$y[1],-1*pks$type[1]);
        lastPos    <- c(plo$x[length(plo$x)],plo$y[length(plo$y)],-1*pks$type[length(pks$type)]);
        pks        <- rbind(firstPos, pks, lastPos);
        peak_std   <- c(qstd[1],peak_std,qstd[length(qstd)]);
        k          <- pks$type;
      }
      all_peak_plots <- pks$peak_x;
      all_ident_y    <- pks$peak_y + k*1.65*peak_std;
      all_x1         <- all_peak_plots[1:(length(all_peak_plots)-1)];
      all_x3         <- floor(all_peak_plots[1:(length(all_peak_plots)-1)]);
      all_x2         <- all_peak_plots[2:(length(all_peak_plots))];
      all_x4         <- all_peak_plots[2:length(all_peak_plots)];
      pos1_list <- c(x[1]);
      pos2_list <- c();
      all_x1_len<- length(all_x1);
      for (i in 1:all_x1_len){
        pos1     <- return_cross_pos(plo, all_x1[i], all_x2[i], all_ident_y[i+1], tol);
        j        <- i ;
        while(j > 1 & length(pos1) == 0){
          j    <- j - 1;
          pos1 <- return_cross_pos(plo, all_x1[j], all_x2[i], all_ident_y[i+1], tol);
        }
        pos2     <- return_cross_pos(plo, all_x3[i], all_x4[i], all_ident_y[i], tol);
        j <- i;
        while(j < all_x1_len & length(pos2) == 0){
          j <- j + 1;
          pos2 <- return_cross_pos(plo, all_x3[i], all_x4[j], all_ident_y[i], tol);
        }
        if (length(pos1) == 0) {pos1 <- c(all_x1[1]);respos1 <- pos1;}
        if (length(pos2) == 0) {pos2 <- c(all_x4[length(all_x4)]);respos2 <- pos2;}
        pos1_list   <- c(pos1_list,pos1);
        pos2_list   <- c(pos2_list,pos2);
      }
      pos2_list     <- c(pos2_list,x[length(x)]);
      pks$type      <- ifelse(k==1,"-","+");
      pks$peak_x    <- trunc(pks$peak_x);
      pks$peak_y    <- round(pks$peak_y,4);
      pks <-data.frame(thischr, pks,left_x=c(pos1_list), right_x=c(pos2_list));
      #
      con1 <- which((pks$peak_y >= threshold | pks$peak_y <= -threshold) & pks$peak_y * k < 0);
      pks  <- pks[con1,];
      con2 <- c();
      for (i in 1:length(pks$left_x)){
        ileft  <- pks$left_x[i];
        iright <- pks$right_x[i];
        ifex   <- which((pks$left_x >= ileft & pks$right_x < iright) | (pks$left_x > ileft & pks$right_x <= iright));
        if(length(ifex) > 0) con2 <- c(con2,i);
      }
      if(length(con2) > 0) pks <- pks[-con2,];
      write.table(pks, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
    }
    
  }
  
  #########################################
  # main
  #########################################
  ptm <- proc.time(); # begin to record runtime 
  #
  # read configure file
  cat("Read global parameters from", file_conf, ".\n");
  conf <- read_conf(file_conf);
  #conf <- read_conf(file="F:/BSA_simdata_Paper/BRM-Huang2019/BRM-master/configureExample/designA/BRM_conf.txt");
  # check parameter
  conf <- chk_conf(conf);
  # print parameter
  for(name in names(conf) ){
    cat("Parameter ", name, "=", conf[[name]],  "\n", sep="");
  }
  
  cat("\n");
  # step 1
  af1s <- run_step1(conf,file_chrlen,file_bsa,conf$AF1_File,"AF1");
  af2s <- run_step1(conf,file_chrlen,file_bsa,conf$AF2_File,"AF2");
  afds <- run_step1(conf,file_chrlen,file_bsa,conf$AFD_File,"AFD");
  if (conf$Design == "A"){
    aafs <- run_step1(conf,file_chrlen,file_bsa,conf$AAF_File,"AAF");
  }else if (conf$Design == "BH" || conf$Design == "BL"){
    aafs <- af2s;
    file.copy(conf$AF2_File,conf$AAF_File);
  }else{
    stop("Is design not A or B?");
  }
  
  af1_model <- af1s$fits;
  af1_data  <- af1s$data;
  af2_model <- af2s$fits;
  af2_data  <- af2s$data;
  aaf_model <- aafs$fits;
  aaf_data  <- aafs$data;
  afd_model <- afds$fits;
  afd_data  <- afds$data;
  
  # step 2
  step2s <- run_step2(conf,aaf_data,af1_data,af2_data,conf$STEP2OUT);
  theoretical_threshold <- round(step2s$theo,4);
  var_data              <- step2s$data;
  cat("Theoretical threshold is ±",theoretical_threshold,".\n");
  
  # threshold and title
  cat("Export data to", conf$Result1_File, "\n");
  thr_title <- "#Theoretical threshold is ±";
  out <- data.frame(thr_title, theoretical_threshold);
  write.table(out, conf$Result1_File, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=FALSE);
  
  dat_chr  <- c("#Chr.");
  dat_pos  <- c("Pos.");
  dat_avg1 <- c("AF1-Observed");
  fit_avg1 <- c("AF1-Expected");
  dat_avg2 <- c("AF2-Observed");
  fit_avg2 <- c("AF2-Expected");
  dat_avg3 <- c("AFD-Observed");
  fit_avg3 <- c("AFD-Expected");
  dat_avg4 <- c("AFP-Observed");
  fit_avg4 <- c("AFP-Expected");
  dat_varua <- c("Sample threshold");
  out <- data.frame(dat_chr, dat_pos, dat_avg1, fit_avg1, dat_avg2, fit_avg2, dat_avg3, fit_avg3, dat_avg4, fit_avg4, dat_varua);
  write.table(out, conf$Result1_File, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
  
  cat("Writing\n");
  out <- data.frame(af1_data[,1], af1_data[,2], af1_data[,3], af1_data[,4], af2_data[,3], af2_data[,4], afd_data[,3], afd_data[,4], aaf_data[,3], aaf_data[,4], round(var_data[,5],4));
  write.table(out, conf$Result1_File, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
  cat("\n");
  
  # step 3
  run_step3(conf,afd_data,var_data,conf$Result2_File,afd_model,af1_model,af2_model,theoretical_threshold);
  
  # runtime
  proc.time() - ptm; 
  #
  cat("finished.\n")
  
  #--------------------------
  
  result1<-read.table("./BRM/result_H_L/result1.xls")
  result2<-read.table("./BRM/result_H_L/result2.xls")
  
# Results for the 9 statistical methods ----
  ResultBSA<-simdata[,c(1:12,14,20,21,25,27)]
  ResultBSA$LOD<-LODresult$LOD
  ResultBSA$SmLOD<-LODresult$Smooth_LOD
  
  ResultBRM<-result1[,c(1,2,8)]
  names(ResultBRM)[1]<-"CHROM"
  names(ResultBRM)[2]<-"POS"
  names(ResultBRM)[3]<-"AFDexp"

}
## Plot results
  plot(ResultBSA$deltaSNP~ResultBSA$POS)
  plot(ResultBSA$tricubeDeltaSNP~ResultBSA$POS)
  plot(ResultBSA$G~ResultBSA$POS)
  plot(ResultBSA$Gprime~ResultBSA$POS)  
  plot(ResultBSA$EDm~ResultBSA$POS)
  plot(ResultBSA$ED100_4~ResultBSA$POS) 
  plot(ResultBSA$LOD~ResultBSA$POS) 
  plot(ResultBSA$SmLOD~ResultBSA$POS) 
  plot(ResultBRM$AFDexp~ResultBRM$POS) 
  