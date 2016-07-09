calc.imp.accuracy <- function(truegeno,imputedgeno,mapHD,mapLD,nIID,format,accuracy_type,missgeno){
  ## Date : 24:09:2015
  options(warn=-1)
  #==================================================================================================#
  # Truegeno == filename of true genotype (PLINK or allele dosage or gene content format)
  # Imputedgeno == filename of true genotype (PLINK or allele dosage or gene contetn format)
  # note that both imputedgeno and truegeno should have the same file format
  # mapHD == filename of the HIGH density SNP map file
  # mapLD == filename of the LOW density or imputed genotype SNP map file
  # nIID == number of animals genotyped, you can specify higher number but not less
  # format == The format of the genotype data (either PLINK or gene content/allele dosage format)
  # accuracy_type == The method to use in calculating accuracy (1. simple 2. calus_mulder 3. both)
  # Read: Evaluation of measures of correctness of genotype imputation in the context of genomic prediction: 
  # a review of livestock applications (Calus et al. 2014 - http://dx.doi.org/10.1017/S1751731114001803 )
  # missgeno == the value of the missing genotype in the data file (eg. NA, 5, -9 ...)
  #==================================================================================================#
  
  if (missing(truegeno))
    stop("  Need to specify the TRUE genotype file  ")
  if (missing(imputedgeno))
    stop("  Need to specify the IMPUTED genotype file  ")
  if (missing(mapHD))
    stop("  Specify the map file of the high density or true genotype file ")
  if (missing(mapLD))
    stop("  Specify the map file of the low density or imputed genotype file  ")
  if (missing(nIID))
    stop(" Specify the minimum number of animals in the dataset")
  if (missing(format))
    stop(cat(" Specify a format for your datasets \n Three (3) types of fomat are allowed \n 1. PLINK 'ped' format \n 2. Genecontent or \n 3. allele dosage format"))
  if (missing(accuracy_type))
    stop(cat("Specify the method of computing the accuracies \n Three (3) methods are described \n 1. 'Simple' imputation without centering and scaling genotypes \n 2. 'calus_mulder' this requires centering and scaling genotype file \n 3. 'both' using both methods "))
  if (missing(missgeno))
    stop(" Specify what value or character is used as missing genotype ")
  
  ####################################################################################
  ##### Importing True genotypes ###########
  if(format[1]=="genotype" | format[1]=="dosage"){
    ####### Dealing with True genotypes
    ####### importing and initial editing
    del <- read.table(truegeno,header=F,nrows=500)
    classes <- sapply(del,class)
    rm(del)
    cat('....... Prelim true genotypes.......\n')      
    Xtrue <- read.table(truegeno,header=F,nrows=nIID,colClasses=classes,na.strings=missgeno[1])
    #Xtrue[,2:ncol(Xtrue)] <-sapply(Xtrue[,2:ncol(Xtrue)],function(x){ifelse(x==missgeno,NA,x)})
    cat('....... True genotypes have been imported.......\n')
    cat(' \n')
  }
  #### using tped and tfam format == NOTE this format is very similar to BEAGLE v3 output format
  #### dummy columns can be created for a beagle output file to be used
  else if(format[1]=="tped"){
    ####### Dealing with True genotypes
    ####### importing and initial editing
    Xtrue <- read.table(paste(truegeno,'.tped',sep=''),header=F,na.strings=missgeno[1])
    Xtrue <- (Xtrue[,-1:-4])-1
    seq1 <- seq(1,ncol(Xtrue),2)
    seq2 <- seq(2,ncol(Xtrue),2)
    Xtrue <- Xtrue[,seq1] + Xtrue[,seq2]
    Xtrue <- t(Xtrue)
    Xfam <- read.table(paste(truegeno,'.tfam',sep=''),header=F)
    Xtrue <- cbind.data.frame(data.frame(Xfam[,2]),Xtrue)
    cat('....... True genotypes have been imported.......\n')
    cat(' \n')
  }
  
  #******************* using PLINK file ped and map format 
  else if(format[1]=="plink"){
    ####### Dealing with True genotypes
    ####### importing and initial editing
    del <- read.table(truegeno,header=F,nrows=800)
    classes <- sapply(del,class)
    rm(del)
    cat('....... Prelim true genotypes.......\n')      
    Xtrue <- read.table(truegeno,header=F,nrow=nIID,na.strings=missgeno[1])
    geno <- (Xtrue[,-1:-6]-1) 
    geno <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]
    Xtrue <- cbind(Xtrue[,2],geno)
    cat('....... True genotypes have been imported.......\n')
    cat(' \n')
  } 
  
  ####################################################################################
  ##### Importing Imputed genotypes ###########
  if(format[2]=="genotype" | format[2]=="dosage"){
    ####### Dealing with imputed genotypes
    ####### importing and initial editing
    del <- read.table(imputedgeno,header=F,nrows=800)
    classes <- sapply(del,class)
    rm(del)
    cat('....... Prelim imputed genotypes.......\n')
    Ximp <- read.table(imputedgeno,header=F,nrows=nIID,colClasses=classes,na.strings=missgeno[2])
    cat('...... Imputed genotypes have been imported .......\n')
    cat(' \n')
  }
  #### using tped and tfam format == NOTE this format is very similar to BEAGLE v3 output format
  #### dummy columns can be created for a beagle output file to be used
  else if(format[2]=="tped"){
    ####### Dealing with imputed genotypes
    ####### importing and initial editing
    Ximp <- read.table(paste(imputedgeno,'.tped',sep=''),header=F,na.strings=missgeno[2])
    Ximp <- (Ximp[,-1:-4])-1
    seq1 <- seq(1,ncol(Ximp),2)
    seq2 <- seq(2,ncol(Ximp),2)
    Ximp <- Ximp[,seq1] + Ximp[,seq2]
    Ximp <- t(Ximp)
    XfamI <- read.table(paste(imputedgeno,'.tfam',sep=''),header=F)
    Ximp <- cbind.data.frame(data.frame(XfamI[,2]),Ximp)
    cat('....... imputed genotypes have been imported.......\n')
    cat(' \n')
  }
  
  #******************* using PLINK file ped and map format 
  else if(format[2]=="plink"){
    ####### Dealing with imputed genotypes
    ####### importing and initial editing
    del <- read.table(imputedgeno,header=F,nrows=800)
    classes <- sapply(del,class)
    rm(del)
    cat('....... Prelim imputed genotypes.......\n')      
    Ximp <- read.table(imputedgeno,header=F,nrow=nIID,na.strings=missgeno[2])
    geno <- (Ximp[,-1:-6]-1) 
    geno <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]
    Ximp <- cbind(Ximp[,2],geno)
    rm(geno)      
    cat('...... Imputed genotypes have been imported .......\n')
    cat(' \n')
  } 
  
  #****     importing MAP files     *******#
  # higher density map file
  mapHD <- read.table(mapHD,header=F)
  mapHD$num <- 1:nrow(mapHD)
  cat('...... Higher density map file have been imported .......\n')
  
  # lower density map file
  mapLD <- read.table(mapLD,header=F)
  cat('...... lower density have been imported .......\n')
  
  ### Further editing and preparation of data file 
  idstrue <- data.frame(IID=Xtrue[,1])
  idsimp <- data.frame(IID=Ximp[,1])
  idscommon <- merge(idstrue,idsimp,by=1)
  if(nrow(idscommon)==0){stop("No true genotypes of imputed animals \n")}
  if(nrow(idscommon)==1){cat("true genotypes of ONLY one imputed animal \n SNP accuracies would be NA's \n")}
  
  if(accuracy_type=="simple"){
    cat ('*******************************************************************************************************************\n')
    cat ('*      Imputation acuracy will be computed with the SIMPLE ---- (GENOTYPES ARE NOT CENTRED and SCALED)            *\n')
    cat ('*ALTERNATIVE method requiring centering and scaling genotypes is discussed by Calus et al. 2014 in Animal Journal *\n')
    cat ('*  Paper title "Evaluation of measures of correctness of genotype imputation in the context of genomic prediction *\n')
    cat ('*                                 http://dx.doi.org/10.1017/S1751731114001803                                     *\n')
    cat ('*******************************************************************************************************************\n')
    
    Xtrue <- merge(idscommon,Xtrue,by=1,sort=F)[,-1]
    Ximp <- merge(idscommon,Ximp,by=1,sort=F)[,-1]
    nsnps = ncol(Xtrue)
    anim = nrow(Xtrue)
    iterchecks.snps <- round(nsnps/10,digits=-1)
    iterchecks.anim <- round(anim/10,digits=0)
    
    cat(' \n')
    ##### SNP imputation accuracy
    cat('....... Computation of SNP-specific accuracies started .......\n')
    cat(' \n')
    if (format[2]=='dosage'){
      Accuracy.snps <- matrix(0,nsnps,2)
      colnames(Accuracy.snps) <- c('correl','sqrdcorrel')
      for (j in 1:nsnps){
        datbind <- cbind(Ximp[,j],Xtrue[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])==0){
          dat <- NA
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat^2
        }
        if(j %% iterchecks.snps==0){
          cat(paste('SNP specific accuracy now at ...',j,' ...out of a total of ...',nsnps,sep=''),' \n')
        }
      }
    }
    else {
      Accuracy.snps <- matrix(0,nsnps,4)
      colnames(Accuracy.snps) <- c('correl','sqrdcorrel','PERC',"AError")
      for (j in 1:nsnps){
        datbind <- cbind(Ximp[,j],Xtrue[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])==0){
          dat <- NA
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat^2
        }
        
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.snps[j,3] <- datperc
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.snps[j,4] <- 1-datAER
        
        if(j %% iterchecks.snps==0){
          cat(paste('SNP specific accuracy now at ...',j,' ...out of a total of ...',nsnps,sep=''),' \n')
        }
      }
    }
    Accuracy.snps <- round(Accuracy.snps,digits=4)
    Accuracy.snps <- cbind(mapHD[,c(1,2,4)],Accuracy.snps)
    colnames(Accuracy.snps)[1:3] <- c('CHR','SNPname','Position')
    cat('.... Computation of SNP-specific accuracies finished ....... \n')
    
    #******************************* Animal specific accuracy for all SNPS
    cat(' \n')
    cat('.... Computation of Animal-specific accuracies with ALL SNPs started .......\n')
    cat(' \n')
    Xtrue_T <- t(Xtrue)
    Ximp_T <- t(Ximp)
    
    if(format[2]=='dosage'){
      Accuracy.animALLSNP <- matrix(0,anim,2)
      colnames(Accuracy.animALLSNP) <- c('correl','sqrdcorrel')
      for (l in 1:anim){
        datbind <- cbind(Ximp_T[,l],Xtrue_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,1] <- dat
        Accuracy.animALLSNP[l,2] <- dat^2
        if(l %% iterchecks.anim==0){
          cat(paste('Sample specific accuracy now at ...',l,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    else {
      Accuracy.animALLSNP <- matrix(0,anim,4)
      colnames(Accuracy.animALLSNP) <- c('correl','sqrdcorrel','PERC','AError')
      for (l in 1:anim){
        datbind <- cbind(Ximp_T[,l],Xtrue_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,1] <- dat
        Accuracy.animALLSNP[l,2] <- dat^2
        
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.animALLSNP[l,3] <- datperc
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.animALLSNP[l,4] <- 1-datAER
        
        if(l %% iterchecks.anim==0){
          cat(paste('Sample specific accuracy now at ...',l,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    
    Accuracy.animALLSNP <- round(Accuracy.animALLSNP,digits=4)
    Accuracy.animALLSNP <- cbind.data.frame(idscommon,Accuracy.animALLSNP)
    assign('Accuracy.animALLSNP',Accuracy.animALLSNP)
    cat('.... Computation of Sample-specific accuracies with ALL SNPs finished .......\n')
    
    
    #******** Animal specific accuracy for only imputed SNPS
    map <- merge(mapHD,mapLD,by=2,sort=F)
    map <- sort(map$num)
    Xtrue_T <- Xtrue_T[-map,]
    Ximp_T <- Ximp_T[-map,]
    nsnps = nrow(Xtrue_T)
    
    #********* Sample specific accuracy
    cat(' \n')
    cat('.... Computation of Sample-specific accuracies with ONLY IMPUTED SNPs started .......\n')
    cat(' \n')
    if(format[2]=='dosage'){
      Accuracy.animIMPUTEDSNP <- matrix(0,anim,2)
      colnames(Accuracy.animIMPUTEDSNP) <- c('correl','sqrdcorrel')
      for (m in 1:anim){
        datbind <- cbind(Ximp_T[,m],Xtrue_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,1] <- dat
        Accuracy.animIMPUTEDSNP[m,2] <- dat^2
        if(m %% iterchecks.anim==0){
          cat(paste('Sample specific accuracy now at ...',m,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    else {
      Accuracy.animIMPUTEDSNP <- matrix(0,anim,4)
      colnames(Accuracy.animIMPUTEDSNP) <- c('correl','sqrdcorrel','PERC','AError')
      for (m in 1:anim){
        datbind <- cbind(Ximp_T[,m],Xtrue_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,1] <- dat
        Accuracy.animIMPUTEDSNP[m,2] <- dat^2
        
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.animIMPUTEDSNP[m,3] <- datperc
        
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.animIMPUTEDSNP[m,4] <- 1-datAER
        
        if(m %% iterchecks.anim==0){
          cat(paste('Sample specific accuracy now at ...',m,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    Accuracy.animIMPUTEDSNP <- round(Accuracy.animIMPUTEDSNP,digits=4)
    Accuracy.animIMPUTEDSNP <- cbind.data.frame(idscommon,Accuracy.animIMPUTEDSNP)
    cat('.... Computation of Sample-specific accuracies with ONLY IMPUTED SNPs finished .......\n')
    
    Accuracy <- list(SNP_specific=Accuracy.snps,
                     Sample_specific.ALLSNP=Accuracy.animALLSNP,
                     Sample_specific.impSNP=Accuracy.animIMPUTEDSNP)
    
    cat(' \n')
    cat('.... Four types of accuracies are computed for both SNPs and animals/samples .......\n')
    cat('  1. correlation between true and imputed SNPs............. correl \n')
    cat('  2. squared correlation between true and imputed SNPs..... sqrdcorrel \n')
    cat('  3. Percentage (%) of correctly imputed SNPs.............. PERC \n')
    cat('  4. Allelic error rate (1-allelic error rate) ............ AER \n')
    cat(' ')
    return(Accuracy)
  } 
  
  else if (accuracy_type=="calus_mulder"){
    cat ('*******************************************************************************************************************\n')
    cat ('*  Imputation acuracy will be computed with the recently discussed method of calus et al. 2014 in Animal Journal  *\n')
    cat ('*  ...................        That means GENOTYPES ARE CENTRED and SCALED                   ..................    *\n') 
    cat ('*  Paper title "Evaluation of measures of correctness of genotype imputation in the context of genomic prediction *\n')
    cat ('*                                 http://dx.doi.org/10.1017/S1751731114001803                                     *\n')
    cat ('*******************************************************************************************************************\n')

    Xtrue <- merge(idscommon,Xtrue,by=1,sort=F)[,-1]
    Ximp <- merge(idscommon,Ximp,by=1,sort=F)[,-1]
    Xtrue.calusmulder <- scale(x=Xtrue,center=T,scale=T)
    Ximp.calusmulder <- scale(x=Ximp,center=T,scale=T)
    nsnps = ncol(Xtrue)
    anim = nrow(Xtrue)
    iterchecks.snps <- round(nsnps/10,digits=-1)
    iterchecks.anim <- round(anim/10,digits=0)
    
    cat(' \n')
    ##### SNP imputation accuracy
    cat('....... Computation of SNP-specific accuracies started .......\n')
    cat(' \n')
    if(format[2]=='dosage'){
      Accuracy.snps <- matrix(0,nsnps,2)
      colnames(Accuracy.snps) <- c('correl_CM','sqrdcorrel_CM')
      for (j in 1:nsnps){
        datbind <- cbind(Ximp.calusmulder[,j],Xtrue.calusmulder[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])==0){
          dat <- NA
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat^2
        }
        if(j %% iterchecks.snps==0){
          cat(paste('SNP specific accuracy now at ...',j,' ...out of a total of ...',nsnps,sep=''),' \n')
        }
      }
    } 
    else {
      Accuracy.snps <- matrix(0,nsnps,4)
      colnames(Accuracy.snps) <- c('correl_CM','sqrdcorrel_CM','PERC',"AError")
      for (j in 1:nsnps){
        datbind <- cbind(Ximp.calusmulder[,j],Xtrue.calusmulder[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])==0){
          dat <- NA
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat^2
        }
        
        datbind <- cbind(Ximp[,j],Xtrue[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- datbind[,1]-datbind[,2]
        
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.snps[j,3] <- datperc
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.snps[j,4] <- 1-datAER
        
        if(j %% iterchecks.snps==0){
          cat(paste('SNP specific accuracy now at ...',j,' ...out of a total of ...',nsnps,sep=''),' \n')
        }
      }
    }   
    Accuracy.snps <- round(Accuracy.snps,digits=4)
    Accuracy.snps <- cbind(mapHD[,c(1,2,4)],Accuracy.snps)
    colnames(Accuracy.snps)[1:3] <- c('CHR','SNPname','Position')
    cat('.... Computation of SNP-specific accuracies finished ....... \n')
    
    
    ##########************************* Sample specific accuracy for all SNPS
    cat(' \n')
    cat('.... Computation of Sample-specific accuracies with ALL SNPs started .......\n')
    cat(' \n')
    Xtrue_T <- t(Xtrue)
    Ximp_T <- t(Ximp)
    Xtrue.calusmulder_T <- t(Xtrue.calusmulder)
    Ximp.calusmulder_T <- t(Ximp.calusmulder)
    
    if(format[2]=='dosage'){
      Accuracy.animALLSNP <- matrix(0,anim,2)
      colnames(Accuracy.animALLSNP) <- c('correl_CM','sqrdcorrel_CM')
      for (l in 1:anim){
        datbind <- cbind(Ximp.calusmulder_T[,l],Xtrue.calusmulder_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,1] <- dat
        Accuracy.animALLSNP[l,2] <- dat^2  
        if(l %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',l,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    } 
    else {
      Accuracy.animALLSNP <- matrix(0,anim,4)
      colnames(Accuracy.animALLSNP) <- c('correl_CM','sqrdcorrel_CM','PERC','AError')
      for (l in 1:anim){
        datbind <- cbind(Ximp.calusmulder_T[,l],Xtrue.calusmulder_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,1] <- dat
        Accuracy.animALLSNP[l,2] <- dat^2
        
        datbind <- cbind(Ximp_T[,l],Xtrue_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.animALLSNP[l,3] <- datperc
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.animALLSNP[l,4] <- 1-datAER
        
        if(l %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',l,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    Accuracy.animALLSNP <- round(Accuracy.animALLSNP,digits=4)
    Accuracy.animALLSNP <- cbind.data.frame(idscommon,Accuracy.animALLSNP)
    assign('Accuracy.animALLSNP',Accuracy.animALLSNP)
    cat('.... Computation of Animal-specific accuracies with ALL SNPs finished .......\n')
    
    
    #*************************** Animal specific accuracy for only imputed SNPS
    map <- merge(mapHD,mapLD,by=2,sort=F)
    map <- sort(map$num)
    Xtrue_T <- Xtrue_T[-map,]
    Ximp_T <- Ximp_T[-map,]
    nsnps = nrow(Xtrue_T)
    Xtrue.calusmulder_T <- Xtrue.calusmulder_T[-map,]
    Ximp.calusmulder_T <- Ximp.calusmulder_T[-map,]
    
    #*************************** Animal specific accuracy
    cat(' \n')
    cat('.... Computation of Animal-specific accuracies with ONLY IMPUTED SNPs started .......\n')
    cat(' \n')
    if(format[2]=='dosage'){
      Accuracy.animIMPUTEDSNP <- matrix(0,anim,2)
      colnames(Accuracy.animIMPUTEDSNP) <- c('correl_CM','sqrdcorrel_CM')
      for (m in 1:anim){
        datbind <- cbind(Ximp.calusmulder_T[,m],Xtrue.calusmulder_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,1] <- dat
        Accuracy.animIMPUTEDSNP[m,2] <- dat^2
        if(m %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',m,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    else {
      Accuracy.animIMPUTEDSNP <- matrix(0,anim,4)
      colnames(Accuracy.animIMPUTEDSNP) <- c('correl_CM','sqrdcorrel_CM','PERC','AError')
      for (m in 1:anim){
        datbind <- cbind(Ximp.calusmulder_T[,m],Xtrue.calusmulder_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,1] <- dat
        Accuracy.animIMPUTEDSNP[m,2] <- dat^2
        
        datbind <- cbind(Ximp_T[,m],Xtrue_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.animIMPUTEDSNP[m,3] <- datperc
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.animIMPUTEDSNP[m,4] <- 1-datAER
        
        if(m %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',m,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    
    Accuracy.animIMPUTEDSNP <- round(Accuracy.animIMPUTEDSNP,digits=4)
    Accuracy.animIMPUTEDSNP <- cbind.data.frame(idscommon,Accuracy.animIMPUTEDSNP)
    cat('.... Computation of Animal-specific accuracies with ONLY IMPUTED SNPs finished .......\n')
    
    Accuracy <- list(SNP_specific=Accuracy.snps,
                     Sample_specific.ALLSNP=Accuracy.animALLSNP,
                     Sample_specific.impSNP=Accuracy.animIMPUTEDSNP)
    
    cat(' \n')
    cat('.... Four types of accuracies are computed for both SNPs and animals/samples .......\n')
    cat('  1. correlation between true and imputed SNPs using calus and mulder (calus et al. 2014) ............. correl_CM \n')
    cat('  2. squared correlation between true and imputed SNPs using calus and mulder (calus et al. 2014) ..... sqrdcorrel_CM \n')
    cat('  3. Percentage (%) of correctly imputed SNPs.............. PERC \n')
    cat('  4. Allelic error rate (1-allelic error rate) ............ AER \n')
    cat(' ')
    return(Accuracy)
  }
  else if (accuracy_type=="both"){
    cat ('**********************************************************************************\n')
    cat ('*               Imputation acuracy will be computed with the                     *\n')
    cat ('*   ... 1. SIMPLE ---- (GENOTYPES ARE NOT CENTRED and SCALED)                    *\n')
    cat ('*   ... 2. calus_mulder ---- (GENOTYPES ARE CENTRED and SCALED)                  *\n')  
    cat ('*                           See the Paper                                        *\n')
    cat ('*  Eval. of measures of correctness of genotype imputation in the context of GS  *\n')
    cat ('*                 http://dx.doi.org/10.1017/S1751731114001803                    *\n')
    cat ('**********************************************************************************\n')
    
    Xtrue <- merge(idscommon,Xtrue,by=1,sort=F)[,-1]
    Ximp <- merge(idscommon,Ximp,by=1,sort=F)[,-1]
    Xtrue.calusmulder <- scale(x=Xtrue,center=T,scale=T)
    Ximp.calusmulder <- scale(x=Ximp,center=T,scale=T)
    nsnps = ncol(Xtrue)
    anim = nrow(Xtrue)
    iterchecks.snps <- round(nsnps/10,digits=-1)
    iterchecks.anim <- round(anim/10,digits=0)
    
    cat(' \n')
    #*******************   SNP imputation accuracy
    cat('....... Computation of SNP-specific accuracies started .......\n')
    cat(' \n')
    if(format[2]=='dosage'){
      Accuracy.snps <- matrix(0,nsnps,4)
      colnames(Accuracy.snps) <- c('correl_simple','sqrdcorrel_simple','correl_CM','sqrdcorrel_CM')
      for (j in 1:nsnps){
        datbind <- cbind(Ximp[,j],Xtrue[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])!=0){
          dat <- NA
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat^2
        }
        
        datbind <- cbind(Ximp.calusmulder[,j],Xtrue.calusmulder[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])==0){
          dat <- NA
          Accuracy.snps[j,3] <- dat
          Accuracy.snps[j,4] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,3] <- dat
          Accuracy.snps[j,4] <- dat^2
        }
        if(j %% iterchecks.snps==0){
          cat(paste('SNP specific accuracy now at ...',j,' ...out of a total of ...',nsnps,sep=''),' \n')
        }
      }
    }
    else {
      Accuracy.snps <- matrix(0,nsnps,6)
      colnames(Accuracy.snps) <- c('correl_simple','sqrdcorrel_simple','correl_CM','sqrdcorrel_CM','PERC',"AError")
      for (j in 1:nsnps){
        datbind <- cbind(Ximp[,j],Xtrue[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])==0){
          dat <- NA
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,1] <- dat
          Accuracy.snps[j,2] <- dat^2
        }
        
        datbind <- cbind(Ximp.calusmulder[,j],Xtrue.calusmulder[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        sddat <- as.vector(apply(datbind,2,sd)); 
        if(length(sddat[sddat!=0])==0){
          dat <- NA
          Accuracy.snps[j,3] <- dat
          Accuracy.snps[j,4] <- dat
        } else {
          dat <- cor(datbind[,1],datbind[,2])
          Accuracy.snps[j,3] <- dat
          Accuracy.snps[j,4] <- dat^2
        }
        
        datbind <- cbind(Ximp[,j],Xtrue[,j])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.snps[j,5] <- datperc
        
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.snps[j,6] <- 1-datAER
        
        if(j %% iterchecks.snps==0){
          cat(paste('SNP specific accuracy now at ...',j,' ...out of a total of ...',nsnps,sep=''),' \n')
        }
      }
    }
    
    Accuracy.snps <- round(Accuracy.snps,digits=4)
    Accuracy.snps <- cbind(mapHD[,c(1,2,4)],Accuracy.snps)
    colnames(Accuracy.snps)[1:3] <- c('CHR','SNPname','Position')
    cat('.... Computation of SNP-specific accuracies finished ....... \n')
    
    #***************************** Animal specific accuracy for all SNPS
    cat(' \n')
    cat('.... Computation of Animal-specific accuracies with ALL SNPs started .......\n')
    cat(' \n')
    
    Xtrue_T <- t(Xtrue)
    Ximp_T <- t(Ximp)
    Xtrue.calusmulder_T <- t(Xtrue.calusmulder)
    Ximp.calusmulder_T <- t(Ximp.calusmulder)
    
    if(format[2]=='dosage'){
      Accuracy.animALLSNP <- matrix(0,anim,4)
      colnames(Accuracy.animALLSNP) <- c('correl_simple','sqrdcorrel_simple','correl_CM','sqrdcorrel_CM')
      for (l in 1:anim){
        datbind <- cbind(Ximp_T[,l],Xtrue_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,1] <- dat
        Accuracy.animALLSNP[l,2] <- dat^2
        
        datbind <- cbind(Ximp.calusmulder_T[,l],Xtrue.calusmulder_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,3] <- dat
        Accuracy.animALLSNP[l,4] <- dat^2       
        if(l %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',l,'....out of a total of ...',anim,sep=''),' \n')
        }
      }      
    }
    else {
      Accuracy.animALLSNP <- matrix(0,anim,6)
      colnames(Accuracy.animALLSNP) <- c('correl_simple','sqrdcorrel_simple','correl_CM','sqrdcorrel_CM','PERC','AError')
      for (l in 1:anim){
        datbind <- cbind(Ximp_T[,l],Xtrue_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,1] <- dat
        Accuracy.animALLSNP[l,2] <- dat^2
        
        datbind <- cbind(Ximp.calusmulder_T[,l],Xtrue.calusmulder_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animALLSNP[l,3] <- dat
        Accuracy.animALLSNP[l,4] <- dat^2
        
        datbind <- cbind(Ximp_T[,l],Xtrue_T[,l])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.animALLSNP[l,5] <- datperc
        
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.animALLSNP[l,6] <- 1-datAER
        
        if(l %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',l,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    Accuracy.animALLSNP <- round(Accuracy.animALLSNP,digits=4)
    Accuracy.animALLSNP <- cbind.data.frame(idscommon,Accuracy.animALLSNP)
    assign('Accuracy.animALLSNP',Accuracy.animALLSNP)
    cat('.... Computation of Animal-specific accuracies with ALL SNPs finished .......\n')
    
    
    #*************************** Animal specific accuracy for only imputed SNPS
    map <- merge(mapHD,mapLD,by=2,sort=F)
    map <- sort(map$num)
    Xtrue_T <- Xtrue_T[-map,]
    Ximp_T <- Ximp_T[-map,]
    nsnps = nrow(Xtrue_T)
    Xtrue.calusmulder_T <- Xtrue.calusmulder_T[-map,]
    Ximp.calusmulder_T <- Ximp.calusmulder_T[-map,]
    
    #*************************  Animal specific accuracy
    cat(' \n')
    cat('.... Computation of Animal-specific accuracies with ONLY IMPUTED SNPs started .......\n')
    cat(' \n')
    
    if(format[2]=='dosage'){
      Accuracy.animIMPUTEDSNP <- matrix(0,anim,4)
      colnames(Accuracy.animIMPUTEDSNP) <- c('correl_simple','sqrdcorrel_simple','correl_CM','sqrdcorrel_CM')
      for (m in 1:anim){
        datbind <- cbind(Ximp_T[,m],Xtrue_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,1] <- dat
        Accuracy.animIMPUTEDSNP[m,2] <- dat^2
        
        datbind <- cbind(Ximp.calusmulder_T[,m],Xtrue.calusmulder_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,3] <- dat
        Accuracy.animIMPUTEDSNP[m,4] <- dat^2
        if(m %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',m,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    else {
      Accuracy.animIMPUTEDSNP <- matrix(0,anim,6)
      colnames(Accuracy.animIMPUTEDSNP) <- c('correl_simple','sqrdcorrel_simple','correl_CM','sqrdcorrel_CM','PERC','AError')
      for (m in 1:anim){
        datbind <- cbind(Ximp_T[,m],Xtrue_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,1] <- dat
        Accuracy.animIMPUTEDSNP[m,2] <- dat^2
        
        datbind <- cbind(Ximp.calusmulder_T[,m],Xtrue.calusmulder_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- cor(datbind[,1],datbind[,2])
        Accuracy.animIMPUTEDSNP[m,3] <- dat
        Accuracy.animIMPUTEDSNP[m,4] <- dat^2
        
        datbind <- cbind(Ximp_T[,m],Xtrue_T[,m])
        datbind <- datbind[which(datbind[,2]!='NA'),]
        datbind <- na.omit(datbind)
        dat <- datbind[,1]-datbind[,2]
        datperc <- length(which(dat==0))/length(dat)
        Accuracy.animIMPUTEDSNP[m,5] <- datperc
        
        datAER <- sum(abs(dat))/(2*length(dat))
        Accuracy.animIMPUTEDSNP[m,6] <- 1-datAER
        
        if(m %% iterchecks.anim==0){
          cat(paste('Animal specific accuracy now at ...',m,'....out of a total of ...',anim,sep=''),' \n')
        }
      }
    }
    Accuracy.animIMPUTEDSNP <- round(Accuracy.animIMPUTEDSNP,digits=4)
    Accuracy.animIMPUTEDSNP <- cbind.data.frame(idscommon,Accuracy.animIMPUTEDSNP)
    cat('.... Computation of Animal-specific accuracies with ONLY IMPUTED SNPs finished .......\n')
    
    Accuracy <- list(SNP_specific=Accuracy.snps,
                     Sample_specific.ALLSNP=Accuracy.animALLSNP,
                     Sample_specific.impSNP=Accuracy.animIMPUTEDSNP)
    
    cat(' \n')
    cat('.... Four types of accuracies are computed for both SNPs and animals/samples .......\n')
    cat('  1. correlation between true and imputed SNPs using calus and mulder (calus et al. 2014)............. correl_CM \n')
    cat('  2. squared correlation between true and imputed SNPs using calus and mulder (calus et al. 2014)..... sqrdcorrel_CM \n')
    cat('  3. correlation between true and imputed SNPs............. correl \n')
    cat('  4. squared correlation between true and imputed SNPs..... sqrdcorrel \n')    
    cat('  5. Percentage (%) of correctly imputed SNPs.............. PERC \n')
    cat('  6. Allelic error rate (1-allelic error rate) ............ AER \n')
    cat(' ')
    return(Accuracy)
  }
}
