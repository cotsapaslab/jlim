#' jlim class to keep the result of JLIM single test
#'
#' @slot userIdxBP user specified index SNP
#' @slot actualIdxBP actual found index SNP
#' @slot STAT jlim statistic lambda
#' @slot pvalue permutation pvalue
#' @slot startBP start position of the tested locus
#' @slot endBP end position of the tested locus
#' @slot sectrSampleSize end position of the tested locus
#' @slot sectrGeneName name of the Gene in case of multiple gene in the second trait
#' @slot sectrIndSNPpvalue pvalue of the indexSNP in the second trait.
#' @slot sectrMinpvalue minimum pvalue of in the second trait.
#' @slot sectrSNPWithMinpvalue  SNP with the minimum pvalue in the second trait
#' @slot desc  status of the JLIM test
#' @slot executedPerm number of the executed permutations
#' @slot permmat the permutation matrix used for the JLIM test
#' @export
#' @import methods

setClass("jlim",
         slots = list(
           userIdxBP="numeric",
           actualIdxBP="numeric",
           STAT="numeric",
           pvalue="numeric",
           usedSNPsNo="numeric",
           startBP="numeric",
           endBP="numeric",
           sectrSampleSize="numeric",
           sectrGeneName="character",
           sectrIndSNPpvalue="numeric",
           sectrMinpvalue="numeric",
           sectrSNPWithMinpvalue="numeric",
           desc="character",
           executedPerm="numeric",
           permmat="matrix"
         ))

#getVec.jlim <- function(object){}
setGeneric("getVec.jlim", function(object) standardGeneric("getVec.jlim"))
setMethod("getVec.jlim",
          "jlim",
          function(object) {
            c(object@userIdxBP, object@actualIdxBP, object@STAT, object@pvalue,
              object@usedSNPsNo, object@startBP, object@endBP, object@sectrSampleSize,
              object@sectrGeneName, object@sectrIndSNPpvalue, object@sectrMinpvalue,
              object@executedPerm, object@desc)
          })

#' Run Joint Likelihood Mapping (JLIM) test
#'
#' JLIM tests whether two traits - main and secondary - are driven by shared
#' causal effect or not. The method is described in Chun et al. (
#' http://biorxiv.org/content/early/2016/05/12/053165)
#'
#' @param maintr.file main trait association summary statistics file
#'
#' @param sectr.file secondary trait association summary statistics file
#'
#' @param refld.file reference LD panel file
#'
#' @param secld.file genotype file of secondary trait cohort to gather in-sample LD
#'
#' @param mainld.file genotype/ldmatix file of main trait cohort to gather in-sample LD
#'
#' @param mainld.ldmatrix type of the mainld.file parameter,TRUE=ldmatrix FALSE=genotype
#'
#' @param perm.file permutation file for secondary trait association
#'
#' @param start start position of the analysis window (bp)
#'
#' @param end end position of the analysis window (bp)
#'
#' @param r2res r2 resolution limit of JLIM analysis (default r2=0.8)
#'
#' @param perm.count Number of permutation in case of on-the-fly permutation
#'
#' @param sectr.gene.filter in case, second trait's association file contains data from multiple
#' genes, this parameter specifies that if jlim.test should be ran for all of the genes or not. (default=FALSE)
#'
#' @param geneName which gene out of all genes in the second association file should be tested
#'
#' @param indSNP base pair position of the index SNP
#'
#' @param maintr.col.names a matrix containing list of column names of the main association file to be renamed/replaced
#' CHR <-- "existing name",  BP <-- "existing name",  P <-- "existing name"
#'
#' @param sectr.col.names a matrix containing list of column names of the second association file to be renamed/replaced
#' CHR <-- "existing name",  BP <-- "existing name",  P <-- "existing name" Gene <-- "existing name", VarinatId <-- "existing name"
#'
#' @param resultFileName  output file name
#'
#' @param sectr.sample.size second traits sample size in case of on the fly permutation
#'
#' @param min.SNPs.count minimum number of the common SNPs in the main and second traits to run JLIM test minimum number of the common SNPs in the main and second traits to run JLIM test (recommended 50)..
#'
#' @param min.MAF minimum value of minor allele frequency for the second traits and reference panel(default 0.05)
#'
#' @param min.pvalue threshold for minimum second trait association pvalue to run the analysis (default 0.05)
#'
#' @param rda.file path and prefix of the file name to save JLIM objects in the rda files
#'
#' @param rda.pvalue minimum required pvalue of the JLIM result in order to save the data objects into the output files
#'
#' @recessive.model if it is  not null JLIM will consider the recessive SNPs inside the main trait and run accordingly
#'
#' @return jlim class object
#'
#' @export

jlim.test <- function(maintr.file, sectr.file, refld.file=NULL, secld.file=NULL, perm.file,
                      CHR, start.bp, end.bp, r2res=0.8, withPerm=FALSE,
                      perm.count, sectr.gene.filter=FALSE, geneName, indSNP,
                      maintr.col.names, sectr.col.names, resultFileName, sectr.sample.size,
                      min.SNPs.count, sectr.ref.db, min.MAF, min.pvalue,
                      mainld.file=NULL, mainld.ldmatrix=FALSE,rda.file=NULL,
                      rda.pvalue, recessive.model=FALSE){


  results.allgene <- matrix(ncol=13, nrow=0)
  colnames(results.allgene) <- c("userIdxBP"," actualIdxBP","STAT", "pvalue",
                                 "usedSNPsNo", "startBP","endBP", "sectrSampleSize",
                                 "sectrGeneName","sectrIdxSNPAssocPvalue", "sectrMinAssocPvalue",
                                 "executedPerm" ,"desc")

  assoc1.res <- loadMainTraitsSumStats(maintr.file,remDupBP=FALSE, CHR,
                                       col.names=maintr.col.names, indSNP)
  assoc1_b <-  assoc1.res[[1]]
  assoc1_b.org <-  assoc1.res[[2]]

  if(!is.null(mainld.file))
  {
    mainld.res <- loadMainLD( mainld.file, assoc1_b, ldmatrix=mainld.ldmatrix, start.bp, end.bp)
    ld1_b <- mainld.res[[1]]
    assoc1_b <- mainld.res[[2]]

    if(nrow(ld1_b)!=nrow(assoc1_b)){
      cat("\nNumber of SNPs in the main trait and main ld does not match.\n")
      stop(0)
    }
    assoc1_b.org <-  assoc1_b

  }

  assoc2.genes.res <- loadSecondTraitsSumStats(sectr.file, remDupBP=FALSE, CHR, start.bp, end.bp,
                                               col.names=sectr.col.names, sectr.ref.db = sectr.ref.db,
                                               min.pvalue, indSNP=indSNP)
  assoc2.genes <-  assoc2.genes.res[[1]]
  assoc2.genes.org <-  assoc2.genes.res[[2]]
  if(withPerm)
    permmat <- loadPermFile( perm.file)
  if(!is.null(secld.file)){

    secld.res <- loadSecondLD ( secld.file, assoc2.genes, assoc2.genes.org, start.bp, end.bp)
    ld2_b <- secld.res[[1]]
    ld2.gt <- secld.res[[2]]
    ld2.gtInfo <- secld.res[[3]]
    ld2.maf <- secld.res[[4]]
    assoc2.genes <- secld.res[[5]]

    if(withPerm){
      permmat <- permmat[ ,assoc2.genes.org$BP %in% assoc2.genes$BP]
    }
    assoc2.genes.org <- assoc2.genes
  }

  if(!is.null(refld.file)){
    refLD.res <- loadRefLD( refld.file,assoc1_b, start.bp, end.bp, assoc2.genes.org,
                            assoc2.genes, min.MAF, indSNP )
    refgt.org <- refLD.res[[1]]
    ld0.maf.org <- refLD.res[[2]]

    ld0.res <- processRefLD ( refld.file, assoc1_b, assoc2, refgt.org )
    ld0.org <- ld0.res[[1]]
    refgt0.org <- ld0.res[[2]]
  }else{  ## has to be a separated if after removing the possoble duplicates
    refgt.org <-ld2.gtInfo # gt + info
    refgt0.org <- ld2.gt # only gt
    ld0.org <- ld2_b
    ld0.maf.org <- ld2.maf
  }

  #  if (recessive.model){
  #    recess.res <- recessive.BP.update(assoc1_b, refgt.org, refgt0.org, ld0.org, ld0.maf.org)
  #    assoc1_b <- recess.res[[1]]
  #    refgt.org <- recess.res[[2]]
  #    refgt0.org <- recess.res[[3]]
  #    ld0.org <- recess.res[[4]]
  #    ld0.maf.org <- recess.res[[5]]
  #    recessive.BP <- recess.res[[6]]
  #  }

  #### add remove duplicates here
  if(!exists("ld1_b"))
    ld1_b <- NULL

  mainTrait.Res <- remDuplicateBPMainTrait( assoc=assoc1_b , ld=ld1_b , assocFile=maintr.file)
  assoc1_b <- mainTrait.Res[[1]]
  ld1_b <- mainTrait.Res[[2]]

  Gene.list <-  unique(assoc2.genes$Gene)
  if(sectr.gene.filter){
    if(!(geneName %in% Gene.list)){
      cat ("\nInput gene: ", geneName ," does not exist in the second trait's genes. Existing genes are below:",
           paste(Gene.list, collapse = "\n"))
      q(status=1)
    }else{
      Gene.list <- geneName
    }
  }
  resCounter <- 1
  for (gCounter in 1:length(Gene.list)){
    PL("\nGene name ",Gene.list[gCounter])
    if(gCounter==1 && !withPerm ){
      permmat <- matrix(ncol = nrow(refgt0.org), nrow=0 )
    }
    cat ("\nsectr.sample.size:", sectr.sample.size)
    results.gene <- new("jlim",
                        userIdxBP=indSNP,
                        actualIdxBP=NA_real_,
                        STAT=NA_real_, pvalue=NA_real_,
                        usedSNPsNo=NA_real_,
                        startBP= NA_real_,
                        endBP=NA_real_,
                        sectrSampleSize=sectr.sample.size,
                        sectrGeneName=Gene.list[gCounter],
                        sectrIndSNPpvalue=NA_real_,
                        sectrMinpvalue=NA_real_,
                        sectrSNPWithMinpvalue=NA_real_,
                        executedPerm=0, desc="")
    assoc2 <- assoc2.genes[assoc2.genes$Gene == Gene.list[gCounter],]
    assoc2.org <- assoc2.genes.org[assoc2.genes.org$Gene == Gene.list[gCounter],]
    if(nrow(assoc2[assoc2$BP == indSNP,])<1 ){
      exit_desc ="Second trait does not contain index SNP for the gene"
      results.gene@startBP <- min(assoc2$BP); results.gene@endBP <- max(assoc2$BP)
      cat("\n",exit_desc," :", Gene.list[gCounter] )
      results.gene@desc <- exit_desc
      results.allgene <- rbind (results.allgene, getVec.jlim(results.gene))
      next
    }

    assoc1 <- assoc1_b
    #    if (recessive.model && length(recessive.BP) >= 1){
    #      sectr.recess.SNPs <- assoc2[assoc2$BP %in%  recessive.BP,]
    #      sectr.recess.SNPs$BP <- sectr.recess.SNPs$BP +.1
    #      assoc2 <- rbind.data.frame(sectr.recess.SNPs,assoc2, stringsAsFactors=FALSE)
    #      assoc2 <- assoc2[order(assoc2$BP), ]
    #    }
    assoc2 <- remDuplicateBP(assoc2, sectr.file)
    if (nrow(assoc2)==0){
      exit_desc ="Number of remaining SNPs in the second trait for the gene is zero"
      cat("\n",exit_desc," :", Gene.list[gCounter] )
      results.gene@desc <- exit_desc
      results.allgene <- rbind (results.allgene, getVec.jlim(results.gene))
      next
    }

    if(!exists("ld2_b"))  ld2_b <- NULL
    if(!exists("mainld.ldmatrix"))  mainld.ldmatrix <- NULL

    matchedMarkers <- matchMarkers( assoc1, assoc2, start.bp, end.bp, withPerm,
                                    ld0=ld0.org, ld1=ld1_b, ld2=ld2_b, ld0.maf=ld0.maf.org,
                                    mainld.ldmatrix=mainld.ldmatrix,
                                    refgt=refgt.org, permmat)
    assoc1 <-  matchedMarkers[[1]]
    assoc2 <-  matchedMarkers[[2]]
    refgt <-  matchedMarkers[[3]]
    ld1 <-  matchedMarkers[[4]]
    ld2 <-  matchedMarkers[[5]]
    ld0.maf <-  matchedMarkers[[6]]
    if(withPerm)
      permmat <-  matchedMarkers[[7]]

    if(nrow(assoc2[assoc2$BP == indSNP,])<1 || nrow(assoc1)< min.SNPs.count){
      if (nrow(assoc2[assoc2$BP == indSNP,])<1)
        exit_desc ="Index SNP of second trait is filtered out since it is missing in refrence LD."
      else if (nrow(assoc1)< min.SNPs.count)
        exit_desc ="too few common SNPs to run JLIM"

      if (nrow(assoc1) > 0) {
        results.gene@desc <- exit_desc
        results.gene@startBP= min(assoc1$BP)
        results.gene@endBP= max(assoc1$BP)
      }
      results.allgene <- rbind (results.allgene, getVec.jlim(results.gene))
      next
    }

    jlim.res <- SNPselction ( assoc1, assoc2, ld1=ld1, ld2=ld2, ld0.maf, permmat, r2res = r2res,
                              withPerm, refgt0.org, perm.count, sectr.sample.size, min.SNPs.count,
                              min.pvalue = min.pvalue, indSNP=indSNP)
    jlim.res@userIdxBP <- indSNP
    jlim.res@sectrGeneName <- Gene.list[gCounter]
    permmat=jlim.res@permmat

    cat("\nJLIM results:",colnames(results.allgene),"\n",sep = "   ")
    cat("\n",getVec.jlim(jlim.res))
    results.allgene <- rbind (results.allgene, getVec.jlim(jlim.res))
    resCounter <- resCounter+1

    if(!is.null(rda.file) && (!is.na(jlim.res@pvalue) &&jlim.res@pvalue <= rda.pvalue))
      write.jlim.objets (rda.file, assoc1, assoc2, refgt, ld1, ld2, geneName=jlim.res@sectrGeneName)

  }
  write.table(results.allgene, file= resultFileName,
              quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

  tissueMinPvalue <- as.numeric(results.allgene[(!is.na(results.allgene[,4])),4])

  if(!is.null(rda.file) && length(tissueMinPvalue) >=1 && sum(tissueMinPvalue <= rda.pvalue)>=1 )
    save(permmat, file= paste(rda.file,"_permutaions.rda",sep=""))

  return(results.allgene)
}

write.jlim.objets <- function(rda.file, assoc1, assoc2, refgt, ld1, ld2, geneName) {

  save(assoc1, file= paste(rda.file,"_",geneName,"_mainTrait.rda",sep=""))
  save(assoc2, file= paste(rda.file,"_",geneName,"_secondTrait.rda",sep=""))
  save(refgt, file= paste(rda.file,"_",geneName,"_refgt.rda",sep=""))
  save(ld1, file= paste(rda.file,"_",geneName,"_ld1.rda",sep=""))
  save(ld2, file= paste(rda.file,"_",geneName,"_ld2.rda",sep=""))
}

calcMAF.refPanels <- function(refgt, refgt.all.sel) {

  refgt0 <- refgt[refgt.all.sel, ]
  refgt0 <- refgt0[order(refgt0$POS), ]

  ASSERT(sum(refgt.all.sel) == nrow(refgt0))

  refgt.mat <- refgt0[, 6:ncol(refgt0)]
  refgt.mat <- as.matrix(refgt.mat)

  # ref allele := 1
  refgt.mat0 <-
    t(sapply(1:nrow(refgt.mat),
             function(I) as.numeric(refgt.mat[I, ] == refgt0$REF[I]),
             simplify=TRUE))
  # alt allele := 2
  refgt.mat0[refgt.mat0 != 1] <- 2
  rownames(refgt.mat0) <- refgt0$POS

  refgt.mat <- toGT(refgt.mat0)
  rownames(refgt.mat) <- refgt0$POS

  meanF0 <- apply(refgt.mat, 1, mean)/2
  meanF0 <- pmin(meanF0, 1 - meanF0)

  meanF0
}

colNamesUpdate<- function(fileName, data, col.names){
  if(hasArg(col.names)){
    if(length(col.names)>0){
      for (l in 1:nrow(col.names)){
        if(sum(colnames(data)==col.names[l,2])==1){
          colnames(data)[colnames(data)==col.names[l,2]] <- col.names[l,1]

        }else{
          cat ("\nColumn with the name: ", col.names[l,2], "does not exist in file", fileName,"\n")
          stop()
        }
      }
    }
  }
  return(data)
}

readTraitStat <- function(fileName, CHR, col.names){

  assoc <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE)
  catE(paste("\nDim of the input file:",dim(assoc), collapse=" "))
  if(ncol(assoc)>=1){
    catE(paste("\nSample of the rows:"))
    catE(paste(colnames(assoc), collapse=" "))
    catE(paste(assoc[1,], collapse=" "))
  }
  assoc <- colNamesUpdate(fileName, data = assoc, col.names)
  if("CHR" %in% colnames(assoc) && "BP" %in% colnames(assoc))
  {
    assoc <- assoc[(assoc$CHR ==CHR | assoc$CHR ==0 ), ]

    if ("TEST" %in% colnames(assoc)) {
      MSG("\nRemove NON-ADD TESTs ", sum(assoc$TEST != "ADD"))
      assoc <- assoc[assoc$TEST == "ADD", ]
    }
    if ("Z" %in% colnames(assoc)) {
    } else if ("STAT" %in% colnames(assoc)) {
      assoc <- cbind(assoc, Z=assoc$STAT)
    } else if ("T" %in% colnames(assoc)) {
      assoc <- cbind(assoc, Z=assoc$T)
    } else if ("P" %in% colnames(assoc)) {
      assoc <- cbind(assoc, Z=PtoZ(assoc$P))
      assoc <- assoc[!is.na(assoc$P), ]
      if (sum(assoc$P == 0) > 0) {
        cat(paste("\nPvalue = 0 for", sum(assoc$P == 0),
                  "SNP(s) in ", fileName,
                  ": check for potential numerical underflow"))
      }
    } else {
      cat(paste("\nCannot find association statistics in:", fileName),"\n")
      stop()
    }
  }else{
    cat(paste("\nCannot find BP or CHR in the input file: ", fileName),"\n")
    stop()
  }
  return(assoc)
}

remDuplicateBPMainTrait <-  function(assoc, ld, assocFile){

  if (sum(is.na(assoc$P)) > 0) {
    cat("\nNumber of NA P_values:", sum(is.na(assoc$P)),"in", assocFile ,"has been removed.")
    if(!is.null(ld))
      ld <- ld[!is.na(assoc$P), !is.na(assoc$P)]
    assoc <- assoc[!is.na(assoc$P), ]
  }

  if (sum(duplicated(assoc$BP)) > 0) {
    cat(paste("\nMultiple alleles are prohibited at the same BP:", assocFile))

    BP_frqs <- data.frame(table(assoc$BP))
    duplicate_BP <- BP_frqs[BP_frqs$Freq > 1,1]
    if(!is.null(ld))
      ld <- ld[!(assoc$BP %in% duplicate_BP), !(assoc$BP %in% duplicate_BP)]
    assoc <- assoc[ !(assoc$BP %in% duplicate_BP) ,]
    cat("\n",length(duplicate_BP)," SNP's at positions ", duplicate_BP,"has been removed.")
  }
  outputlist <- list(assoc, ld)
  return(outputlist)

}

remDuplicateBP <-  function(assoc ,fileName){
  if (sum(is.na(assoc$P)) > 0) {
    cat("\nNumber of NA P_values:", sum(is.na(assoc$P)),"in", fileName ,"has been removed.")
    assoc <- assoc[!is.na(assoc$P), ]
  }

  if (sum(duplicated(assoc$BP)) > 0) {
    cat(paste("\nMultiple alleles are prohibited at the same BP:", fileName))

    BP_frqs <- data.frame(table(assoc$BP))
    duplicate_BP <- BP_frqs[BP_frqs$Freq > 1,1]
    assoc <- assoc[ !(assoc$BP %in% duplicate_BP) ,]

    cat("\n",length(duplicate_BP)," SNP's at positions ", duplicate_BP,"has been removed.")
  }
  return(assoc)
}

loadTraitsSumStatsParquet<- function(fileName, CHR, start.bp, end.bp, col.names){
  if("arrow" %in% rownames(installed.packages())){
    if(!("arrow" %in% (.packages()))){
      library(arrow)
    }
  }
  else{
    cat("\npackage arrow is required to read parquet files.")
  }

  assoc <- read_parquet(fileName, as_tibble = TRUE)
  if(ncol(assoc)==0){
    cat("\nNumber of the fetched variants from GTEx data source is Zero.\n")
    q(status=1)
  }
  if(length(col.names)==0)
    col.names <- matrix(c("P","variantId","Gene",
                          "pval_nominal","variant_id","phenotype_id"),
                        byrow = FALSE, nrow = 3,ncol = 2)

  assoc <- colNamesUpdate(fileName, data = assoc, col.names)

  if(sum(colnames(assoc)=="variantId")==1){
    BP <- sapply(assoc$variantId, function(x) substring(x,gregexpr(pattern ='_',x)[[1]][1]+1,gregexpr(pattern ='_',x)[[1]][2]-1))
    assoc <- cbind.data.frame(BP=as.numeric(BP), assoc, stringsAsFactors=FALSE)
  }
  assoc <- assoc[assoc$BP >= start.bp & assoc$BP<= end.bp, ]

  if(sum(colnames(assoc)=="variantId")==1){
    CHR <- sapply(assoc$variantId, function(x) substring(x,1,gregexpr(pattern ='_',x)[[1]][1]-1))
    ref <- sapply(assoc$variantId, function(x) substring(x,gregexpr(pattern ='_',x)[[1]][2]+1,gregexpr(pattern ='_',x)[[1]][3]-1))
    alt <- sapply(assoc$variantId, function(x) substring(x,gregexpr(pattern ='_',x)[[1]][3]+1,gregexpr(pattern ='_',x)[[1]][4]-1))
    assoc <- cbind.data.frame ( CHR, ref, alt, assoc, stringsAsFactors=FALSE)
  }
  assoc <- assoc[assoc$CHR ==CHR, ]
  assoc <- assoc [!is.na(assoc$P),]

  if ("TEST" %in% colnames(assoc)) {
    MSG("\n# of removed NON-ADD TESTs: ", sum(assoc$TEST != "ADD"))
    assoc <- assoc[assoc$TEST == "ADD", ]
  }
  if ("Z" %in% colnames(assoc)) {
  } else if ("STAT" %in% colnames(assoc)) {
    assoc <- cbind(assoc, Z=assoc$STAT)
  } else if ("T" %in% colnames(assoc)) {
    assoc <- cbind(assoc, Z=assoc$T)
  } else if ("P" %in% colnames(assoc)) {
    assoc <- cbind(assoc, Z=PtoZ(assoc$P))

    if (sum(assoc$P == 0) > 0) {
      cat(paste("\nPvalue = 0 for", sum(assoc$P == 0),
                "SNP(s) in ", fileName,
                ": check for potential numerical underflow"))
      assoc <- assoc[!is.na(assoc$P), ]
    }
  } else {
    cat(paste("\nCannot find association statistics in", fileName),"\n")
    stop()
  }

  return(assoc)
}

eQTLCatalogueExtract <- function(fileName, region,col.names){

  if("dplyr" %in% rownames(installed.packages())){
    if(!("dplyr" %in% (.packages()))){
      library(dplyr)
    }
  }
  else{
    cat("\nPackage dplyr is required to read and process eQTLCatalogue")
  }

  if("seqminer" %in% rownames(installed.packages())){
    if(!("seqminer" %in% (.packages()))){
      library(seqminer)
    }
  }
  else{
    cat("\nPackage seqminer is required to read and process eQTLCatalogue.")
  }

  cat("\nWorking dirctory: ",getwd())
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = fileName, tabixRange = region, stringsAsFactors = FALSE)
  if (exists("fetch_table") ){
    if( nrow(fetch_table)==0){
      cat("\nNumber of the fetched SNPs from eQTLCaltalogue for the input region: ",region," is Zero.\n")
      q(status=1)
    }
  }else{
    cat("\nCould not fetch any data from eQTLCaltalogue check the connection!\n")
    q(status=1)
  }

  if(all(grepl(pattern = "^http://", fileName) | grepl(pattern = "^ftp://", fileName) )){
    column_names = colnames(readr::read_tsv(fileName, n_max = 1))
  }else
    column_names = colnames(read.csv(fileName,sep = "\t",nrows=1))

  colnames(fetch_table) = column_names
  assoc = dplyr::filter(fetch_table, type%in% c("SNP"))
  catE (paste("\nColname length:",length(col.names)))

  if(length(col.names)==0)
    col.names <- matrix(c("CHR","P","BP","Gene",
                          "chromosome","pvalue","position","molecular_trait_id"),
                        byrow=FALSE, nrow=4, ncol= 2)

  assoc <- colNamesUpdate(fileName, data = assoc, col.names)

  catE(paste("\nColname parameter: ",col.names))
  catE(paste("\nColname summary stats: ",colnames(assoc)))

  assoc <- cbind(assoc, Z=PtoZ(assoc$P))
  assoc <- assoc[!is.na(assoc$P), ]

  return(assoc)
}

loadMainTraitsSumStats<- function(fileName, remDupBP=TRUE, CHR,
                                  col.names,indSNP){

  cat ("\nReading the main trait summary stat from: ", fileName )
  assoc <- readTraitStat(fileName, CHR, col.names)

  if(nrow(assoc)!=0)
    rownames(assoc) <- 1:nrow(assoc)
  else{
    cat("\nNumber of the related SNPs is zero in file: ", fileName,"\n")
    stop()
  }
  if(nrow(assoc[assoc$BP == indSNP,])<1 ){
    cat("\nIndex SNP is missing in the main trait.\n")
    stop()
  }
  assoc.org <- assoc
  if(remDupBP){
    assoc <- remDuplicateBP( assoc, fileName )
  }
  outputlist <- list(assoc, assoc.org)
  return(outputlist)
}

loadSecondTraitsSumStats<- function(fileName, remDupBP=TRUE, CHR, start.bp, end.bp,
                                    col.names, sectr.ref.db, min.pvalue, indSNP){

  if("tools" %in% rownames(installed.packages())){
    if(!("tools" %in% (.packages()))){
      library(tools)
    }
  }
  else{
    cat("\nThe required package tools does not exit.")
  }
  cat("\nReading the second trait  summary stat from: ", fileName )
  if(length(sectr.ref.db)!=0 && (sectr.ref.db=="eqtlcatalogue" || sectr.ref.db=="eqtlcatalogue:remote")){
    region <-cbind( paste(CHR,":",start.bp,"-",end.bp,sep=""))
    assoc <- eQTLCatalogueExtract (fileName, region, col.names)
  }else{
    if (file_ext(fileName)=="parquet")
      assoc <-loadTraitsSumStatsParquet(fileName, CHR, start.bp, end.bp, col.names )
    else
      assoc <- readTraitStat(fileName, CHR, col.names)
  }

  if(nrow(assoc)!=0)
    rownames(assoc) <- 1:nrow(assoc)
  else{
    cat("\nNumber of the related SNPs is zero in file: ", fileName,"\n")
    stop()
  }

  if(nrow(assoc[assoc$BP == indSNP,])<1 ){
    cat("\nIndex SNP is missing in the second trait.")
    #    stop()
  }

  if(!c("Gene") %in% names(assoc))
    assoc <- cbind.data.frame(assoc, Gene=rep("NoGeneName",nrow(assoc)),
                              stringsAsFactors=FALSE)

  assoc.org <- assoc

  if(remDupBP){
    assoc <- remDuplicateBP( assoc, fileName )
  }
  outputlist <- list(assoc, assoc.org)
  return(outputlist)
}

processRefLD <- function( refld.file, assoc1, assoc2, refgt){
  # process the refLD
  if (length(grep(".hap$", refld.file)) > 0 || length(grep(".ped$", refld.file)) > 0) {
    ld.res <- process.refLD.impute(refgt, refgt$POS %in% assoc1$BP, NULL)
    ld0 <- ld.res[[1]]
    refgt0 <- ld.res[[2]]

  } else {  # refgt is from ref Panels
    ld.res <- process.refPanels(refgt, refgt$POS %in% assoc1$BP, NULL)

    ld0 <- ld.res[[1]]
    refgt0 <- ld.res[[2]]

  }
  return(list(ld0,refgt0))
}

process.refPanels <- function(refgt, refgt.all.sel, af.check=NULL) {

  refgt0 <- refgt[refgt.all.sel, ]
  refgt0 <- refgt0[order(refgt0$POS), ]

  ASSERT(sum(refgt.all.sel) == nrow(refgt0))

  refgt.mat <- refgt0[, 6:ncol(refgt0)]
  refgt.mat <- as.matrix(refgt.mat)

  ASSERT(sum(refgt.mat == ".") == 0)
  # ref allele := 1
  refgt.mat <-
    t(sapply(1:nrow(refgt.mat),
             function(I) as.numeric(refgt.mat[I, ] == refgt0$REF[I]),
             simplify=TRUE))
  # alt allele := 2
  refgt.mat[refgt.mat != 1] <- 2
  refgt.mat0 <- toGT(refgt.mat)
  rownames(refgt.mat0) <- refgt0$POS
  rownames(refgt.mat) <- refgt0$POS
  if (!is.null(af.check)) {
    meanF0 <- apply(refgt.mat0, 1, mean)/2
    meanF0 <- pmin(meanF0, 1 - meanF0)

    af.check.tab <- cbind(refgt0[, 1:3], meanF0, af.check)
    af.check.tab <-
      af.check.tab[order(abs(meanF0 - af.check), decreasing=TRUE), ]
  }

  mode(refgt.mat0) = "numeric"
  gt0 <- stdGT2(t(refgt.mat0))
  ld0 <- (t(gt0) %*% gt0)
  var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v)))
  for (I in 1:nrow(ld0)) {
    for (J in 1:nrow(ld0)) {
      ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
    }
  }

  return(list(ld0,refgt.mat))
}
SNPselction <- function( assoc1, assoc2, ld1, ld2, ld0.maf, permmat, r2res,
                         withPerm, refgt0, perm.count, sectr.sample.size,
                         min.SNPs.count, min.pvalue, indSNP){

  assoc1.sel <- assoc1$BP[assoc1$P <= 0.1]
  assoc2.sel <- assoc2$BP[assoc2$P <= 0.1]
  #  PL("Pvalue-Threshold", min.pvalue)

  markers.t <- union(assoc1.sel, assoc2.sel)
  PL("# Markers to run JLIM on:", length(markers.t))

  ASSERT(length(markers.t) > 1)
  ASSERT(assoc1$BP[which.min(assoc1$P)] %in% markers.t)

  ld1.t <- ld1[colnames(ld1) %in% markers.t, colnames(ld1)  %in% markers.t]
  ld2.t <- ld2[colnames(ld2) %in% markers.t, colnames(ld2)  %in% markers.t]

  assoc1.t <- assoc1[assoc1$BP %in% markers.t, ]
  assoc2.t <- assoc2[assoc2$BP %in% markers.t, ]

  if (!is.null(ld0.maf)) {
    ld0.maf.t <- ld0.maf[assoc1$BP %in% markers.t]

    if (sum(is.na(ld0.maf.t)) > 0) {
      cat(paste("\nNA in ref LD matrix:", refld.file))
      stop()
    }

    if (sum(ld0.maf.t == 0) > 0) {
      cat(paste("\nMonomorphic SNPs in ref LD matrix:", refld.file, ": BP=",
                paste(names(ld0.maf.t)[ld0.maf.t == 0], collapse=", ")),"\n")
      stop()
    }
  }

  if (sum(is.na(ld1.t)) > 0) {
    cat(paste("\nNA or monomorphic SNP in ref LD matrix:", refld.file),"\n")
    stop()
  }
  best1 <- which.max(abs(assoc1.t$Z))
  sectrIndSNPpvalue <- assoc2$P[assoc2$BP==indSNP]
  sectrMinpvalue <- min(assoc2$P)
  sectrSNPWithMinpvalue <-assoc2$BP[ assoc2$P==min(assoc2$P)][1]
  jlim.res <- new("jlim",
                  userIdxBP=assoc1.t$BP[best1],
                  actualIdxBP=assoc1.t$BP[best1],
                  STAT=NA_real_, pvalue=NA_real_,
                  usedSNPsNo=nrow(assoc1.t),
                  startBP= min(assoc1.t$BP),
                  endBP= max(assoc1.t$BP),
                  sectrSampleSize=sectr.sample.size,
                  sectrGeneName="",
                  sectrIndSNPpvalue=sectrIndSNPpvalue,
                  sectrMinpvalue=sectrMinpvalue,
                  sectrSNPWithMinpvalue=sectrSNPWithMinpvalue,
                  desc="", executedPerm=0,
                  permmat=permmat)

  # check the number of remaining snps in the assoc1
  if(nrow(assoc1.t) < min.SNPs.count ){
    jlim.res@desc <- "too few SNPs to run JLIM"
    return(jlim.res)
  }

  if(min(abs(assoc2.t$P))> min.pvalue)
  {
    jlim.res@desc <- paste("Second trait min pvalue is higher than threshold:",min.pvalue,sep=' ')
    return(jlim.res)
  }

  lambda.t <- calc.stat(assoc1.t, assoc2.t, ld1.t, ld2.t, r2res)

  #########################################################################
  if(withPerm){
    NULLDIST <- perm.test(assoc1, assoc2, permmat,ld1, ld2, min.pvalue,
                          r2res, lambda.t)
    executedPerm <- nrow(permmat)
  }else{
    perm.res <- perm.test2(assoc1, assoc2, ld1, refgt0, min.pvalue,
                           r2res, lambda.t, perm.count, sectr.sample.size,
                           permmat.acc=permmat)
    NULLDIST <-  perm.res[[1]]
    permmat <-  perm.res[[2]]
    executedPerm <-  perm.res[[3]]
  }

  #  if(mean(!is.na(NULLDIST)) > 0.9){
  #    cat( "mean(!is.na(NULLDIST)) > 0.9\n")
  #  }
  permP <- sum(NULLDIST >= lambda.t, na.rm=TRUE)/sum(!is.na(NULLDIST))

  jlim.res@desc <- "executed"
  jlim.res@STAT=lambda.t
  jlim.res@pvalue=permP
  jlim.res@executedPerm=executedPerm
  jlim.res@permmat=permmat
  return(jlim.res)
}

loadMainLD<- function( file, assoc, ldmatrix, start.bp, end.bp){

  cat (paste("\nReading main ld file from :", file,"\n",sep = "" ))

  if(ldmatrix){
    catE(paste("\nMain ld file is a ld-matrix. "))
    ld  <- read.delim(file=file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    if( !(nrow(ld) == nrow(assoc) )){
      cat("\nNumber of the SNPs in main trait and its ld matrix does not match.\n")
      stop(0)
    }else
      colnames(ld) <- assoc$BP
    ld <- ld[assoc$BP >= start.bp & assoc$BP<= end.bp ,assoc$BP >= start.bp & assoc$BP<= end.bp]
  }else{
    if (length(grep(".hap$", file)) > 0) {
      legfile <-
        paste(substr(file, 1, nchar(file) - 4), ".leg", sep='')

      ASSERT(file.exists(legfile))

      ld.alinfo <- read.delim(file=legfile,
                              header=TRUE, sep=" ", stringsAsFactors=FALSE)

      ld.alinfo <- ld.alinfo[, 1:4]
      colnames(ld.alinfo)[1:4] <- c("SNP", "POS", "Allele0", "Allele1")

      ld.hap <- read.delim(file=file,
                           header=FALSE, sep=" ", stringsAsFactors=FALSE)

      ld <- cbind(ld.alinfo, ld.hap, stringsAsFactors=FALSE)

      #keep only biallelic variants
      #refgt <- refgt[!(nchar(refgt$Allele0) >1 | nchar(refgt$Allele1)>1),]

    }else if ((length(grep(".ped.gz$", file)) > 0) ||
              (length(grep(".ped$", file)) > 0)) {

      cat("\nDerive LD for main trait from plink ped...")

      ld <- load.plink.ped(file)
      ld <- as.data.frame(ld)
      ld <- cbind( CHROM="", POS=as.numeric(assoc$BP), REF="", ALT="",ld)

    } else if ((length(grep("dosage.txt$", file)) > 0) ||
               (length(grep("dosage.txt.gz$", file)) > 0) ||
               (length(grep("dosage.gz$", file)) > 0)) {

      cat("\nDerive LD for main trait from imputed genotypes...")
      dosage.res <- load.isLD1.dosage2(file, assoc$BP)
      ld <- dosage.res[[1]]
      ld.gt <- dosage.res[[2]]
      ld.gtInfo <- dosage.res[[3]]

    } else {
      PL("\nERROR: Unsupported mainld ", file)
      ASSERT(FALSE)
    }
    if (((length(grep("dosage.txt$", file)) == 0) &&
         (length(grep("dosage.txt.gz$", file)) == 0) && (length(grep("dosage.gz$", file)) == 0))){

      if (sum(duplicated(ld$POS)) > 0) {
        dup.POS <- ld$POS[duplicated(ld$POS)]

        MSG("\nRemoving multiple duplicate snps at same BP in mainld: BP =",
            paste(dup.POS, collapse=", "))
        ld <- ld[!(ld$POS %in% dup.POS), ]
      }
      assoc <- assoc[assoc$BP %in% ld$POS, ]


      mainld.res <- process.refLD.impute(ld, ld$POS %in% assoc$BP, NULL)
      ld <- mainld.res[[1]]
      ld.gt <- mainld.res[[2]]
      ld <- ld[colnames(ld) >= start.bp & colnames(ld) <= end.bp, colnames(ld) >= start.bp & colnames(ld) <= end.bp ]

    }
  }
  assoc <- assoc[assoc$BP >= start.bp & assoc$BP<= end.bp  ,]
  outputlist <- list(ld, assoc)
  return (outputlist)
}


loadSecondLD<- function( file, assoc, assoc.org, start.bp, end.bp){
  PL("\nLoading second LD file:", file)

  if (length(grep(".hap$", file)) > 0) {
    legfile <-
      paste(substr(file, 1, nchar(file) - 4), ".leg", sep='')

    ASSERT(file.exists(legfile))

    ld.alinfo <- read.delim(file=legfile,
                            header=TRUE, sep=" ", stringsAsFactors=FALSE)

    ld.alinfo <- ld.alinfo[, 1:4]
    colnames(ld.alinfo)[1:4] <- c("SNP", "POS", "Allele0", "Allele1")

    ld.hap <- read.delim(file=file,
                         header=FALSE, sep=" ", stringsAsFactors=FALSE)

    ld <- cbind(ld.alinfo, ld.hap, stringsAsFactors=FALSE)

    #keep only biallelic variants
    #refgt <- refgt[!(nchar(refgt$Allele0) >1 | nchar(refgt$Allele1)>1),]

  }else if ((length(grep(".ped.gz$", file)) > 0) ||
            (length(grep(".ped$", file)) > 0)) {

    cat("\nDerive LD for secondary trait from plink ped...")

    ld <- load.plink.ped(file)
    ld <- as.data.frame(ld)
    ld <- cbind( CHROM="", POS=as.numeric(assoc$BP), REF="", ALT="",ld)

  } else if ((length(grep("dosage.txt$", file)) > 0) ||
             (length(grep("dosage.txt.gz$", file)) > 0) ||
             (length(grep("dosage.gz$", file)) > 0)) {

    cat("\nDerive LD for second trait from imputed genotypes...")

    dosage.res <- load.isLD1.dosage2(file, assoc$BP)
    ld <- dosage.res[[1]]
    ld.gt <- dosage.res[[2]]
    ld.gtInfo <- dosage.res[[3]]

  }else{
    PL("ERROR: Unsupported secld", file)
    ASSERT(FALSE)
  }
  if (((length(grep("dosage.txt$", file)) == 0) &&  (length(grep("dosage.txt.gz$", file)) == 0) && (length(grep("dosage.gz$", file)) == 0))){
    if(nrow(ld)!=nrow(assoc)){
      cat("\nNumber of SNPs in the second traits and reference panels does not match.\n")
      stop()
    }

    if (sum(duplicated(ld$POS)) > 0) {
      dup.POS <- ld$POS[duplicated(ld$POS)]

      MSG("\nRemoving multiple duplicate snps at same BP in second traits ld: BP =",
          paste(dup.POS, collapse=", "))
      ld <- ld[!(ld$POS %in% dup.POS), ]
    }
    assoc <- assoc[assoc$BP %in% ld$POS, ]

    secld.res <- process.refLD.impute(ld, ld$POS %in% assoc$BP, NULL)
    ld <- secld.res[[1]]
    ld.gt <- secld.res[[2]]
    ld.gtInfo <- secld.res[[3]]
  }
  #calculate MAF for .hap file
  meanF0 <- apply(ld.gt, 1, mean)/2
  ld.maf <- pmin(meanF0, 1 - meanF0)

  outputlist <- list(ld, ld.gt, ld.gtInfo, ld.maf, assoc)
  return (outputlist)
}

loadPermFile<- function( perm.file)
{
  PL("Loading permutation data", perm.file)
  permmat <- load.mperm(perm.file)
  return(permmat)
}

matchMarkers<- function( assoc1, assoc2,  start.bp, end.bp, withPerm,
                         ld0, ld1, ld2, ld0.maf,
                         mainld.ldmatrix, refgt=NULL, permmat){

  catE (paste("\nBefore SNP matching #sample/SNPs in first trait:",nrow(assoc1)," second trait:",nrow(assoc2),
              " ref panel:",nrow(refgt)," start bp:", start.bp," end bp:", end.bp," offline perm:", withPerm))

  refgt_b <- refgt ; assoc1_b <- assoc1 ; assoc2_b <- assoc2
  refgt <- refgt[refgt$POS >= start.bp & refgt$POS <= end.bp, ]
  assoc1 <- assoc1[assoc1$BP >= start.bp & assoc1$BP <= end.bp, ]
  assoc2 <- assoc2[assoc2$BP >= start.bp & assoc2$BP <= end.bp, ]

  if(withPerm){
    if(ncol(permmat) == nrow(assoc2_b)){
      permmat <- permmat[, assoc2_b$BP >= start.bp & assoc2_b$BP <= end.bp]
    }else{
      cat("\nNumber of SNPs in the second trait and permutaion data is not equal.\n")
      stop(0)
    }
  }

  # common markers
  all.sel <- intersect(intersect(assoc1$BP, assoc2$BP), refgt$POS)
  all.sel <- setdiff(all.sel, assoc1$BP[is.na(assoc1$Z)])
  all.sel <- setdiff(all.sel, assoc2$BP[is.na(assoc2$Z)])

  refgt <- refgt[refgt$POS %in% all.sel, ]
  assoc1 <- assoc1[assoc1$BP %in% all.sel, ]

  if(withPerm){
    permmat <- permmat[, assoc2$BP %in% all.sel]
  }
  assoc2 <- assoc2[assoc2$BP %in% all.sel, ]

  ASSERT(nrow(refgt) == nrow(assoc1))
  ASSERT(nrow(assoc1) == nrow(assoc2))

  if(withPerm)
    ASSERT(ncol(permmat) == nrow(assoc2))

  if(!is.null(refgt)){
    ld0 <- ld0[refgt_b$POS %in% refgt$POS,refgt_b$POS %in% refgt$POS]
    ld0.maf <- ld0.maf[refgt_b$POS %in% refgt$POS]
  }
  if(!is.null(ld1)){
    ld1 <- ld1[colnames(ld1) %in% assoc1$BP, colnames(ld1) %in% assoc1$BP]
  }

  if(!is.null(ld2)){
    ld2 <- ld2[colnames(ld2) %in% assoc2$BP, colnames(ld2) %in% assoc2$BP]
  }

  # sort assoc1 by BP
  assoc1 <- assoc1[order(assoc1$BP), ]
  if(!is.null(ld1))
    ld1 <- ld1[order(colnames(ld1),decreasing=FALSE), order(colnames(ld1),decreasing=FALSE) ]

  # It should have been pre-ordered
  ASSERT(all(sort(assoc2$BP) == assoc2$BP))

  # Matching BP
  ASSERT(all(assoc1$BP == assoc2$BP))
  ASSERT(all(assoc1$BP == refgt$POS))

  if(is.null(ld1))
    ld1 <- ld0

  if(is.null(ld2))
    ld2 <- ld0

  reslist <-  list(assoc1, assoc2,refgt, ld1,ld2, ld0.maf,  permmat)

  catE(paste("After SNP Matching #sample/SNPs in first trait:",nrow(assoc1)," second trait:",nrow(assoc2),
             " ref panel:",nrow(refgt)," start bp:", start.bp," end bp:", end.bp," offline perm:", withPerm))

  return( reslist )

}

load.plink.ped <- function(file) {
  ped <- read.table(file, header=FALSE, stringsAsFactors=FALSE,
                    colClasses="character")

  gt <- ped[, 7:ncol(ped)]
  gt <- as.matrix(gt)

  # Make sure all sites are biallelic
  al1 <- gt[, seq(1, ncol(gt), by=2)]
  al2 <- gt[, seq(2, ncol(gt), by=2)]

  colnames(al2) <- colnames(al1)
  al <- rbind(al1[1,], al2[1,])

  for (i in 2:nrow(al1)){
    al <- rbind(al, al1[i,], al2[i,])
  }

  #if gt is not in 1/2 format
  if(sum(apply(al, 2, function(x) !grepl("\\D", x)))==0){
    # calculate allele count in each column of gt
    alcount <- apply(al, 2, function (X) length(setdiff(union(X, NULL), "0")))
    ASSERT(sum(alcount > 2) == 0)

    al.recode <-
      matrix(apply(al, 2, function (X) as.numeric(X != X[1])),
             byrow=FALSE, nrow=nrow(al), ncol=ncol(al))
    al.recode[al == "0"] <- NA
    gt2 <- as.matrix(t(al.recode))
    storage.mode(gt2) <- "numeric"
  }else{
    #if gt is in 1/2 format
    gt2 <- as.matrix(t(al))
    storage.mode(gt2) <- "numeric"
    gt2 <- gt2-1
  }

  return(gt2)
}

load.refLD.ped <- function(file) {
  ped <- read.table(file, header=FALSE, stringsAsFactors=FALSE,
                    colClasses="character")
  gt <- ped[, 7:ncol(ped)]
  al1 <- gt[, seq(1, ncol(gt), by=2)]
  al2 <- gt[, seq(2, ncol(gt), by=2)]
  colnames(al2) <- colnames(al1)
  al <- rbind(al1[1,], al2[1,])

  for (i in 2:nrow(al1)){
    al <- rbind(al, al1[i,], al2[i,])
  }
  refgt <- as.matrix(t(al))
  storage.mode(refgt) <- "numeric"
  refgt <- refgt-1
  return(refgt)
  # Make sure all sites are biallelic

  # select SNPs
  al1 <- al1[, snpset.flag]
  al2 <- al2[, snpset.flag]

  alcount <- apply(al, 2, function (X) length(setdiff(union(X, NULL), "0")))
  ASSERT(sum(alcount > 2) == 0)

  # Recode
  # Sample 1, allele 1 is the reference allele (:= 0), the other one is
  # alt (:= 1), missing gt is NA.
  ref <- al1[1, ]

  al.recode <-
    matrix(apply(al, 2, function (X) as.numeric(X != X[1])),
           byrow=FALSE, nrow=nrow(al), ncol=ncol(al))
  al.recode[al == "0"] <- NA

  # To genotype
  N <- nrow(al.recode) / 2
  gt.recode <-
    matrix(apply(al.recode, 2, function (X) (X[1:N] + X[(N+1):(2*N)])),
           byrow=FALSE, nrow=N, ncol=ncol(al))

  ASSERT(sum(!is.na(gt.recode) & gt.recode != 0 & gt.recode != 1 &
               gt.recode != 2) == 0)

  gt <- t(gt.recode)

}

loadRefLD<- function(refld.file, assoc1, start.bp, end.bp, assoc2.genes.org,
                     assoc2.genes, min.MAF, indSNP ){
  ld0.maf <- NULL
  PL("\nLoading reference LD ", refld.file)
  if (length(grep(".hap$", refld.file)) > 0) {
    legfile <-
      paste(substr(refld.file, 1, nchar(refld.file) - 4), ".leg", sep='')

    ASSERT(file.exists(legfile))

    refgt.alinfo <- read.delim(file=legfile,
                               header=TRUE, sep=" ", stringsAsFactors=FALSE)

    refgt.alinfo <- refgt.alinfo[, 1:4]
    colnames(refgt.alinfo)[1:4] <- c("SNP", "POS", "Allele0", "Allele1")

    refgt.hap <- read.delim(file=refld.file,
                            header=FALSE, sep=" ", stringsAsFactors=FALSE)

    refgt <- cbind(refgt.alinfo, refgt.hap, stringsAsFactors=FALSE)

    #keep only biallelic variants
    #refgt <- refgt[!(nchar(refgt$Allele0) >1 | nchar(refgt$Allele1)>1),]

  }else if ((length(grep(".ped.gz$", refld.file)) > 0) ||
            (length(grep(".ped$", refld.file)) > 0)) {

    cat("\nDerive LD for secondary trait from plink ped...")

    refgt <- load.plink.ped(file=refld.file)
    refgt <- as.data.frame(refgt)

    if(nrow(refgt)!=nrow(assoc2.genes)){
      cat("\nNumber of SNPs in the second traits and reference panels does not match.\n")
      stop()
    }
    refgt <- cbind( CHROM="", POS=as.numeric(assoc2.genes$BP), REF="", ALT="",refgt)

  } else {
    #  reference gt 1000 genome
    refgt <-
      read.delim(file=refld.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

    if (refgt[1,5]=="GT")
      colnames(refgt)[1:5] <- c("CHROM", "POS", "REF", "ALT", "GT")
    else{
      colnames(refgt)[1:5] <- c("CHROM", "POS","ID", "REF", "ALT")
      refgt <- refgt[,-c(6,7)]

      # exclude SVs
      if (sum(substr(refgt$ID, 1, 3) == "esv") > 0) {
        MSG("INFO: exclude esv in refgt",
            sum(substr(refgt$ID, 1, 3) == "esv"))
        refgt <- refgt[substr(refgt$ID, 1, 3) != "esv", ]
      }
    }
    refgt <- refgt[complete.cases(refgt[ , 6:ncol(refgt)]),]
  }

  refgt <- refgt[refgt$POS >= start.bp & refgt$POS <= end.bp, ]
  refgt <- refgt[refgt$POS %in% assoc1$BP,]

  if(nrow(refgt[refgt$POS == indSNP,])<1 ){
    cat("\nIndex SNP is missing in the reference panel file.\n")
    q(status=1)
  }

  # filter rare variant from refLD
  if (length(grep(".hap$", refld.file)) == 0 & length(grep(".ped$", refld.file)) == 0) {
    ld0.maf <- calcMAF.refPanels(refgt, rep(TRUE, nrow(refgt)))
  }else {
    #calculate MAF for .hap file
    meanF0 <- apply(refgt[,5:ncol(refgt)], 1, mean)/2
    ld0.maf <- pmin(meanF0, 1 - meanF0)
  }
  if(length(grep(".hap$", refld.file)) == 0){
    refgt <- refgt[ld0.maf >= min.MAF, ]
    ld0.maf <- ld0.maf[ld0.maf >= min.MAF]
  }
  # exclude overlapping variants
  if (sum(duplicated(refgt$POS)) > 0) {
    dup.POS <- refgt$POS[duplicated(refgt$POS)]

    MSG("\nRemoving multiple common alleles at same BP in refgt: BP =",
        paste(dup.POS, collapse=", "))

    ld0.maf <- ld0.maf[!(refgt$POS %in% dup.POS)]
    refgt <- refgt[!(refgt$POS %in% dup.POS), ]
  }

  if (sum(!(assoc1$BP %in% refgt$POS)) > 0) {
    MSG("\nMarkers missing in RefLD",
        sum(!(assoc1$BP %in% refgt$POS)))
  }

  reslist <- list( refgt, ld0.maf )
  return ( reslist )
}


calc.stat <- function (assoc1, assoc2, ld0, ld2, R2thr) {

  logP1 <- (abs(assoc1$Z) ** 2) / 2
  logP2 <- (abs(assoc2$Z) ** 2) / 2

  ASSERT(sum(is.na(assoc1$Z)) == 0)
  ASSERT(sum(is.na(assoc2$Z)) == 0)

  ### PEAK SELECTION (local)
  maxI1 <- which.max(logP1)
  relP1 <- exp(logP1 - max(logP1))
  postP1 <- exp(logP1) / sum(exp(logP1))
  local <- which(ld0[maxI1, ]**2 >= R2thr)

  gap <- 0
  sumRelP1 <- sum(relP1[local])

  for (I in local) {
    gap <- gap +
      relP1[I] * (logP2[I] - max(logP2[ld2[I, ]**2 < R2thr]))
  }

  gap <- gap / sumRelP1
  gap

  #  PL("--GAP", gap)
  gap
}

# FOR SIMULATED NULL DIST
perm.test <- function (assoc1, assoc2, permmat, ld0, ld2,
                       min.pvalue, R2thr, lambda.t) {

  thresholdingZ <- PtoZ(0.1)

  logP1 <- (abs(assoc1$Z) ** 2) / 2

  # SNP selection part 1
  markers.p1 <- abs(assoc1$Z) >= thresholdingZ

  # peak selection
  maxlogP1 <- max(logP1)
  maxI1 <- which.max(logP1)
  relP1 <- exp(logP1 - maxlogP1)

  simNo <- 1
  NULLGAP <- c()
  ASSERT(ncol(permmat) == nrow(ld0))
  ASSERT(ncol(permmat) == nrow(ld2))

  for (simNo in 1:nrow(permmat)) {

    assoc2n.Z <- permmat[simNo, ]

    ASSERT(sum(is.na(assoc2n.Z)) == 0)

    # SNP selection part 2
    markers.p <-
      markers.p1 | (abs(assoc2n.Z) >= thresholdingZ)

    ASSERT(markers.p[maxI1]) # always include maxI1

    # need multiple markers
    if (sum(markers.p) > 1) {

      local <- intersect(which(ld0[maxI1, ]**2 >= R2thr),
                         which(markers.p))

      logP2n <- (abs(assoc2n.Z) ** 2) / 2

      gap <-
        sum(relP1[local] *
              sapply(local, function (I)
                (logP2n[I] -
                   max(logP2n[(ld2[I, ]**2 < R2thr) & markers.p])),
                simplify=TRUE))

      NULLGAP[simNo] <- gap / sum(relP1[local])
    }

    # Adaptive
    if (simNo == 100 && sum(is.na(NULLGAP)) == 0) {

      pv.partial <-
        sum(NULLGAP >= lambda.t) / length(NULLGAP)

      if (pv.partial > 0.1) {
        break
      }
    }
  }

  NULLGAP
}

recessive.BP.update <- function (assoc1_b, refgt.org, refgt0.org, ld0.org, ld0.maf.org){

  recessive.BP <- assoc1_b$BP[duplicated(assoc1_b$BP)]
  if(length(recessive.BP) >= 1 ){
    assoc1_b$BP[duplicated(assoc1_b$BP)] <- assoc1_b$BP[duplicated(assoc1_b$BP)]+.1

    refgt.org.recess.SNPs <- refgt.org[refgt.org$POS %in%  recessive.BP,]
    refgt.org.recess.SNPs$POS <- refgt.org.recess.SNPs$POS +.1
    refgt.new.POS <- c(refgt.org.recess.SNPs$POS,refgt.org$POS)

    refgt0.org.recess.SNPs <- refgt0.org[refgt.org$POS %in%  recessive.BP,]
    rownames(refgt0.org.recess.SNPs) <- paste(rownames(refgt0.org.recess.SNPs), ".1", sep = "")
    refgt0.org <- rbind.data.frame(refgt0.org.recess.SNPs,refgt0.org, stringsAsFactors=FALSE)
    refgt0.org <- refgt0.org[order(refgt.new.POS), ]

    ld0.maf.org <-c(ld0.maf.org,ld0.maf.org[refgt.org$POS %in%  recessive.BP])
    ld0.maf.org <- ld0.maf.org[order(refgt.new.POS)]
    ld0.org <-cbind.data.frame(ld0.org, ld0.org[,refgt.org$POS %in%  recessive.BP], stringsAsFactors=FALSE)
    ld0.org <-rbind.data.frame(ld0.org, ld0.org[refgt.org$POS %in%  recessive.BP,], stringsAsFactors=FALSE)

    ld0.org <- ld0.org[order(refgt.new.POS), ]
    ld0.org <- ld0.org[,order(refgt.new.POS) ]

    refgt.org <- rbind.data.frame(refgt.org.recess.SNPs,refgt.org, stringsAsFactors=FALSE)
    refgt.org <- refgt.org[order(refgt.org$POS), ]
  }
  reslist <-  list(assoc1_b, refgt.org, refgt0.org, ld0.org, ld0.maf.org, recessive.BP)
  return(reslist)
}

standardize.hap <- function (refgt){
  # pre-standardize haplotypes to the mean of 0 and variance of 1/2
  # refgt: unstandardized, 2n x m (haplotypes)
  mode(refgt) = "numeric"
  refgt <- refgt-1
  allele.freq <- apply(refgt, 1, mean)
  for (I in 1:nrow(refgt)) {
    refgt[I, ] <- ((refgt[I, ] - allele.freq[I]) /
                     sqrt(2*allele.freq[I]*(1-allele.freq[I])))
  }
  return( refgt )
}

perm.gen.gt <- function (refgt, nperm, sectr.sample.size)
{
  # random y
  y <- matrix(rnorm(sectr.sample.size * nperm), ncol=sectr.sample.size, nrow=nperm)
  permmat <- matrix(0, nrow=nperm, ncol=nrow(refgt))
  refgt <- t(refgt)
  for (IP in 1:nperm) {
    # sample sectr.sample.size individuals out of refgt
    sampledgt <- refgt[sample(1:nrow(refgt), sectr.sample.size, replace=TRUE),]
    permmat[IP, ] <- y[IP, , drop=FALSE] %*% sampledgt
  }
  permmat <- permmat / sqrt(sectr.sample.size)
  return(permmat)
}

perm.test2 <- function (assoc1, assoc2, ld0, refgt,
                        min.pvalue, R2thr, lambda.t, perm.count, sectr.sample.size, permmat.acc ) {
  NULLGAP <- c()
  # In bundles of 10K
  nperm <- 10000
  maxIT <- ceiling(perm.count / nperm)
  snpIds <- rownames(refgt)
  refgt <- toGT(refgt)
  refgt <- stdGT2(t(refgt))
  refgt <- t(refgt)
  rownames(refgt) <- snpIds

  # iterations in bundles of 10K (nperm)
  for (IT in 1:maxIT) {

    if(nrow(permmat.acc)< IT*nperm ){
      permmat <- perm.gen.gt (refgt, nperm, sectr.sample.size)

      permmat.acc <- rbind(permmat.acc, permmat)
    }else{
      cat("\nRange: ",((IT-1)*nperm+1)," ",(IT*nperm),"dim(permmat.acc): ",dim(permmat.acc))
      permmat <-  permmat.acc[c(((IT-1)*nperm+1):(IT*nperm)),]
    }
    permmat <- permmat[,rownames(refgt) %in% assoc2$BP]
    thresholdingZ <- PtoZ(0.1)

    logP1 <- (abs(assoc1$Z) ** 2) / 2

    markers.p1 <- abs(assoc1$Z) >= thresholdingZ

    ## peak selection
    maxlogP1 <- max(logP1)
    maxI1 <- which.max(logP1)
    relP1 <- exp(logP1 - maxlogP1)
    ld0 <- data.matrix(ld0)
    for (simNo in 1:nrow(permmat)) {
      assoc2n.Z <- permmat[simNo, ]
      #    cat ("\r Permutation counter:",(IT-1) * 10000, " + ",simNo)
      ASSERT(sum(is.na(assoc2n.Z)) == 0)

      ## SNP selection part 2
      markers.p <-
        markers.p1 | (abs(assoc2n.Z) >= thresholdingZ)

      ASSERT(markers.p[maxI1]) # always include maxI1

      ## need multiple markers
      if (sum(markers.p) > 1) {

        local <- intersect(which(ld0[maxI1, ]**2 >= R2thr),
                           which(markers.p))

        logP2n <- (abs(assoc2n.Z) ** 2) / 2

        gap <- sum(relP1[local] *
                     sapply(local, function (I)
                       (logP2n[I] -
                          max(logP2n[(ld0[I, ]**2 < R2thr) & markers.p])),
                       simplify=TRUE))

        NULLGAP[(IT - 1) * nperm + simNo] <- gap / sum(relP1[local])
      }
    }

    ## Adaptive
    pv.partial <- sum(NULLGAP >= lambda.t) / sum(!is.na(NULLGAP))
    # after 10K perms
    if (IT == 1 && pv.partial > 0.1 ) {
      reslist <-  list(NULLGAP, permmat.acc, IT*nperm)
      return( reslist )
    }
    # after 100k perms
    if (IT == 10 && pv.partial > 0.01) {
      reslist <-  list(NULLGAP, permmat.acc, IT*nperm)
      return( reslist )
    }

  }
  reslist <-  list(NULLGAP, permmat.acc, IT*nperm)
  return( reslist )
}
