#' Joint Likelihood Mapping (JLIM) test on a single locus
#'
#' command-line interface to to run JLIM test on a single locus.
#' JLIM tests whether two traits - main and secondary - are driven by shared
#' causal effect or not. The method is described in Chun et al. (
#' http://biorxiv.org/content/early/2016/05/12/053165)
#'
#'Command line arguments:
#' \itemize{
#' \item \code{--maintr-file}: main trait association summary statistics file names including path.
#' \item \code{--sectr-file}: second trait association summary statistics file.
#' \item \code{--maintr-ld-file}: main trait ld file contains genotype or pre-calculated ld-matrix of the first cohort.
#'  It has to contain identical variants as the main trait association file. currently genotype file could be of .ped or .hap file.
#' \item \code{--maintr-ld-format}: in case maintr-ld-file  is filled with a ldmatrix file name, this argument should be assigned with ldmatrix.
#' \item \code{--sectr-ld-file}: second trait ld file contains genotype  of the second cohort.
#'  It has to contain identical variants as the second trait association file.
#' directory containing trait 2 association data
#' \item \code{--ref-ld}: a directory address that contains the complete reference panels (already constructed for the non-Finish European from  1000 genome). It could also be an already chosen reference panel that covers the locus.
#' \item \code{--index-snp}: this option accept the value in chr:BP format for index SNP.
#' Based on the specified BP, JLIM calculate an interval for its test.
#' \item \code{--manual-window-boundary}: if provided JLIM will override and readjust it's already calculated interval based on the index SNP.
#' \item \code{--perm-file}: assigning this parameter JLIM assumed offline option unless explicitly specifying the other format.
#' \item \code{--max-perm}: in case of on-the-fly permutation, JLIM expect the maximum number of the permutation
#' \item \code{--sectr-ref-db}: currently JLIM apart from the normal association file specified by the user
#' it can accept the second trait association file from two different file source (eQTL catalog or GTEx data).
#' in case user likes to test his primary trait against eQTL data catalog, he has to assign eQTL to this parameter.
#' for GTEx data source JLIM can handle it based the file extension.
#' eQTL has .tsv indexed file format which need .tbi could be read directly read drom it's original file source,
#' GTEx has .parquat file format which can only be read locally]
#' JLIM automatically extract and assign the column names of these two formats but to cover different and new data sources, user has to specifies column names in the below arguments.
#' \item \code{--sectr-colname-variant-id}: in case of reading GTEx .parquat format user has to specify the name of the variant-id column.
#' \item \code{--sectr-sample-size}: in on-the-fly permutation mode, to increase the accuracy of the results user has to specify the sample size of second cohorts that second association file is produced from.
#' \item \code{--sectr-gene-filter}: second trait's gene name to filter in case it contains multiple genes.
#' \item \code{--output-file}: name of the output file
#' \item \code{--maintr-colname-chr}: main trait's column name which contains chromosome Id
#' \item \code{--maintr-colname-bp}: main trait's column name which contains base pair position of variant.
#' \item \code{--maintr-colname-p}: main trait's column name which contains pvalue.
#' \item \code{--sectr-colname-chr}: second trait's column name which contains chromosome Id
#' \item \code{--sectr-colname-gene}: second trait's column name which contains gene Id
#' \item \code{--sectr-colname-bp}: second trait's column name which contains base pair position of variant.
#' \item \code{--sectr-colname-p}: second trait's column name which contains pvalue.
#' \item \code{--r2-resolution}: r2 resolution limit of JLIM analysis (default r2=0.8)
#' \item \code{--window-size}: the locus size to run JLIM test (default 200k)
#' \item \code{--min-MAF}: minimum value of minor allele frequency for the second traits and reference panel(default  0.05)
#' \item \code{--min-SNPs-count}: minimum number of the common SNPs in the main and second traits to run JLIM test (default 50).
#' \item \code{--min-pvalue}: pvalue minimum threshold for second trait association pvalue to run the
#' analysis (default 0.05)
#' \item \code{--save-to-rda-file}: path and prefix of the file name to save JLIM objects in the rda files
#' \item \code{--save-to-rda-p}: minimum required pvalue of the JLIM result in order to save the data objects into the output files
#' \item \code{--recessive-model}: if it is not null JLIM will consider the recessive SNPs inside the main trait and run accordingly
#' }
#'
#' \code{--maintr-file --sectr-file  --index-snp} are required
#' arguments, and the rest of the arguments could be optional or mandatory based on the requstedpermutation mode.
#' @export
#' @import getopt
#' @import arrow
#' @import tools
#' @import seqminer
#' @import dplyr
run.jlim <-  function() {
  library('getopt')
  spec = matrix(c(
    'verbose' , 'v', 0, "character",
    'help' , 'h', 0, "character",
    'maintr-file', 'm', 1, "character",
    'sectr-file' , 's', 1, "character",
    'index-snp' , 'i', 1, "character",
    'manual-window-boundary','c', 1, "character",
    'ref-ld' , 'r', 1, "character",
    'sectr-ld-file' , 't', 1, "character",
    'maintr-ld-file', 'q', 1, "character",
    'maintr-ld-format', 'u', 0, "character",
    'perm-file' , 'p', 1, "character",
    'max-perm' , 'k', 1, "character",
    'r2-resolution' , 'z', 1, "character",
    'window-size' , 'w', 1, "numeric",
    'sectr-sample-size', 'n', 1, "numeric",
    'sectr-gene-filter','g',1,"character",
    'output-file', 'j', 1, "character",
    'maintr-colname-chr', '1', 1, "character",
    'maintr-colname-bp', '2', 1, "character",
    'maintr-colname-p', '3', 1, "character",
    'sectr-colname-chr', '4', 1, "character",
    'sectr-colname-gene', '5', 1, "character",
    'sectr-colname-bp', '6', 1, "character",
    'sectr-colname-p', '7', 1, "character",
    'sectr-colname-variant-id', '8', 1, "character",
    'sectr-ref-db','d',1,"character",
    'min-MAF', 'e', 1, "numeric",
    'min-pvalue', 'f', 1, "numeric",
    'min-SNPs-count', 'l', 1, "numeric",
    'save-to-rda-file','a',1,"character",
    'save-to-rda-p','b',1,"numeric"
#    'recessive-model', '9', 0, "character"
  ), byrow=TRUE, ncol=4)

  default.min.SNPs.count <- 50
  default.window.size <-200000
  default.max.perm <-  10000
  default.r2res <- 0.8
  default.min.MAF <- 0.05
  default.min.pvalue <- 0.05
  default.rda.pvalue <- 0.0001
  default.recessive.model <- F
  args <- commandArgs(TRUE)

  if (length(args) <= 1) {
    cat("\n",getopt(spec, usage=TRUE, command="run.jlim() "))
    q(status=1)
  }
  args <- args[2:length(args)]
  opt <- getopt(spec, opt=args)

  if ( !is.null(opt$help) ) {
    cat ("\nLook at the following link for the details: https://github.com/cotsapaslab/jlim/tree/2.0.1")
    q(status=1)
  }

  if ( !is.null(opt$verbose) )
    verbose <<- T
  else
    verbose <<- F

  if (length(args) < 3) {
    cat("\nRequired options are missing: --maintr-file --sectr-file  --index-snp")
    q(status=1)
  }

  if (is.null(opt[["maintr-file"]]) || is.null(opt[["sectr-file"]])
      || is.null(opt[["index-snp"]]) ) {
    cat("\n",getopt(spec, usage=TRUE, command="jlim.Run.sh"))
    cat("\nRequired options: --maintr-file --sectr-file  --index-snp")
    q(status=1)
  }

  maintr.file <- opt[["maintr-file"]]
  sectr.file <- opt[["sectr-file"]]
  refld.file <- NULL
  ref.LD <- NULL
#######################
  if ( !is.null(opt[["index-snp"]]) ) {
    leadSNP = strsplit( opt[["index-snp"]], ':', fixed=TRUE)[[1]]
    if(!is.na(leadSNP[1])){
      CHR <- as.numeric(gsub("chr", "", tolower(leadSNP[1])))
    }else{
      cat("\nRequired options: --index-snp has to be given in the format: CHR#:bp")
      q(status=1)
    }
    if(!is.na(leadSNP[2])){
      if(!grepl("\\D", leadSNP[2])){
        indexSNP <- as.numeric( leadSNP[2])
      }else{
        cat("\nIndex SNP BP is not a valid number")
        q(status=1)
      }
    }else{
      cat("\nIndex SNP BP does not exist.")
      q(status=1)
    }
  }

  if(!is.null(opt[["window-size"]])){
    window.size <- opt[["window-size"]]
  }else{
    window.size<- default.window.size
  }

  if ( !is.null(opt[["manual-window-boundary"]]) ) {
    window.boundary = strsplit( opt[["manual-window-boundary"]], '-', fixed=TRUE)[[1]]
    if(!grepl("\\D",window.boundary[1])){
      start.bp <- as.numeric(window.boundary[1])
    }else{
      cat("\n--manual-window-boundary has to be given in the format: startBP-endBP")
      q(status=1)
    }
    if(!grepl("\\D",window.boundary[2])){
      end.bp <- as.numeric(window.boundary[2])
    }else{
      cat("\n--manual-window-boundary has to be given in the format: startBP-endBP")
      q(status=1)
    }
  }else{
    start.bp = max(1, indexSNP - (window.size %/% 2))
    end.bp = indexSNP + (window.size %/% 2)
  }

  if (!is.numeric(start.bp) || !is.numeric(end.bp)) {
    cat(paste("\nInterval boundaries are not valid numbers", start.bp," - ",end.bp))
    q(status=1)
  }else if(start.bp > end.bp){
    cat(paste("\nInvalid interval", start.bp," - ",end.bp))
    q(status=1)
  }
#######################
  if ( !is.null(opt[["sectr-sample-size"]]) ) {
    sectr.sample.size <- opt[["sectr-sample-size"]]
  }else
    sectr.sample.size<-NULL

  if ( is.null(opt[["ref-ld"]]) ) {
    cat("\nMissing options: --ref-ld ")
  } else{
    ref.LD <- opt[["ref-ld"]]
    if (is.na(file.info(ref.LD)$isdir)) {
      cat(paste("\nFile or directory does not exist:", ref.LD))
      q(status=1)
    }else if (!file.info(ref.LD)$isdir) {
      refld.file <- ref.LD
      catE(paste("\nInput ref file:", refld.file ))
    }else{
      refld.file <- find.panel(ref.LD, start.bp, end.bp, CHR)
      catE(paste("\nSelected ref panel:", refld.file ))
    }
  }

  if( !is.null(opt[["sectr-ref-db"]]) ){
    sectr.ref.db <-opt[["sectr-ref-db"]]
    cat ("\nSecond trait is extracted from data source: ", sectr.ref.db)
    sectr.ref.db <-tolower( sectr.ref.db )

    if(sectr.ref.db %in%c("gtex.v8.eur","eqtlcatalogue:remote","eqtlcatalogue")){
      if(sectr.ref.db!="eqtlcatalogue"){
        sectr.res <- sectr.sample.size.lookup(sectr.file ,sectr.ref.db, sectr.sample.size)
        sectr.file <-  sectr.res[[2]]
        if(is.null(sectr.sample.size))
          sectr.sample.size <-  sectr.res[[1]]
      }
    }else{
      cat("\n--sectr-ref-db options should be one of the following data sources (case-insensitive): GTEx.V8.EUR, eQTLCatalogue or eQTLCatalogue:remote.")
      q(status=1)
    }
  }else
    sectr.ref.db<-c()

#### if on the fly permutation is requested
if (is.null(opt[["perm-file"]])){
    secld.file <-NULL
    withPermFile <- FALSE
    if ( !is.null(opt[["max-perm"]]) ){

      max.perm <- get.MaxPerm(opt[["max-perm"]])
    }else{
      max.perm <- default.max.perm #  assign the default number of perm
    }
  } else{
    withPermFile <-TRUE
    perm.file <- opt[["perm-file"]]
    sectr.sample.size <- NA_real_
  }
  ##########
  if(!is.null(opt[["sectr-ld-file"]]))
    secld.file <- opt[["sectr-ld-file"]]

  #############

  if( !is.null(maintr.file))  fileExists(maintr.file)
  if(!is.null(sectr.file) && !all(grepl(pattern = "^http://", sectr.file) || grepl(pattern = "^ftp://", sectr.file) ))
    fileExists(sectr.file)

  if( withPermFile) fileExists(perm.file)
  if( !is.null(secld.file)) fileExists(secld.file)

### check the parameter with the default option
  if(!is.null(opt[["r2-resolution"]])){
    r2res <- opt[["r2-resolution"]]
  }else{
    r2res <- default.r2res  # set default value for r2-resolution
  }

  if(!is.null(opt[["min-pvalue"]])){
    min.pvalue <- as.numeric(opt[["min-pvalue"]])
  }else{
    min.pvalue<- default.min.pvalue
  }
  cat("\nPvalue threshold: ",min.pvalue)

  if(!is.null(opt[["min-MAF"]])){
    min.MAF <- as.numeric(opt[["min-MAF"]])
  }else{
    min.MAF<- default.min.MAF
  }

  if(!is.null(opt[["min-SNPs-count"]])){
    min.SNPs.count <- as.numeric(opt[["min-SNPs-count"]])
  }else{
    min.SNPs.count<- default.min.SNPs.count
  }

  PL("##### Input Parameters #####")
  PL("primary trait sumarry stat",maintr.file )
  PL("secondary trait summary stat", sectr.file)
  PL("reference LD", ref.LD)
  PL("CHR ", CHR)
  PL("index SNP BP", indexSNP)
  PL("second trait sample size", sectr.sample.size)

  if(!withPermFile){
    PL("permutation type", "on the fly perm")
  }else {
    PL("permutation type", "offline perm")
    PL("permutation file", perm.file)
    PL("second ld file", secld.file)
  }

  PL("start (BP)", start.bp)
  PL("end (BP)", end.bp)
  PL("r2 resolution", r2res)
  PL("window-size", window.size)
  if ( !is.null(opt[["max-perm"]]) )
    PL("max.perm", max.perm)

  if ( !is.null(opt[["sectr-colname-gene"]]) ){
    sectr.gene.col.name <- opt[["sectr-colname-gene"]]
  }

  if ( !is.null(opt[["sectr-gene-filter"]]) ) {
    geneNameToFilter <- opt[["sectr-gene-filter"]]
    sectr.gene.filter <- TRUE
  }else{
    sectr.gene.filter <- FALSE
  }
  PL ("second trait Gene filter(False=All Genes) ", sectr.gene.filter)

  if ( !is.null(opt[["output-file"]]) ){
    resultFileName <- opt[["output-file"]]
    if(!file.exists(dirname(resultFileName) ) )
    {
      cat("\nOutput file's path is not valid: ", dirname(resultFileName))
      q(status=1)
    }
  }else{
    resultFileName <- paste(getwd(),"/jlimResult.txt",sep="")
    cat("\n--output-file missing: Writing to ", resultFileName)
  }
  PL ("output file ", resultFileName)

  maintr.col.names <- vector()
  if ( !is.null(opt[["maintr-colname-chr"]]) ){
    maintr.col.names <- rbind(maintr.col.names, c("CHR",opt[["maintr-colname-chr"]]))
  }
  if ( !is.null(opt[["maintr-colname-bp"]]) ){
    maintr.col.names <- rbind(maintr.col.names, c("BP",opt[["maintr-colname-bp"]]))
  }
  if ( !is.null(opt[["maintr-colname-p"]]) ){
    maintr.col.names <- rbind(maintr.col.names, c("P",opt[["maintr-colname-p"]]))
  }

  sectr.col.names <- vector()
  if ( !is.null(opt[["sectr-colname-chr"]]) ){
    sectr.col.names <- rbind(sectr.col.names,c("CHR",opt[["sectr-colname-chr"]]))
  }
  if ( !is.null(opt[["sectr-colname-gene"]]) ){
    sectr.col.names <-rbind(sectr.col.names,c("Gene",opt[["sectr-colname-gene"]]))
  }
  if ( !is.null(opt[["sectr-colname-bp"]]) ){
    sectr.col.names <-rbind(sectr.col.names,c("BP",opt[["sectr-colname-bp"]]))
  }
  if ( !is.null(opt[["sectr-colname-p"]]) ){
    sectr.col.names <-rbind(sectr.col.names,c("P",opt[["sectr-colname-p"]]))
  }
  if ( !is.null(opt[["sectr-colname-variant-id"]]) ){
    sectr.col.names <-rbind(sectr.col.names,c("variantId",opt[["sectr-colname-variant-id"]]))
  }

  if ( !is.null(opt[["parse-sectr-variant-id"]]) ){
    parse.sectr.variant.id <-opt[["parse-sectr-variant-id"]]
  }


  if ( !is.null(opt[["maintr-ld-file"]]) ){
    mainld.file <- opt[["maintr-ld-file"]]

    if ( !is.null(opt[["maintr-ld-format"]]) &&
         (tolower(opt[["maintr-ld-format"]])=="ldmatrix"))
      mainld.ldmatrix <- TRUE
    else
      mainld.ldmatrix <- FALSE
  }else{
    mainld.file <- NULL
    mainld.ldmatrix <- FALSE
  }

  PL ("Main trait LD file ", mainld.file)
#  PL ("Is main trait LD file a genotype file ", !mainld.ldmatrix)
  if(!is.null(opt[["save-to-rda-p"]])){
    rda.pvalue <- as.numeric(opt[["save-to-rda-p"]])
  }else{
    rda.pvalue <- default.rda.pvalue
  }

  if ( !is.null(opt[["save-to-rda-file"]]) ){
    rda.file <- opt[["save-to-rda-file"]]
  }else{
    rda.file <- NULL
  }

  if(!is.null(rda.file)){
    PL (".Rda files path and prefix ", rda.file)
    PL ("Min pvalue to write rda files for each run ", rda.pvalue)
  }
#  if (!is.null(opt[["recessive-model"]])){
#    if ( tolower(opt[["recessive-model"]])=="t" || tolower(opt[["recessive-model"]])=="true")
#      recessive.model<- T
#    else
#      recessive.model<- F
#  }else
#    recessive.model <- default.recessive.model
#  PL ("recessive model: ", recessive.model)

  results.allgenes <-
    jlim.test(maintr.file=maintr.file, sectr.file=sectr.file, refld.file=refld.file,
              secld.file=secld.file, perm.file=perm.file,
              CHR, start.bp=start.bp, end.bp=end.bp,
              r2res=r2res, withPermFile, perm.count=max.perm,
              sectr.gene.filter, geneName=geneNameToFilter, indSNP=indexSNP,
              maintr.col.names=maintr.col.names, sectr.col.names=sectr.col.names,
              resultFileName, sectr.sample.size, min.SNPs.count=default.min.SNPs.count,
              sectr.ref.db=sectr.ref.db, min.MAF= min.MAF , min.pvalue= min.pvalue,
              mainld.file=mainld.file, mainld.ldmatrix= mainld.ldmatrix, rda.file=rda.file,
              rda.pvalue=rda.pvalue, recessive.model=recessive.model)

}

# fetch the sample size for eqtlcatalogue from the ftp sever and also GTEx as it is hard-coded
sectr.sample.size.lookup <-function( sectr.file ,sectr.ref.db, sectr.sample.size ) {

  if (tolower(sectr.ref.db) =="eqtlcatalogue:remote"){
    tabix_paths = read.delim(
      "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv",
      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    tabix_paths <- tabix_paths[tabix_paths$quant_method=="ge",]

    tissueNames <- sapply(tabix_paths$ftp_path, function(x) substring(x,regexpr("\\/[^\\/]*$",x) +1,gregexpr(pattern ='.all.tsv.gz',x)[[1]][1]-1))

    tabix_paths <-  cbind.data.frame(tabix_paths, name=tolower(tissueNames))
    rownames(tabix_paths) <- NULL

    if(sum(tabix_paths$name == tolower(sectr.file))==1){
      sectr.sample.size <- tabix_paths$sample_size[tabix_paths$name == tolower(sectr.file)]
      sectr.file <- tabix_paths$ftp_path[tabix_paths$name == tolower(sectr.file)]
      cat ("\nSecond trait sample size fetched from online eQTLCatalogue: ", sectr.sample.size)
      cat ("\nSecond trait ftp path fetched from online eQTLCatalogue: ", sectr.file)
    }else{
      cat ("\nThe input file name does not match any entry in the eqtlcatalogue.")
      q(status=1)
    }
  }
  if (sectr.ref.db %in%c("gtex.v8.eur")){
    gtex.sample.size <- readr::read_delim(
      '  tissue                 |   sample_size
       # ----------------------------
      Adipose_Subcutaneous      |   479
      Adipose_Visceral_Omentum  |   393
      Adrenal_Gland             |   194
      Artery_Aorta              |   329
      Artery_Coronary           |   175
      Artery_Tibial	            | 476
      Brain_Amygdala	          | 119
      Brain_Anterior_cingulate_cortex_BA24	| 135
      Brain_Caudate_basal_ganglia	| 172
      Brain_Cerebellar_Hemisphere	| 157
      Brain_Cerebellum	          | 188
      Brain_Cortex	              | 183
      Brain_Frontal_Cortex_BA9  	| 157
      Brain_Hippocampus	          | 150
      Brain_Hypothalamus	        | 156
      Brain_Nucleus_accumbens_basal_ganglia	| 181
      Brain_Putamen_basal_ganglia	          | 153
      Brain_Spinal_cord_cervical_c-1      	| 115
      Brain_Substantia_nigra	              | 100
      Breast_Mammary_Tissue               	| 329
      Cells_Cultured_fibroblasts	          | 403
      Cells_EBV-transformed_lymphocytes     | 113
      Colon_Sigmoid	            | 266
      Colon_Transverse	        | 294
      Esophagus_Gastroesophageal_Junction	  | 275
      Esophagus_Mucosa	        | 411
      Esophagus_Muscularis	    | 385
      Heart_Atrial_Appendage  	| 316
      Heart_Left_Ventricle	    | 327
      Kidney_Cortex	            | 65
      Liver	                    | 178
      Lung                    	| 436
      Minor_Salivary_Gland	    | 114
      Muscle_Skeletal	          | 588
      Nerve_Tibial	            | 438
      Ovary   	  | 138
      Pancreas	  | 243
      Pituitary	  | 219
      Prostate	  | 181
      Skin_Not_Sun_Exposed_Suprapubic	| 430
      Skin_Sun_Exposed_Lower_leg	    | 508
      Small_Intestine_Terminal_Ileum	| 141
      Spleen	    | 179
      Stomach	    | 260
      Testis	    | 272
      Thyroid	    | 482
      Uterus      | 107
      Vagina	    | 120
      Whole_Blood	| 558',
      trim_ws = TRUE, comment="#", delim="|")

    gtex.sample.size$tissue <- tolower(gtex.sample.size$tissue)
    tissueName <- tolower(substring(sectr.file,regexpr("all_associations_",sectr.file)+17,
                                    gregexpr(pattern ='.v8.EUR.allpairs',sectr.file)[[1]][1]-1))
    if(sum(gtex.sample.size$tissue == tissueName)==1){
      sectr.sample.size <- gtex.sample.size$sample_size [gtex.sample.size$tissue == tissueName]
      cat ("\nSecond trait sample size for the associated GTEx.V8.EUR is: ",sectr.sample.size)
    }else if(is.null(sectr.sample.size)){
      cat ("\nSecond trait sample size is missing.")
      cat ("\nThe input file name does not match any entry in the GTEx.V8.EUR")
      q(status=1)
    }

  }
  return(list(sectr.sample.size ,sectr.file))
}

# parse the input value for the max.perm( 1k -> 10^3 or  1M-> 10^6)
get.MaxPerm<-function(inputMaxPerm){
  if(nchar(gsub("[^a-z]", "", tolower(inputMaxPerm)))>0) # if inputMaxPerm contains character
  {
    chrPart <- sub("\\d*", "", tolower(inputMaxPerm))
    digitPart <- sub("\\D.*", "", tolower(inputMaxPerm))
    if(nchar(chrPart) == 1 &&(chrPart=="m" ||chrPart=="k")){
      chrPartTransformed <- if(chrPart=="m")"1e6"else"1e3"
      max.perm <- as.numeric(digitPart)*as.numeric(chrPartTransformed)
    }else{
      cat("\nUnsupported format for max permutation number")
      stop()
    }
  }else{
    max.perm <- as.numeric(inputMaxPerm)
  }
  return(max.perm)
}


# find the proper ref-panel from the list of the files in the specified path
find.panel <-function( refLD.dir, start.bp, end.bp, CHR ) {

  ref.LD.list <- vector()
  ref.LD.list0 <-
    list.files(refLD.dir,
               paste("[^.]+.*.txt.gz", sep=''))


  if (length(ref.LD.list0) == 0) {
    PL("No ref panel exist in the folder ", tr2subdir)
    next
  }
  ref.LD.list <- gsub("_", ".", ref.LD.list0)

  panels.List<-
    sapply(ref.LD.list,
           function(str) strsplit(str, '.', fixed=TRUE),
           simplify=TRUE)

  panels <- as.data.frame(matrix(unlist(panels.List), ncol  = 5, byrow = TRUE))
  panels[,2] <- as.numeric(as.character(panels[,2]))
  panels[,3] <- as.numeric(as.character(panels[,3]))

  panel.se <- which(panels[,1]==paste("chr",CHR,sep="") & panels[,2] <= start.bp & panels[,3]>=end.bp)

  if(length(panel.se)==0){
    cat("\nThere is not any matching panel for the given interval: ", start.bp," - ",end.bp )
    q(status=1)
  }else if(length(panel.se)==1){
    panel.toUse <- panel.se
  }else{
    panel.toUse <- panel.se[1]
  }

  return(paste(refLD.dir,ref.LD.list0[panel.toUse],sep=""))
}

