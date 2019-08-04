#' jlim class to contain the processed input data and results of JLIM test
#'
#' @slot assoc1 data.frame of trait 1 with markers acutally used
#' @slot assoc2 data.frame of trait 2 with markers acutally used
#' @slot ld0 reference LD of markers actually used
#' @slot ld2 in-sample LD of trait 2 cohort only with markers actually used
#' @slot idx2BP BP position of most associated SNP in assoc1
#' @slot STAT jlim statistic lambda
#' @slot p.value permutation p-value
#' 
#' @export
#' @import methods
setClass("jlim",
         slots = list(
             # internal data
             assoc1="data.frame",
             assoc2="data.frame",
             ld0="matrix",
             ld2="matrix",

             # results
             idx2BP="numeric",
             STAT="numeric",
             p.value="numeric"
             ))

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
#' @param refld.file file of reference LD panel
#'
#' @param secld.file genotype file of secondary trait cohort to gather
#' in-sample LD
#'
#' @param perm.file permutation file for secondary trait association
#'
#' @param start position at the start of analysis window (bp)
#' 
#' @param end position at the end of analysis window (bp)
#'
#' @param r2res r2 resolution limit of JLIM analysis (default r2=0.8)
#'
#' @param mainld.file genotype file of main trait cohort to gather
#' in-sample LD if provided. This is not needed in general.
#'
#' @return jlim class object
#' @export
jlim.test <- function(maintr.file, sectr.file, refld.file,
                      secld.file, perm.file,
                      start, end, r2res=0.8, mainld.file=NULL) {

    ###########################################################################
    # START
    # Load data files (assoc1, assoc2, refld, and permutation)

    # Load trait 1
    assoc1 <- read.table(maintr.file, header=TRUE, stringsAsFactors=FALSE)

    # Note: PtoZ doesn't account for direction of effect
    if ("Z" %in% colnames(assoc1)) {
        # use provided Z 
    } else if ("STAT" %in% colnames(assoc1)) {
        assoc1 <- cbind(assoc1, Z=assoc1$STAT)
    } else if ("T" %in% colnames(assoc1)) {
        assoc1 <- cbind(assoc1, Z=assoc1$T)
    } else if ("P" %in% colnames(assoc1)) {
        assoc1 <- cbind(assoc1, Z=PtoZ(assoc1$P))
        
        if (sum(assoc1$P == 0) > 0) {
            cat(paste("P-value = 0 for", sum(assoc1$P == 0),
                      "SNP(s) in ", maintr.file,
                      ": check for potential numerical underflow\n"))
            stop()
        }
    } else {
        cat(paste("Cannot find association statistics in", maintr.file))
        stop()
    }

    # Save table containing the full set of markers in their original order
    rownames(assoc1) <- 1:nrow(assoc1)
    assoc1.org <- assoc1


    ###
    # Load ref LD
    PL("Loading reference LD:", refld.file)
    
    if (length(grep(".hap$", refld.file)) > 0) {

        # .leg
        legfile <-
            paste(substr(refld.file, 1, nchar(refld.file) - 4), ".leg", sep='')

        ASSERT(file.exists(legfile))
    
        refgt.alinfo <- 
            read.delim(file=legfile,
                       header=TRUE, sep=" ", stringsAsFactors=FALSE)
    
    #    refgt.altype <- refgt.alinfo[, 5]
        refgt.alinfo <- refgt.alinfo[, 1:4]
        colnames(refgt.alinfo)[1:4] <- c("SNP", "POS", "Allele0", "Allele1")
        #cat(dim(refgt.alinfo))
        #cat("\n")

        # .hap
        refgt.hap <- 
            read.delim(file=refld.file,
                       header=FALSE, sep=" ", stringsAsFactors=FALSE)
        dim(refgt.hap)

        refgt <-
            cbind(refgt.alinfo, refgt.hap, stringsAsFactors=FALSE)

    } else {
        # 1kg
    
        refgt <-
            read.delim(file=refld.file, header=FALSE, sep="\t",
                       stringsAsFactors=FALSE)
        colnames(refgt)[1:7] <-
            c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")

        # exclude SVs (waste memory space)
        if (sum(substr(refgt$ID, 1, 3) == "esv") > 0) {
            MSG("INFO: exclude esv in refgt",
                sum(substr(refgt$ID, 1, 3) == "esv"))
            
            refgt <- refgt[substr(refgt$ID, 1, 3) != "esv", ]
        }
    }

    if (sum(!(assoc1$BP %in% refgt$POS)) > 0) {
        MSG("Markers missing in RefLD",
              sum(!(assoc1$BP %in% refgt$POS)))
    }


    ###
    # Load trait 2
    assoc2 <- read.table(sectr.file, header=TRUE, stringsAsFactors=FALSE)

    # covariates from plink
    if ("TEST" %in% colnames(assoc2)) {
        MSG("Remove NON-ADD TESTs: ", sum(assoc2$TEST != "ADD"))
        assoc2 <- assoc2[assoc2$TEST == "ADD", ]
    }

    if ("Z" %in% colnames(assoc2)) {
        # Use provided Z    
    } else if ("STAT" %in% colnames(assoc2)) {
        assoc2 <- cbind(assoc2, Z=assoc2$STAT)
    } else if ("T" %in% colnames(assoc2)) {
        assoc2 <- cbind(assoc2, Z=assoc2$T)
    } else {
        PL("CAN'T FIND STATISTIC", sectr.file)
        ASSERT(FALSE)
    }

    # Save table containing the full set of markers in their original order
    rownames(assoc2) <- 1:nrow(assoc2)
    assoc2.org <- assoc2


    ###
    # Load permutation

    PL("Loading permutation data", perm.file)

    permmat <- load.mperm(perm.file)

    ASSERT(ncol(permmat) == nrow(assoc2))

    #PL("assoc1", paste(dim(assoc1), collapse=", "))
    #PL("Ref-LD dim", paste(dim(refgt), collapse=", "))    
    #PL("assoc2", paste(dim(assoc2), collapse=", "))

    ##########################################################################
    # Match and align markers across tables

    if (sum(is.na(assoc1$P)) > 0) {
        MSG("NA P1:", sum(is.na(assoc1$P)))
        assoc1 <- assoc1[!is.na(assoc1$P), ]
    }

    if (sum(is.na(assoc2$P)) > 0) {
        MSG("NA P2", sum(is.na(assoc2$P)))
        assoc2 <- assoc2[!is.na(assoc2$P), ]
    }

    # filter rare from refgt
    ld0.maf <- NULL
    
    if (length(grep(".hap$", refld.file)) > 0) {
        # FIXME: filter by MAF
    } else {
        ld0.maf <- calcMAF.refLD.vcf(refgt, rep(TRUE, nrow(refgt)))

        refgt <- refgt[ld0.maf >= 0.05, ]
        ld0.maf <- ld0.maf[ld0.maf >= 0.05]
    }
    
    # exclude overlapping variants
    if (sum(duplicated(refgt$POS)) > 0) {
        dup.POS <- refgt$POS[duplicated(refgt$POS)]

        MSG("Removing multiple common alleles at same BP in refgt: BP =",
            paste(dup.POS, collapse=", "))
        
        ld0.maf <- ld0.maf[!(refgt$POS %in% dup.POS)]
        refgt <- refgt[!(refgt$POS %in% dup.POS), ]
    }

    if (sum(duplicated(assoc1$BP)) > 0) {
        dup.POS <- assoc1$BP[duplicated(assoc1$BP)]

        MSG("Removing multiple common alleles at same BP in trait 1: BP =",
            paste(dup.POS, collapse=", "))

        assoc1 <- assoc1[!(assoc1$BP %in% dup.POS), ]
    }

    if (sum(duplicated(assoc2$BP)) > 0) {        
#        PL("duplicated POS in tr2", sum(duplicated(assoc2$BP)))
#        permmat <- permmat[, !duplicated(assoc2$BP)]    
#        assoc2 <- assoc2[!duplicated(assoc2$BP), ]

        cat(paste("multiple alleles are prohibited at the same BP:",
                  sectr.file))
        stop()
    }

    # enforce the boundary of tested region
    ASSERT(start < end)

    refgt <- refgt[refgt$POS >= start & refgt$POS <= end, ]
    assoc1 <- assoc1[assoc1$BP >= start & assoc1$BP <= end, ]
    permmat <- permmat[, assoc2$BP >= start & assoc2$BP <= end]
    assoc2 <- assoc2[assoc2$BP >= start & assoc2$BP <= end, ]

    ASSERT(ncol(permmat) == nrow(assoc2))

    # common markers
    all.sel <- intersect(intersect(assoc1$BP, assoc2$BP), refgt$POS)
    all.sel <- setdiff(all.sel, assoc1$BP[is.na(assoc1$Z)])
    all.sel <- setdiff(all.sel, assoc2$BP[is.na(assoc2$Z)])

    refgt <- refgt[refgt$POS %in% all.sel, ]
    assoc1 <- assoc1[assoc1$BP %in% all.sel, ]
    permmat <- permmat[, assoc2$BP %in% all.sel]
    assoc2 <- assoc2[assoc2$BP %in% all.sel, ]
    
    ASSERT(nrow(refgt) == nrow(assoc1))
    ASSERT(ncol(permmat) == nrow(assoc2))
    ASSERT(nrow(assoc1) == nrow(assoc2))

    # sort assoc1 by BP
    assoc1 <- assoc1[order(assoc1$BP), ]

    # Don't sort assoc2 as it corresponds to permmat order
    # It should have been pre-ordered
    ASSERT(all(sort(assoc2$BP) == assoc2$BP))

    # Matching BP
    ASSERT(all(assoc1$BP == assoc2$BP))
    ASSERT(all(assoc1$BP == refgt$POS))


    ##########################################################################
    # Load LD (ref or in-sample LD)

    ###
    # ref LD

    if (length(grep(".hap$", refld.file)) > 0) {
        # refgt table is in .hap format
        ld0 <- process.refLD.impute(refgt, refgt$POS %in% assoc1$BP, NULL)

    } else {
        # refgt is from 1kg 
        ld0 <- process.refLD.vcf(refgt, refgt$POS %in% assoc1$BP, NULL)
    }

    ASSERT(nrow(ld0) == nrow(assoc1))
    ASSERT(nrow(ld0) == nrow(assoc2))
    ASSERT(sum(assoc2$BP != rownames(ld0)) == 0)


    ###
    # Secondary LD
    if (secld.file == ".") {

        # NOT RECOMMENDED TO USE
        ASSERT(FALSE)

#        ld2 <- ld0
    
    } else if (secld.file == "MPERM") {

        # NOT RECOMMENDED TO USE
        ASSERT(FALSE)
        
#        cat("Derive LD for secondary trait from mperm...\n")
    
#        pm <- permmat
#        pm <- pm ** 2
#        ld2 <- cor(pm)
#        ld2[ld2 < 0] <- 0
#        ld2 <- sqrt(ld2)

#        ASSERT(nrow(ld0) == nrow(ld2))

#        if (FALSE) {
#            plot(as.vector(abs(ld0)), as.vector(ld2), cex=0.5, xlim=c(0, 1),
#                 ylim=c(0, 1))
#            abline(0, 1, col="red")
#            abline(h=0.8, col="red")
#            abline(v=0.8, col="red")
#        }

    } else if ((length(grep("dosage.txt$", secld.file)) > 0) ||
               (length(grep("dosage.txt.gz$", secld.file)) > 0) ||
               (length(grep("dosage.gz$", secld.file)) > 0)) {
        # internal allelic dosage format
        
        cat("Derive LD for secondary trait from imputed dosages...\n")

        ld2 <- load.isLD1.dosage2(secld.file, all.sel)
        
        ASSERT(nrow(ld0) == nrow(ld2))    
    
    } else if ((length(grep(".ped.gz$", secld.file)) > 0) ||
               (length(grep(".ped$", secld.file)) > 0)) {

        cat("Derive LD for secondary trait from plink ped...\n")

        ld2 <- load.isLD1.ped(secld.file,
            1:nrow(assoc2.org) %in% as.numeric(rownames(assoc2)))

        ASSERT(nrow(ld0) == nrow(ld2))

        #plot(as.vector(abs(ld0)), as.vector(ld2), cex=0.5, xlim=c(0, 1),
        #     ylim=c(0, 1))
        #abline(0, 1, col="red")
        #abline(h=0.8, col="red")
        #abline(v=0.8, col="red")

    } else {
        PL("ERROR: Unsupported secld", secld.file)
        ASSERT(FALSE)
    }

    # MAIN LD
    # Use in-sample LD for main trait: power gain is limited
    if (mainld.file == ".") {

        # Default choice
        PL("Use Ref LD for main trait...\n")

    } else if ((length(grep("dosage.txt$", mainld.file)) > 0) ||
               (length(grep("dosage.txt.gz$", mainld.file)) > 0) ||
               (length(grep("dosage.gz$", mainld.file)) > 0)) {
        
        cat("Derive LD for main trait from imputed genotypes...\n")

        ld0 <- load.isLD1.dosage2(mainld.file, all.sel)
        ld0.maf <- NULL

        ASSERT(nrow(ld0) == nrow(ld1))

    } else if ((length(grep(".ped.gz$", mainld.file)) > 0) ||
               (length(grep(".ped$", mainld.file)) > 0)) {
    
        cat("Derive LD for main trait from plink ped...\n")

        ld1 <- load.isLD1.ped(mainld.file,
            1:nrow(assoc1.org) %in% as.numeric(rownames(assoc1)))

        ASSERT(nrow(ld0) == nrow(ld1))

        ld0 <- ld1
        ld0.maf <- NULL

    } else {
        PL("ERROR: Unsupported mainld", mainld.file)
        ASSERT(FALSE)
    }


    ##########################################################################
    # SNP selection step

    thresholdingP <- 0.1

    assoc1.sel <- assoc1$BP[assoc1$P <= thresholdingP]
    assoc2.sel <- assoc2$BP[assoc2$P <= thresholdingP]

    PL("Passoc-Threshold", thresholdingP)

    markers.t <- union(assoc1.sel, assoc2.sel)
    PL("Markers", length(markers.t))

    ASSERT(length(markers.t) > 1)
    ASSERT(assoc1$BP[which.min(assoc1$P)] %in% markers.t)

    ld0.t <- ld0[assoc1$BP %in% markers.t, assoc1$BP %in% markers.t]
    ld2.t <- ld2[assoc2$BP %in% markers.t, assoc2$BP %in% markers.t]

    assoc1.t <- assoc1[assoc1$BP %in% markers.t, ]
    assoc2.t <- assoc2[assoc2$BP %in% markers.t, ]

    if (!is.null(ld0.maf)) {
        ld0.maf.t <- ld0.maf[assoc1$BP %in% markers.t]

        if (sum(is.na(ld0.maf.t)) > 0) {
            cat(paste("NA in ref LD matrix:", refld.file))
            stop()
        }

        if (sum(ld0.maf.t == 0) > 0) {
            cat(paste("monomorphic SNPs in ref LD matrix:", refld.file, ": BP=",
                      paste(names(ld0.maf.t)[ld0.maf.t == 0], collapse=", ")))
            stop()
        }
    }
    
    if (sum(is.na(ld0.t)) > 0) {
        cat(paste("NA or monomorphic SNP in ref LD matrix:", refld.file))
        stop()
    }


    #########################################################################
    # Calculate test stats

    best1 <- which.max(abs(assoc1.t$Z))

    #if (idx2BP != assoc1.t$BP[best1]) {
    #    PL("Warning: idx2BP != most associated SNP for tr1",
    #         assoc1.t$BP[best1])
    #}


    lambda.t <- calc.stat(assoc1.t, assoc2.t, ld0.t, ld2.t, r2res)


    #########################################################################
    # Calculate permutation p

    NULLDIST <-
        perm.test(assoc1, assoc2, permmat,
                  ld0, ld2, thresholdingP, r2res, lambda.t)

    ASSERT(mean(!is.na(NULLDIST)) > 0.9)

    permP <- sum(NULLDIST >= lambda.t, na.rm=TRUE)/sum(!is.na(NULLDIST))

    return(new("jlim",
               assoc1=assoc1.t,
               assoc2=assoc2.t,
               ld0=ld0.t,
               ld2=ld2.t,
               idx2BP=assoc1$BP[best1],
               STAT=lambda.t,
               p.value=permP))
}
