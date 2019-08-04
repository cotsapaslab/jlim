
# Default window size
# +/- 100 kb
WINDOW.SIZE <- 100000

##' command-line interface to generate config file.
##'
##' This function scans through index SNPs listed in provided indexSNP file, 
##' and lists up all trait1-trait2 combinations to examine later by
##' jlim.test(). Typically, trait 1 is GWAS association and trait 2 is eQTLs.
##' Regarding the file formats, please see accompanied example test cases 
##' in the distribution.
##'
##' Usage:
##' \code{jlim_gencfg.sh}
##'
##' Command line arguments:
##' \itemize{
##' \item \code{--tr1-name}: trait 1 name (no white space allowed)
##' \item \code{--tr1-dir}: directory containing trait 1 association summary statistics
##' \item \code{--idxSNP-file}: indexSNP file
##' \item \code{--tr2-dir}: directory containing trait 2 association data
##' \item \code{--refld-dir}: reference LD panel (extracted by fecth.refld0.EUR.sh)
##' \item \code{--out}: output config file name
##' \item \code{--p-tr2-cutoff}: trait 2 association p-value cut-off to run the
##' analysis (default 0.01)
##' \item \code{--reassign-window}: add this option to readjust analysis window when provided index SNPs are not the most associated SNPs in trait 1 association summary data.
##' }
##'
##' \code{--tr1-name --tr1-dir --tr2-dir --idxSNP-file --refld-dir} are required
##' arguments, and the rest are optional.
##' @export
##' @import getopt
jlim.gencfg <- function () {
    
    spec = matrix(c(
        'tr1-name', 'n', 1, 'character',
        'tr1-dir', 'm', 1, 'character',
        'idxSNP-file', 'x', 1, 'character',
        'tr2-dir', 's', 1, 'character',
        'refld-dir', 'l', 1, 'character',

        'tr2-genotype-filetype', 'g', 1, 'character',
        'out', 'o', 1, 'character',
        'p-tr2-cutoff', 'p', 1, 'numeric',
        'maf-cutoff', 'f', 1, 'numeric'
#        , 'reassign-window', 'W', 0, 'logical'
        
        ), byrow=TRUE, ncol=4)

    args <- commandArgs(TRUE)

    ASSERT(args[1] == "ARGSTART")

    if (length(args) == 1) {

        cat(getopt(spec, usage=TRUE, command="jlim_gencfg.sh"))
        cat("Required options: --tr1-name --tr1-dir --tr2-dir --idxSNP-file --refld-dir\n")
        q(status=1)
    }

    args <- args[2:length(args)]

    opt <- getopt(spec, opt=args)

    if (is.null(opt[["tr1-dir"]]) || is.null(opt[["idxSNP-file"]]) ||
        is.null(opt[["tr2-dir"]]) || is.null(opt[["refld-dir"]]) ||
        is.null(opt[["tr1-name"]])) {

        cat(getopt(spec, usage=TRUE, command="jlim_gencfg.sh"))
        cat("Required options: --tr1-name --tr1-dir --tr2-dir --idxSNP-file --refld-dir\n")
        q(status=1)
    }

    tr1name <- opt[["tr1-name"]]
    tr1dir <- opt[["tr1-dir"]]
    idxSNPfile <- opt[["idxSNP-file"]]
    tr2dir <- opt[["tr2-dir"]]
    reflddir <- opt[["refld-dir"]]

    cfgfile <- opt[["out"]]
    minP2co <- opt[["p-tr2-cutoff"]]
    minMAF <- opt[["maf-cutoff"]]
    tr2gttype <- opt[["tr2-genotype-filetype"]]

    if (!file.exists(idxSNPfile)) {
        cat(paste("File does not exist:", idxSNPfile, "\n"))
        q(status=1)
    }

    if (!file.exists(tr1dir) || !file.info(tr1dir)$isdir) {
        cat(paste("Directory does not exist:", tr1dir, "\n"))
        q(status=1)
    }

    if (!file.exists(tr2dir) || !file.info(tr2dir)$isdir) {
        cat(paste("Directory does not exist:", tr2dir, "\n"))
        q(status=1)
    }

    if (!file.exists(reflddir) || !file.info(reflddir)$isdir) {
        cat(paste("Directory does not exist:", reflddir, "\n"))
        q(status=1)
    }

    # default
    if (is.null(cfgfile)) {
        cfgfile <- "jlim.cfg.tsv"
    }

    if (is.null(minP2co)) {
        minP2co <- 0.01
    }

    if (is.null(minMAF)) {
        minMAF <- 0.05
    }

    if (is.null(tr2gttype)) {
        tr2gttype <- "ped"
    }

    if (tr2gttype != "ped" && tr2gttype != "dosage") {
        cat(paste("tr2-genotype-filetype should be either ped or dosage:",
                  tr2gttype, "\n"))
        q(status=1)
    }

    # Print out set-up
    PL("trait 1 name", tr1name)    
    PL("trait 1 summary data directory", normalizePath(tr1dir))
    PL("trait 1 index SNP file", normalizePath(idxSNPfile))    
    PL("trait 2 data directory", normalizePath(tr2dir))
    PL("reference LD directory", normalizePath(reflddir))
    PL("output config file", cfgfile)
    
    #####
    # START
    idxSNPtab <- read.delim(idxSNPfile, header=TRUE, sep="\t", 
                            stringsAsFactors=FALSE)

    if (!all(c("CHR", "SNP", "BP", "STARTBP", "ENDBP") %in%
             colnames(idxSNPtab))) {

        idxSNPtab <- read.delim(idxSNPfile, header=TRUE, sep=" ",
                                na.strings='.', stringsAsFactors=FALSE)
    }

    if (!all(c("CHR", "SNP", "BP", "STARTBP", "ENDBP") %in%
             colnames(idxSNPtab))) {

        cat(paste("Invalid index SNP file:", idxSNPfile, "\n"))
        cat("required columns: CHR SNP BP STARTBP ENDBP\n")
    }


    df.idxsnp <- c()
    df.idxchr <- c()
    df.idxbp <- c()
    df.idx2bp <- c()
    df.idxP <- c()
    df.idx2P <- c()
    df.minP2 <- c()
    df.intvstart <- c()
    df.intvend <- c()
    df.mainassoc <- c()
    df.secassoc <- c()
    df.tr2name <- c()
    df.tr2subID <- c()    
    df.perm <- c()
    df.ld <- c()
    df.ld2 <- c()

    for (idx in 1:nrow(idxSNPtab)) {

        idxSNP <- idxSNPtab$SNP[idx]
        
        chrom <- idxSNPtab$CHR[idx]        
        idxSNP.bp <- idxSNPtab$BP[idx]

        # interval start to end        
        startpos <- idxSNPtab$STARTBP[idx]
        endpos <- idxSNPtab$ENDBP[idx]

        if (!is.numeric(startpos) || !is.numeric(endpos)) {
            cat(paste("Invalid interval for", idxSNP, "\n"))
            q(status=1)
        }
        
        ASSERT(startpos < endpos)

        # interval identifier
        locusname <- paste(chrom, startpos, endpos, sep=".")

        PL("Scanning", locusname)

        # reference LD
        refgt0.file <- paste(reflddir, "/locus.", locusname, ".txt.gz", sep='')
        refgt0 <-
            read.delim(file=refgt0.file,
                       header=FALSE, sep="\t", stringsAsFactors=FALSE)
        colnames(refgt0)[1:7] <-
            c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
        
        # exclude SVs
        if (sum(substr(refgt0$ID, 1, 3) == "esv") > 0) {
            refgt0 <- refgt0[substr(refgt0$ID, 1, 3) != "esv", ]
        }

        # filter rare from refgt
        ld0.maf <- calcMAF.refLD.vcf(refgt0, rep(TRUE, nrow(refgt0)))

        refgt0 <- refgt0[ld0.maf >= 0.05, ]
        ld0.maf <- ld0.maf[ld0.maf >= 0.05]
        
        # exclude overlapping variants
        if (sum(duplicated(refgt0$POS)) > 0) {
            dup.POS <- refgt0$POS[duplicated(refgt0$POS)]

            MSG("Removing multiple common alleles at same BP in refgt: BP =",
                paste(dup.POS, collapse=", "))
            
            ld0.maf <- ld0.maf[!(refgt0$POS %in% dup.POS)]
            refgt0 <- refgt0[!(refgt0$POS %in% dup.POS), ]
        }
        
        # eQTL gene list
        tr2subdir <- paste(tr2dir, "/locus.", locusname, sep='')

        if (!file.exists(tr2subdir) || !file.info(tr2subdir)$isdir) {
            cat(paste("directory does not exist:", tr2subdir, "\n"))
            q(status=1)
        }
        
        prefix2list <-
            list.files(tr2subdir,
                       paste("[^.]+.*.assoc.linear.gz", sep=''))

        #print(prefix2list)

        if (length(prefix2list) == 0) {
            PL("No secondary trait available", tr2subdir)
            next
        }

        genelist <- gsub(".assoc.linear.gz", "", prefix2list, fixed=TRUE)
        
        tr2name <-
            sapply(genelist,
                   function(str) strsplit(str, '.', fixed=TRUE)[[1]][1],
                   simplify=TRUE)
        
        genelist <-
            sapply(genelist,
                   function(str) {
                       v <- strsplit(str, ".", fixed=TRUE)[[1]]
                       paste(v[2:length(v)], collapse=".")
                   },
                   simplify=TRUE)

        # main trait path
        tr1path <- paste(tr1dir, "/", tr1name, '.', locusname, ".txt", sep="")
        
        assoc1 <- read.table(tr1path, header=TRUE, stringsAsFactors=FALSE)
        
        assoc1 <- assoc1[!is.na(assoc1$P), ]

        # Estimate Z statistics from P-values unless provided
        if (!("Z" %in% colnames(assoc1))) {
            assoc1 <- cbind(assoc1, Z=PtoZ(assoc1$P))

            if (sum(assoc1$P == 0) > 0) {
                cat(paste("P-value = 0 for", sum(assoc1$P == 0),
                          "SNP(s) in the window around",
                          idxSNP, ": check for potential numerical underflow\n"))
                q(status=1)
            }
        }
        
        dim(assoc1)
            
        # exclude overlapping variants
        if (sum(duplicated(assoc1$BP)) > 0) {
            dup.POS <- assoc1$BP[duplicated(assoc1$BP)]

            MSG("Removing multiple common alleles at same BP in trait 1: BP =",
                paste(dup.POS, collapse=", "))

            assoc1 <- assoc1[!(assoc1$BP %in% dup.POS), ]
        }
        
        idxP <- assoc1[assoc1$BP == idxSNP.bp, "P"]

        if (sum(assoc1$BP == idxSNP.bp) == 0) {
            #PL("Index SNP missing in main trait", idxSNP)
            idxP <- NA
        }
            
        assoc1 <- assoc1[assoc1$BP >= startpos & assoc1$BP <= endpos, ]
        ASSERT(nrow(assoc1) > 0)        
        assoc1.startbp <- min(assoc1$BP)
        assoc1.endbp <- max(assoc1$BP)
        
        Ign <- 1
        assoc1.org <- assoc1

        for (Ign in 1:length(prefix2list)) {

            eQTL <- tr2name[Ign]
            gene <- genelist[Ign]
                
            assoc1 <- assoc1.org

            # sec trait
            tr2file <- prefix2list[Ign]
            tr2path <- paste(tr2subdir, "/", tr2file, sep="")

            PL("Scanning secondary trait", tr2file)

            permfile <- paste(c(eQTL, gene, "mperm.dump.all.gz"), collapse='.')
            permpath <- paste(tr2subdir, "/", permfile, sep='')
            
            ld2path <- paste(tr2subdir, "/", eQTL, '.', tr2gttype, '.gz',
                             sep='')
                             
            
            assoc2 <-
                read.table(tr2path, header=TRUE, stringsAsFactors=FALSE)
            
            assoc2 <- assoc2[!is.na(assoc2$P), ]
            
            assoc2 <- assoc2[assoc2$BP >= assoc1.startbp &
                             assoc2$BP <= assoc1.endbp, ]
            
            if ("TEST" %in% colnames(assoc2)) {
                #PL("Remove NON-ADD TESTs: ", sum(assoc2$TEST != "ADD"))
                assoc2 <- assoc2[assoc2$TEST == "ADD", ]
            }

            if ("Z" %in% colnames(assoc2)) {
                # use provided Z
            } else if ("STAT" %in% colnames(assoc2)) {
                assoc2 <- cbind(assoc2, Z=assoc2$STAT)
            } else if ("T" %in% colnames(assoc2)) {
                assoc2 <- cbind(assoc2, Z=assoc2$T)
            } else {
                cat(paste("Cannot find association statistics in", tr2path))
                q(status=1)
            }

            if (sum(assoc2$BP == idxSNP.bp) != 1) {
                #PL("Index SNP is missing in assoc2", idxSNP)
            }

            minP2 <- min(assoc2$P)
                
            if (min(assoc2$P) > minP2co) {
                PL("min{assoc P of secondary trait} is high", min(assoc2$P))
                next
            }

            if (sum(duplicated(assoc2$BP)) > 0) {
                cat(paste("multiple alleles are prohibited at the same BP:",
                          tr2path))
                q(status=1)
            }

            assoc1.sel <- assoc1$BP[assoc1$P <= 0.01]
            assoc2.sel <- assoc2$BP[assoc2$P <= 0.01]

            all.sel <- union(assoc1.sel, assoc2.sel)
            all.sel <- intersect(all.sel, intersect(assoc1$BP, assoc2$BP))
            all.sel <- intersect(all.sel, refgt0$POS)
            all.sel <- intersect(all.sel, assoc1$BP[!is.na(assoc1$Z)])
            all.sel <- intersect(all.sel, assoc2$BP[!is.na(assoc2$Z)])    
            #PL("Markers", length(all.sel))

            if (length(all.sel) < 2) {
                PL("Too few usable markers", length(all.sel))
                next
            }

            assoc1 <- assoc1[assoc1$BP %in% all.sel, ]
            assoc2 <- assoc2[assoc2$BP %in% all.sel, ]

            assoc1 <- assoc1[order(assoc1$BP), ]
            assoc2 <- assoc2[order(assoc2$BP), ]        

            ASSERT(sum(assoc1$BP != assoc2$BP) == 0)

            refgt <- refgt0[refgt0$POS %in% all.sel, ]
            
            ASSERT(sum(assoc1$BP != refgt$POS) == 0)
             
            refgt.maf <- calcMAF.refLD.vcf(refgt, rep(TRUE, nrow(refgt)))

            # Reassign index SNP
            best1 <- which.max(abs(assoc1$Z))

            idx2bp <- assoc1$BP[best1]
            idx2P <- assoc1$P[best1]

            if (refgt.maf[best1] < minMAF) {
                PL(paste("index SNP ( BP =", idx2bp,
                         ") is below specified min MAF in Ref LD panel"),
                   refgt.maf[best1])
                next
            }

            # Reconfigure analysis window around new index SNP 
            assoc1 <- assoc1[(assoc1$BP >= (idx2bp - WINDOW.SIZE)) & 
                             (assoc1$BP <= (idx2bp + WINDOW.SIZE)), ]
            assoc2 <- assoc2[(assoc2$BP >= (idx2bp - WINDOW.SIZE)) &
                             (assoc2$BP <= (idx2bp + WINDOW.SIZE)), ]
            
            assoc1.startbp <- max(assoc1.startbp, idx2bp - WINDOW.SIZE)
            assoc1.endbp <- min(assoc1.endbp, idx2bp + WINDOW.SIZE)

            if (nrow(assoc1) < 2) {
                PL("Too few usable markers", nrow(assoc1))
                next
            }
            
            ASSERT(sum(assoc1$BP < assoc1.startbp |
                       assoc1$BP > assoc1.endbp)
                   == 0)
            ASSERT(sum(assoc2$BP < assoc1.startbp |
                       assoc2$BP > assoc1.endbp)
                   == 0)
            
            minP2 <- min(assoc2$P)
                
            if (minP2 > minP2co) {
                PL("min{assoc P of secondary trait} is high", min(assoc2$P))
                next
            }

            PL("Accept", tr2path)
            
            df.idxsnp <- append(df.idxsnp, idxSNP)
            df.idxchr <- append(df.idxchr, chrom)
            df.idxbp <- append(df.idxbp, idxSNP.bp)
            df.idx2bp <- append(df.idx2bp, idx2bp)
            df.idxP <- append(df.idxP, idxP)
            df.idx2P <- append(df.idx2P, idx2P)
            df.tr2name <- append(df.tr2name, eQTL)
            df.tr2subID <- append(df.tr2subID, gene)
            df.minP2 <- append(df.minP2, minP2)
            df.intvstart <- append(df.intvstart, assoc1.startbp)
            df.intvend <- append(df.intvend, assoc1.endbp)
            df.mainassoc <- append(df.mainassoc, tr1path)
            df.secassoc <- append(df.secassoc, tr2path)
            df.perm <- append(df.perm, permpath)
            df.ld <- append(df.ld, refgt0.file)
            df.ld2 <- append(df.ld2, ld2path)
        }

    }
    
    cfg <- 
        data.frame(maintrID=tr1name,
                   chrom=df.idxchr, idxSNP=df.idxsnp,
                   idxBP=df.idxbp, idxP=df.idxP,
                   idx2BP=df.idx2bp, idx2P=df.idx2P,
                   start=df.intvstart, end=df.intvend,                   
                   sectrID=df.tr2name, sectrsubID=df.tr2subID, 
                   minP2=df.minP2,
                   maintr=df.mainassoc, sectr=df.secassoc,
                   refld=df.ld, mainld=".", secld=df.ld2,
                   perm=df.perm,
                   stringsAsFactors=FALSE)

    write.table(cfg, file=cfgfile,
                quote=FALSE, sep="\t", row.names=FALSE,
                col.names=TRUE)

    q(status=0)
}
