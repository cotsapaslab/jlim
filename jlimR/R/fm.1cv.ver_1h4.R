
load.isLD1.flag <- function(file1, snpset.flag) {
    
    # cohort 1 
    gt1 <- read.GT0(file1)
    gt1 <- gt1[, snpset.flag]

    ASSERT(ncol(gt1) == sum(snpset.flag))
        
    # IN-SAMPLE LD        
    gt0 <- gt1
    gt0 <- stdGT(gt0)

    PL("mean gt0", range(apply(gt0, 2, mean)))
    PL("sd gt0", range(apply(gt0, 2, sd)))

    ld0 <- (t(gt0) %*% gt0)
    var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v)))
    for (I in 1:nrow(ld0)) {
        for (J in 1:nrow(ld0)) {
            ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
        }
    }

    ld0
}

load.isLD1.ped <- function(file1, snpset.flag) {
    ped <- read.table(file1, header=FALSE, stringsAsFactors=FALSE,
                      colClasses="character")

    gt <- ped[, 7:ncol(ped)]
    gt <- as.matrix(gt)

    # Make sure all sites are biallelic
    al1 <- gt[, seq(1, ncol(gt), by=2)]
    al2 <- gt[, seq(2, ncol(gt), by=2)]

    # select SNPs
    al1 <- al1[, snpset.flag]
    al2 <- al2[, snpset.flag]

    al <- rbind(al1, al2)
    
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

    # IN-SAMPLE LD        
    gt0 <- gt.recode
    gt0 <- stdGT(gt0)

    PL("mean gt0", range(apply(gt0, 2, mean)))
    PL("sd gt0", range(apply(gt0, 2, sd)))

    ld0 <- (t(gt0) %*% gt0)
    var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v)))
    for (I in 1:nrow(ld0)) {
        for (J in 1:nrow(ld0)) {
            ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
        }
    }

    ld0
}

load.isLD1.ped2 <- function(file1, snpset.flag) {
    ped <- read.table(file1, header=FALSE, stringsAsFactors=FALSE)

    gt <- ped[, 7:ncol(ped)]
    gt <- as.matrix(gt)

    al1 <- gt[, seq(1, ncol(gt), by=2)]
    al2 <- gt[, seq(2, ncol(gt), by=2)]

    # select SNPs
    al1 <- al1[, snpset.flag]
    al2 <- al2[, snpset.flag]

    al <- rbind(al1, al2)

    # Make sure all sites are biallelic
    alcount <- apply(al, 2, function (X) length(setdiff(unique(X), 0)))
    ASSERT(sum(alcount > 2) == 0)

    # Recode
    # missing gt is NA.
    al.recode <- 
        matrix(apply(al, 2, function (X) as.numeric(X)),
               byrow=FALSE, nrow=nrow(al), ncol=ncol(al))
    al.recode[al == 0] <- NA
    
    ASSERT(sum(!is.na(al.recode) & al.recode != 1 &
               al.recode != 2) == 0)

    # To genotype
    N <- nrow(al.recode) / 2
    gt.recode <-
        matrix(apply(al.recode, 2,
                     function (X) (as.numeric(X[1:N] == 2) +
                                   as.numeric(X[(N+1):(2*N)] == 2))),
               byrow=FALSE, nrow=N, ncol=ncol(al))

    ASSERT(sum(!is.na(gt.recode) & gt.recode != 0 & gt.recode != 1 &
               gt.recode != 2) == 0)

    # IN-SAMPLE LD        
    gt0 <- gt.recode
    gt0 <- stdGT(gt0)

    PL("mean gt0", range(apply(gt0, 2, mean, na.rm=TRUE)))
    PL("sd gt0", range(apply(gt0, 2, sd, na.rm=TRUE)))

    ld0 <- cor(gt0, use="pairwise.complete.obs")

    ld0
}

get.dosage.altaf <- function(dosage.file, pheno.file, bpset) {
    # read in samples IDs (pheno1$fam)
    pheno1 <- read.delim(pheno.file, header=FALSE, stringsAsFactors=FALSE,
                         sep=" ")
    dim(pheno1)
    colnames(pheno1) <- c("fam", "ind", "expr")
    
    dosage1 <- read.delim(dosage.file, header=TRUE, stringsAsFactors=FALSE,
                          sep="\t")
    dim(dosage1)
    gt1 <- dosage1[, 5:ncol(dosage1)]
    rownames(gt1) <- dosage1$clone # SNP
    
    # reorder SNPs        
    gt1 <- gt1[dosage1$Start %in% bpset, ]
    dim(gt1)
    gt1 <- gt1[, colnames(gt1) %in% pheno1$fam]
    gt1 <- as.matrix(gt1)
    dim(gt1)

    meanF0 <- 1 - apply(gt1, 1, mean)/2
    meanF0
}    

load.isLD1.dosage <- function(dosage.file, pheno.file, bpset) {
    # read in samples IDs (pheno1$fam)
    pheno1 <- read.delim(pheno.file, header=FALSE, stringsAsFactors=FALSE,
                         sep=" ")
    dim(pheno1)
    colnames(pheno1) <- c("fam", "ind", "expr")
    
    dosage1 <- read.delim(dosage.file, header=TRUE, stringsAsFactors=FALSE,
                          sep="\t")
    dim(dosage1)
    gt1 <- dosage1[, 5:ncol(dosage1)]
    rownames(gt1) <- dosage1$clone # SNP
    
    # reorder SNPs        
    gt1 <- gt1[dosage1$Start %in% bpset, ]
    dim(gt1)
    gt1 <- gt1[, colnames(gt1) %in% pheno1$fam]
    gt1 <- as.matrix(gt1)
    gt1 <- t(gt1)
    dim(gt1)

    
    # IN-SAMPLE LD        
    gt0 <- gt1
    gt0 <- stdGT(gt0)

    PL("mean gt0", range(apply(gt0, 2, mean)))
    PL("sd gt0", range(apply(gt0, 2, sd)))

    ld0 <- (t(gt0) %*% gt0)
    var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v)))
    for (I in 1:nrow(ld0)) {
        for (J in 1:nrow(ld0)) {
            ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
        }
    }

    ld0
}

load.isLD1.dosage2 <- function(dosage.file, bpset) {
    
    dosage1 <- read.delim(dosage.file, header=TRUE, stringsAsFactors=FALSE,
                          sep="\t")
    dim(dosage1)
    gt1 <- dosage1[, 5:ncol(dosage1)]
    rownames(gt1) <- dosage1$clone # SNP
    
    gt1 <- gt1[dosage1$Start %in% bpset, ]
    dim(gt1)
    gt1 <- as.matrix(gt1)
    gt1 <- t(gt1)
    dim(gt1)

    # IN-SAMPLE LD        
    gt0 <- gt1
    gt0 <- stdGT(gt0)

    PL("mean gt0", range(apply(gt0, 2, mean)))
    PL("sd gt0", range(apply(gt0, 2, sd)))

    ld0 <- (t(gt0) %*% gt0)
    var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v)))
    for (I in 1:nrow(ld0)) {
        for (J in 1:nrow(ld0)) {
            ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
        }
    }

    ld0
}

process.refLD.impute <- function(refgt, refgt.all.sel, af.check=NULL) {

    refgt0 <- refgt[refgt.all.sel, ]

    if ("BP" %in% colnames(refgt0)) {
        refgt0 <- refgt0[order(refgt0$BP), ]
    } else if ("POS" %in% colnames(refgt0)) {
        refgt0 <- refgt0[order(refgt0$POS), ]
    }
    
    ASSERT(sum(refgt.all.sel) == nrow(refgt0))

    PL("RefLD cols", paste(paste(colnames(refgt0)[1:4], collapse=" "),
                           colnames(refgt0)[5], sep=" / "))

    refgt.mat <- refgt0[, 5:ncol(refgt0)]
    refgt.mat <- as.matrix(refgt.mat)
    
    ASSERT(sum(refgt.mat != 0 & refgt.mat != 1) == 0)

    # CONVERT
    # ref allele := 0
    # alt allele := 1
    # TO
    # ref allele := 1
    # alt allele := 2
    refgt.mat <- refgt.mat + 1
    refgt.mat <- toGT(refgt.mat)
    rownames(refgt.mat) <- refgt0$BP

    if (!is.null(af.check)) {
        meanF0 <- apply(refgt.mat, 1, mean)/2
        plot(meanF0, af.check, xlab="Alt AF (RefGT)",
             ylab="Alt AF (In-sample)")
    }

    gt0 <- stdGT(t(refgt.mat))
    
    ld0 <- (t(gt0) %*% gt0)
    var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v)))
    for (I in 1:nrow(ld0)) {
        for (J in 1:nrow(ld0)) {
            ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
        }
    }

    ASSERT(sum(refgt.all.sel) == nrow(ld0)) 

    ld0
}

process.refLD.vcf <- function(refgt, refgt.all.sel, af.check=NULL) {

    refgt0 <- refgt[refgt.all.sel, ]
    refgt0 <- refgt0[order(refgt0$POS), ]
    
    ASSERT(sum(refgt.all.sel) == nrow(refgt0))

    refgt.mat <- refgt0[, 8:ncol(refgt0)]
    refgt.mat <- as.matrix(refgt.mat)
    
    ASSERT(sum(refgt.mat == ".") == 0)

    # ref allele := 1
    refgt.mat <-
        t(sapply(1:nrow(refgt.mat),
                 function(I) as.numeric(refgt.mat[I, ] == refgt0$REF[I]),
                 simplify=TRUE))
    # alt allele := 2
    refgt.mat[refgt.mat != 1] <- 2
    refgt.mat <- toGT(refgt.mat)
    rownames(refgt.mat) <- refgt0$POS
    
    if (!is.null(af.check)) {
#        meanF0 <- apply(refgt.mat, 1, mean)/2
#        plot(meanF0, af.check, xlab="Alt AF (RefGT)",
#             ylab="Alt AF (In-sample)")

        meanF0 <- apply(refgt.mat, 1, mean)/2
        meanF0 <- pmin(meanF0, 1 - meanF0)
        plot(meanF0, af.check, xlab="MAF (RefLD)",
             ylab="MAF (In-sample)")

        af.check.tab <- cbind(refgt0[, 1:3], meanF0, af.check)
        af.check.tab <-
            af.check.tab[order(abs(meanF0 - af.check), decreasing=TRUE), ]
        print(af.check.tab[1:10, ])
    }

    gt0 <- stdGT(t(refgt.mat))
    
    ld0 <- (t(gt0) %*% gt0)
    var0 <- apply(gt0, 2, function(v) sqrt(sum(v*v)))
    for (I in 1:nrow(ld0)) {
        for (J in 1:nrow(ld0)) {
            ld0[I, J] <- ld0[I, J]/var0[I]/var0[J]
        }
    }

    ASSERT(sum(refgt.all.sel) == nrow(ld0))     

    ld0
}

calcMAF.refLD.vcf <- function(refgt, refgt.all.sel) {

    refgt0 <- refgt[refgt.all.sel, ]
    refgt0 <- refgt0[order(refgt0$POS), ]
    
    ASSERT(sum(refgt.all.sel) == nrow(refgt0))

    refgt.mat <- refgt0[, 8:ncol(refgt0)]
    refgt.mat <- as.matrix(refgt.mat)
    
    ASSERT(sum(refgt.mat == ".") == 0)

    # ref allele := 1
    refgt.mat <-
        t(sapply(1:nrow(refgt.mat),
                 function(I) as.numeric(refgt.mat[I, ] == refgt0$REF[I]),
                 simplify=TRUE))
    # alt allele := 2
    refgt.mat[refgt.mat != 1] <- 2
    refgt.mat <- toGT(refgt.mat)
    rownames(refgt.mat) <- refgt0$POS
    
#        meanF0 <- apply(refgt.mat, 1, mean)/2
#        plot(meanF0, af.check, xlab="Alt AF (RefGT)",
#             ylab="Alt AF (In-sample)")

    meanF0 <- apply(refgt.mat, 1, mean)/2
    meanF0 <- pmin(meanF0, 1 - meanF0)
#        plot(meanF0, af.check, xlab="MAF (RefLD)",
#             ylab="MAF (In-sample)")

    meanF0
}

load.mperm <- function (permfile) {
    permtab <- read.delim(permfile, sep=" ", header=FALSE,
                          stringsAsFactors=FALSE)
    # perm row 0
    permtab <- permtab[permtab[, 1] >= 1, ]

    # last column NA
    if (sum(!is.na(permtab[, ncol(permtab)])) == 0) {
        permtab <- permtab[, 1:(ncol(permtab)-1)]
    }

    # remove run no
    as.matrix(permtab[, 2:ncol(permtab)])
}

calc.stat <- function (assoc1, assoc2, ld0, ld2, R2thr) {

    # FIXED for 1h4
    logP1 <- (abs(assoc1$Z) ** 2) / 2
    logP2 <- (abs(assoc2$Z) ** 2) / 2

    ASSERT(sum(is.na(assoc1$Z)) == 0)
    ASSERT(sum(is.na(assoc2$Z)) == 0)

    ### PEAK SELECTION (local)
    maxI1 <- which.max(logP1)
    relP1 <- exp(logP1 - max(logP1))
    postP1 <- exp(logP1) / sum(exp(logP1))
    local <- which(ld0[maxI1, ]**2 >= R2thr)

    cat(paste("nonpeakPP1.ld0: ", sum(postP1[-local]), "\n", sep=''))

    gap <- 0
    sumRelP1 <- sum(relP1[local])
    
    for (I in local) {
        gap <- gap +
            relP1[I] * (logP2[I] - max(logP2[ld2[I, ]**2 < R2thr]))
        
        cat(paste("*", 
                  relP1[I],
                  postP1[I],
                  ld0[maxI1, I]**2,
                  logP2[I] - max(logP2[ld2[I, ]**2 < R2thr]),
                  "Farthest:",
                  max(ld2[I, (ld2[I, ]**2 < R2thr)]**2), 
                  "\n")
              )
    }

    gap <- gap / sumRelP1
    gap        

    PL("--GAP", gap)        

    gap
}

# FOR SIMULATED NULL DIST
perm.test <- function (assoc1, assoc2, permmat, ld0, ld2,
                       thresholdingP, R2thr, lambda.t) {

    thresholdingZ <- PtoZ(thresholdingP)

    # FIXED in version 1h4    
    logP1 <- (abs(assoc1$Z) ** 2) / 2

    # SNP selection part 1
    markers.p1 <- abs(assoc1$Z) >= thresholdingZ
    
    # peak selection
    maxlogP1 <- max(logP1)
    maxI1 <- which.max(logP1)
    relP1 <- exp(logP1 - maxlogP1)
    
    # iterate
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

            # FIXED in 1h4        
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

