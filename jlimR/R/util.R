ASSERT <- function(test) {
    if (!test) {
        stop(paste("ASSERT fail:", deparse(substitute(test))))
    }
}

PL <- function(pre, mesg="") {
    cat(paste(pre, ": ", mesg, "\n", sep=""))
}

MSG <- function(pre, mesg) {
    #cat(paste(pre, ": ", mesg, "\n", sep=""))
}

Zscore <- function(assoc, N) {
    F_avg <- (assoc$F_A + assoc$F_U)/2
    z <- sqrt(N/2)*(assoc$F_A - assoc$F_U)/sqrt(F_avg*(1-F_avg))
    z[assoc$A1 == 1] <- -z[assoc$A1 == 1]
    z
}               

toGT <- function(gt) {
    M <- ncol(gt)/2
    gt.alt <- (gt != 1) + 0 # nonref
    matrix(apply(gt.alt, 1, function(v) (v[seq(1, 2*M, by=2)] + v[seq(2, 2*M, by=2)])), 
           byrow=TRUE, ncol=M)
}

stdGT <- function(gt) {
    f <- apply(gt, 2, mean)/2
    for (I in 1:ncol(gt)) {
        gt[, I] <- (gt[, I] - 2*f[I]) / sqrt(2*f[I]*(1-f[I]))
    }
    gt
}

read.GT0 <- function(filepath) {
    ped <- read.table(filepath, header=FALSE, stringsAsFactors=FALSE)

    samples <- 1:nrow(ped)

    gt <- ped[samples, 7:ncol(ped)]
    gt <- as.matrix(gt)

    # missing genotype
    PL("GT 0", sum(gt == 0))    
    gt[gt == 0] <- NA         
    
    gt <- toGT(gt)
    PL("GT DIM", dim(gt))

    gt
}

to.refal <- function (assoc) {
    flip <- assoc$A1 == 1
    
    assoc[flip, "F_A"] <- 1 - assoc[flip, "F_A"]
    assoc[flip, "F_U"] <- 1 - assoc[flip, "F_U"]
    assoc[flip, "OR"] <- 1 / assoc[flip, "OR"]
    assoc[flip, "Z"] <- -assoc[flip, "Z"]

    assoc
}

# without OR info
PtoZ <- function(pv) {
    -qnorm(pv/2)
}
