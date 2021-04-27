catE<- function(input) {
  if(verbose)
     cat(paste("\n",input,sep=""))
}

fileExists<- function(fileName) {
  if (!file.exists(fileName) ){
    cat(paste("\nFile does not exist:", fileName))
    q(status=1)
  }
}

ASSERT <- function(test) {
    if (!test) {
        stop(paste("ASSERT fail:", deparse(substitute(test))))
    }
}

splitRefPanels <- function(refgt){

    gt<- refgt[,c(10:ncol(refgt))]
    gt2<-gt[, rep(1:ncol(gt), each=2)]
    for (i in 1:ncol(gt2) ){
        if(i%%2==0){
            gt2[,i] <- sapply(gt2[,i], function(x) strsplit(x,"|", fixed = TRUE)[[1]][2])
        }else{
            gt2[,i] <- sapply(gt2[,i], function(x) strsplit(x,"|", fixed = TRUE)[[1]][1])
        }
    }
    refgt_final<-  cbind.data.frame(refgt[,c(1:7)],gt2)

    return(refgt_final)
}

toOldRef <- function(refgt){
    refgt<-  refgt[,-c(8)]
    gt<- refgt[,c(10:ncol(refgt))]
    gt2<-gt[, rep(1:ncol(gt), each=2)]
    for (i in 1:ncol(gt2) ){
        if(i%%2==0){
            gt2[,i] <- sapply(gt2[,i], function(x) strsplit(x,"|", fixed = TRUE)[[1]][2])
        }else{
            gt2[,i] <- sapply(gt2[,i], function(x) strsplit(x,"|", fixed = TRUE)[[1]][1])
        }
    }
    for (i in 1:nrow(gt2) )
        gt2[i,] <- sapply(gt2[i,], function(x,R=refgt[i,4], A=refgt[i,5]) ifelse(x ==0, x<- R, x<- A))

    refgt_final<-  cbind.data.frame(refgt[,c(1:7)],gt2)

    return(refgt_final)
}

expTest <- function(x,y) {
    expr <- quote(x + y)
    print(expr)
    print(eval(expr) )
    print( quote(x))
    print (substitute(x - y))
    print( deparse(substitute(x - y)))
}

PL <- function(pre, mesg="") {
    if(!is.null(mesg))
      cat(paste("\n",pre, ": ", mesg, sep=""))
    else
      cat(paste("\n",pre, " is NA/Null.",sep=""))
}

MSG <- function(pre, mesg) {
    cat(paste("\n", pre, ": ", mesg,  sep=""))
}

Zscore <- function(assoc, N) {
    F_avg <- (assoc$F_A + assoc$F_U)/2
    z <- sqrt(N/2)*(assoc$F_A - assoc$F_U)/sqrt(F_avg*(1-F_avg))
    z[assoc$A1 == 1] <- -z[assoc$A1 == 1]
    z
}

# alleles ref/Alt = 1/2  -> 0/1
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
stdGT2 <- function(gt) {
    f <- apply(gt, 2, mean)/2
    for (I in 1:ncol(gt)) {
        gt[, I] <- (gt[, I] - 2*f[I])
    }

    #    sd0 <- apply(gt, 2, function(v) sqrt(mean(v*v)))
    sd0 <- apply(gt, 2, sd)
    for (I in 1:ncol(gt)) {
        gt[, I] <- gt[, I] / sd0[I]
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

log10PtoZ <- function(pv) {
    -qnorm(pv * log(10) - log(2), log.p=TRUE)
}
