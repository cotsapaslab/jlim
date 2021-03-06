% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jlim.R
\docType{class}
\name{jlim-class}
\alias{jlim-class}
\title{jlim class to keep the result of JLIM single test}
\description{
jlim class to keep the result of JLIM single test
}
\section{Slots}{

\describe{
\item{\code{userIdxBP}}{user specified index SNP}

\item{\code{actualIdxBP}}{actual found index SNP}

\item{\code{STAT}}{jlim statistic lambda}

\item{\code{pvalue}}{permutation pvalue}

\item{\code{startBP}}{start position of the tested locus}

\item{\code{endBP}}{end position of the tested locus}

\item{\code{sectrSampleSize}}{end position of the tested locus}

\item{\code{sectrGeneName}}{name of the Gene in case of multiple gene in the second trait}

\item{\code{sectrIndSNPpvalue}}{pvalue of the indexSNP in the second trait.}

\item{\code{sectrMinpvalue}}{minimum pvalue of in the second trait.}

\item{\code{sectrSNPWithMinpvalue}}{SNP with the minimum pvalue in the second trait}

\item{\code{desc}}{status of the JLIM test}

\item{\code{executedPerm}}{number of the executed permutations}

\item{\code{permmat}}{the permutation matrix used for the JLIM test}
}}

