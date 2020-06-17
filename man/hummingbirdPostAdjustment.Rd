\name{hummingbirdPostAdjustment}
\alias{hummingbirdPostAdjustment}
\title{
Post Adjustment algorithm for the output of the EM
}
\description{
This function adjusts HMM output such that each detected DMR has a minimum length and a minimum number of CpGs in each DMR.
}
\usage{
hummingbirdPostAdjustment(experimentInfo, emInfo, minCpGs, minLength, maxGap)
}
\arguments{
    \item{experimentInfo}{
A SummarizedExperiment object containing the four assays: normM, normUM, abnormM, abnormUM and the CpG position information: pos.
}
    \item{emInfo}{
The output GenomicRanges object of the hummingbirdEM function.
}
    \item{minCpGs}{
The minimum number of CpGs contained in a DMR. Default value: 10.
}
    \item{minLength}{
The minimum length of a DMR. Default value: 500.
}
    \item{maxGap}{
The maximum gap between any two CpGs. Default value: 300.
}
}
\value{
\item{postAdjInfo }{ A list of two GenomicRanges objects, the obsPostAdj and the DMRs.
The obsPostAdj contains the transformed observations for each bin and the adjusted predictions for each bin (0-normal, 1-hyper, 2-hypo).
The DMRs contains the detected regions based on the user-defined arguments, providing the start genomic position, the end genomic position, the length, the direction and the number of CpGs.
}
}

\examples{
library(GenomicRanges)
library(SummarizedExperiment)
data(exampleHummingbird)
emInfo <- hummingbirdEM(experimentInfo = exampleSE, binSize = 40)
postAdjInfo <- hummingbirdPostAdjustment(experimentInfo = exampleSE, emInfo = emInfo, minCpGs = 10, minLength = 100, maxGap = 300)
}

