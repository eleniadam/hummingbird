\name{hummingbirdPostAdjustment}
\alias{hummingbirdPostAdjustment}
\title{
Post Adjustment algorithm for the output of the EM
}
\description{
This function adjusts HMM output. It enables three additional requirements on
DMRs: 1) the minimum length of a DMR, 2) the minimum number of CpGs in a DMR,
and 3) the maximum distance (in base pairs) between any two adjacent CpGs in a
DMR.
}
\usage{
hummingbirdPostAdjustment(experimentInfoControl, experimentInfoCase, emInfo,
                            minCpGs, minLength, maxGap)
}
\arguments{
    \item{experimentInfoControl}{
A SummarizedExperiment object containing the input data for the control group:
The two assays: normM, normUM and the CpG position information: pos.
}
    \item{experimentInfoCase}{
A SummarizedExperiment object containing the input data for the case group:
The two assays: abnormM, abnormUM and the CpG position information: pos.
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
\item{postAdjInfo }{ A list of two GenomicRanges objects, the DMRs and the
obsPostAdj.
Specifically, the DMRs contains the detected regions based on the
user-defined arguments (minLength, minCpGs, and maxGap). It contains the
refined DMRs with the start genomic position, the end genomic position, length
of the region, direction of predicted methylation change ("0" indicates no
significant change, "1" indicates predicted hyper-methylation, and "2"
indicates predicted hypo-methylation), and the number of CpGs. The obsPostAdj
object contains methylation status of each CpG site.
}
}

\examples{
library(GenomicRanges)
library(SummarizedExperiment)
data(exampleHummingbird)
emInfo <- hummingbirdEM(experimentInfoControl = exampleSEControl,
                        experimentInfoCase = exampleSECase, binSize = 40)
postAdjInfo <- hummingbirdPostAdjustment(
                                    experimentInfoControl = exampleSEControl,
                                    experimentInfoCase = exampleSECase,
                                    emInfo = emInfo, minCpGs = 10,
                                    minLength = 100, maxGap = 300)
}

