\name{hummingbirdPostAdjustmentInternal}
\alias{hummingbirdPostAdjustmentInternal}
\title{
Post Adjustment algorithm (internal)
}
\description{
Post Adjustment algorithm for the output of the EM.
This function adjusts HMM output such that each detected DMR
has a minimum length and a minimum number of CpGs in each DMR.
}
\usage{
hummingbirdPostAdjustmentInternal(em, pos, minCpGs, minLength, maxGap)
}
\arguments{
    \item{em}{
The output of the hummingbirdEMinternal function, specifically the obs 
object.
}
    \item{pos}{
The CpG position information.
}
    \item{minCpGs}{
The minimum number of CpGs contained in a DMR.
}
    \item{minLength}{
The minimum length of a DMR.
}
    \item{maxGap}{
The maximum gap between any two CpGs.
}
}
\details{
Users do not need to call this function directly.
This is a low-level function used by the higher-level function in the 
hummingbird package, the hummingbirdPostAdjustment.
}
\value{
\item{DMRs }{The detected regions based on the user-defined arguments 
(minLength, minCpGs, and maxGap). 

It contains the (numbered) refined DMRs with the 
start genomic position, the end genomic position, length of the region, 
direction of predicted methylation change ("0" indicates no significant 
change, "1" indicates predicted hyper-methylation, and "2" indicates 
predicted hypo-methylation) and the number of CpGs.
}
\item{obsPostAdj }{The methylation status of each CpG site.
}
}

\seealso{
Users may call the \code{\link{hummingbirdPostAdjustment}} function.
}
\examples{
library(GenomicRanges)
library(SummarizedExperiment)

# Load sample dataset containing input data 
data(exampleHummingbird)

# Run the EM (internal) function
hmmbirdEMinternalOutput <- hummingbirdEMinternal(
    normM = assays(exampleSEControl)[["normM"]],
    normUM = assays(exampleSEControl)[["normUM"]],
    abnormM = assays(exampleSECase)[["abnormM"]],
    abnormUM = assays(exampleSECase)[["abnormUM"]],
    pos = pos, binSize = 40)
                    
# Run the Post Adjustment (internal) function
hmmbirdPAinternalOutput <- hummingbirdPostAdjustmentInternal(
    em = hmmbirdEMinternalOutput$obs, 
    pos = pos, minCpGs = 10, minLength = 100, maxGap = 300)
}
