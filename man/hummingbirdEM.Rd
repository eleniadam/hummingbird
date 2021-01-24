\name{hummingbirdEM}
\alias{hummingbirdEM}
\title{
EM Algorithm for Fitting the Hidden Markov Model
}
\description{
This function reads input data, sets initial values, executes the
Expectation-Maximization (EM) algorithm for the Bayesian HMM and infers the
best sequence of methylation states.
}
\usage{
hummingbirdEM(experimentInfoControl, experimentInfoCase, binSize)
}
\arguments{
    \item{experimentInfoControl}{
A SummarizedExperiment object containing the input data for the control group:
The two assays: normM, normUM and the CpG position information: pos. 

normM is a matrix containing the methylated read count data of the normal 
group and normUM is a matrix containing the unmethylated read count data of 
the normal group. Each column of a matrix represents a replicate and each 
row represents a CpG position.
}
    \item{experimentInfoCase}{
A SummarizedExperiment object containing the input data for the case group:
The two assays: abnormM, abnormUM and the CpG position information: pos. 

abnormM is a matrix containing the methylated read count data of the abnormal 
group and abnormUM is a matrix containing the unmethylated read count data of 
the abnormal group. Each column of a matrix represents a replicate and each 
row represents a CpG position.
}
    \item{binSize}{
The size of a bin. Default value is: 40.
}

}

\value{
A GenomicRanges object that contains the start and end positions of each bin,
the distance between the current bin the bin ahead of it, the average
methylation rate of normal and abnormal groups and the predicted direction of
methylation change ("0" means a predicted normal bin; "1" means a predicted
hypermethylated bin; "2" means a predicted hypomethylated bin).
}
\examples{
library(GenomicRanges)
library(SummarizedExperiment)
data(exampleHummingbird)

# CpG position vector
pos <- pos[,1]
# Assays for the normal group
assaysControl <- list(normM = normM, normUM = normUM)
# Assays for the case group
assaysCase <- list(abnormM = abnormM, abnormUM = abnormUM)
# SummarizedExperiment objects
exampleSEControl <- SummarizedExperiment(assaysControl, 
                                        rowRanges = GPos("chr29", pos))
exampleSECase <- SummarizedExperiment(assaysCase, 
                                    rowRanges = GPos("chr29", pos))

emInfo <- hummingbirdEM(experimentInfoControl = exampleSEControl,
                        experimentInfoCase = exampleSECase, binSize = 40)
}

