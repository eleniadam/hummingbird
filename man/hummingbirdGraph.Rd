\name{hummingbirdGraph}
\alias{hummingbirdGraph}
\title{
Observations and Predictions Graphs
}
\description{
This function generates observation and prediction graphs for a user specified
region.
}
\usage{
hummingbirdGraph(experimentInfoControl, experimentInfoCase, postAdjInfoDMRs,
                    coord1, coord2)
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
    \item{postAdjInfoDMRs}{
The DMRs GenomicRanges object output of the hummingbirdPostAdjustment function.
}
    \item{coord1}{
The start coordinate of the genomic region to plot.
}
    \item{coord2}{
The end coordinate of the genomic region to plot.
}
}
\value{
The function outputs two graphs: The Observations graph and the Predictions
graph.
The observation figure shows bin-wise average methylation rate for case and
control groups. The prediction figure shows bin-wise prediction, where "0"
denotes a predicted normal bin; "1" denotes a predicted hypermethylated bin;
and "2" denotes a predicted hypomethylated bin.
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
hummingbirdGraph(experimentInfoControl = exampleSEControl,
                    experimentInfoCase = exampleSECase,
                    postAdjInfoDMRs = postAdjInfo$DMRs,
                    coord1 = 107991, coord2 = 108350)
}

