\name{hummingbirdGraph}
\alias{hummingbirdGraph}
\title{
Methylation Level and Predictions Graphs
}
\description{
Generates the Methylation Level and Predictions graphs for the user-specified regions. The Methylation Level graph is the sample averages from two comparison groups. The Predictions graph is the sample average difference along with predictions, where "0" indicates no significant change, "1" indicates predicted hyper-methylation, and "2" indicates predicted hypo-methylation.
}
\usage{
hummingbirdGraph(experimentInfo, postAdjInfoDMRs, coord1, coord2)
}
\arguments{
    \item{experimentInfo}{
A SummarizedExperiment object containing the four assays: normM, normUM, abnormM, abnormUM and the CpG position information: pos.
}
    \item{postAdjInfoDMRs}{
The DMRs GenomicRanges object output of the hummingbirdPostAdjustment function.
}
    \item{coord1}{
The start coordinate of the genomic region to look into.
}
    \item{coord2}{
The end coordinate of the genomic region to look into.
}
}
\value{
The function outputs two graphs: The Methylation Level graph and the Predictions graph.
}
\examples{
hummingbirdGraph(experimentInfo = exampleSE, postAdjInfoDMRs = postAdjInfo$DMRs, coord1 = 107991, coord2 = 108350)

}

