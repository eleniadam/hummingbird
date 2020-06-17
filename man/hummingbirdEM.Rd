\name{hummingbirdEM}
\alias{hummingbirdEM}
\title{
EM Algorithm for Fitting the Hidden Markov Model
}
\description{
This function reads in methylated and unmethylated read count data, transforms it into logarithm bin-wise data, sets up initial values and implements the EM algorithm to estimate HMM parameters and find the best sequence of hidden states based on model fitting.
}
\usage{
hummingbirdEM(experimentInfo, binSize)
}
\arguments{
    \item{experimentInfo}{
A SummarizedExperiment object containing the four assays: normM, normUM, abnormM, abnormUM and the CpG position information: pos.
}
    \item{normM}{
A matrix containing the methylated read count data of the normal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{normUM}{
A matrix containing the unmethylated read count data of the normal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{abnormM}{
A matrix containing the methylated read count data of the abnormal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{abnormUM}{
A matrix containing the unmethylated read count data of abnormal group. Each column of the matrix represents a replicate and each row represents a CpG position.
}
    \item{pos}{
The CpG positions. It must be in the structure of a vector.
}
    \item{binSize}{
The size of a bin. Default value is: 40.
}

}

\value{
\item{emInfo }{
A GenomicRanges object containing the start and end positions of each bin, the distance between the current bin the bin ahead of it, the average methylation rate of abnormal and normal groups, the initial direction, the transformed observations for each bin.
}
}
\examples{
library(GenomicRanges)
library(SummarizedExperiment)
data(exampleHummingbird)
emInfo <- hummingbirdEM(experimentInfo = exampleSE, binSize = 40)
}

