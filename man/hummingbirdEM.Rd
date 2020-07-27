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
}
    \item{experimentInfoCase}{
A SummarizedExperiment object containing the input data for the case group:
The two assays: abnormM, abnormUM and the CpG position information: pos.
}
    \item{normM}{
A matrix containing the methylated read count data of the normal group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
}
    \item{normUM}{
A matrix containing the unmethylated read count data of the normal group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
}
    \item{abnormM}{
A matrix containing the methylated read count data of the abnormal group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
}
    \item{abnormUM}{
A matrix containing the unmethylated read count data of abnormal group. Each
column of the matrix represents a replicate and each row represents a CpG
position.
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
A GenomicRanges object that contains the start and end positions of each bin,
the distance between the current bin the bin ahead of it, the average
methylation rate of normal and abnormal groups and the predicted direction of
methylation change ("0" means a predicted normal bin; "1" means a predicted
hypermethylated bin; "2" means a predicted hypomethylated bin).
}
}
\examples{
library(GenomicRanges)
library(SummarizedExperiment)
data(exampleHummingbird)
emInfo <- hummingbirdEM(experimentInfoControl = exampleSEControl,
                        experimentInfoCase = exampleSECase, binSize = 40)
}

