\name{hummingbirdEMinternal}
\alias{hummingbirdEMinternal}
\title{
EM Algorithm (internal)
}
\description{
Expectation-Maximization Algorithm for Fitting the Hidden Markov Model.
This function reads in methylated and unmethylated read count data,
transforms it into logarithm bin-wise data, sets up initial values
and implements the EM algorithm to estimate HMM parameters and find
the best sequence of hidden states based on model fitting.
}
\usage{
hummingbirdEMinternal(normM, normUM, abnormM, abnormUM, pos, binSize)
}

\arguments{
    \item{normM}{
A matrix containing the methylated read count data of the normal group.
Each column of a matrix represents a replicate and each row represents 
a CpG position.
}
    \item{normUM}{
A matrix containing the unmethylated read count data of the normal group.
Each column of a matrix represents a replicate and each row represents 
a CpG position.
}
    \item{abnormM}{
A matrix containing the methylated read count data of the abnormal group.
Each column of a matrix represents a replicate and each row represents 
a CpG position.
}
    \item{abnormUM}{
A matrix containing the unmethylated read count data of the abnormal group.
Each column of a matrix represents a replicate and each row represents 
a CpG position.
}
    \item{pos}{
The CpG position information.
}
    \item{binSize}{
The size of a bin.
}
}
\details{
Users do not need to call this function directly.
This is a low-level function used by the higher-level function in the 
hummingbird package, the hummingbirdEM.
}
\value{
\item{obs }{For each bin: 
The predicted direction of methylation change ("0" means a predicted normal 
bin; "1" means a predicted hypermethylated bin; "-1" means a predicted 
hypomethylated bin).
The distance between the current bin and the bin ahead of it,
the start and end positions of each bin.
}
\item{normAbnorm }{The average methylation rate of normal and 
abnormal groups.}
}

\seealso{
Users may call the \code{\link{hummingbirdEM}} function.
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
}
