\name{hummingbirdGraph}
\alias{hummingbirdGraph}
\title{
Methylation Level and Predictions Graphs
}
\description{
Generates the Methylation Level and Predictions graphs for the user-specified regions. The Methylation Level graph is the sample averages from two comparison groups. The Predictions graph is the sample average difference along with predictions, where "0" indicates no significant change, "1" indicates predicted hyper-methylation, and "2" indicates predicted hypo-methylation.
}
\usage{
hummingbirdGraph(pos, normM, normUM, abnormM, abnormUM, dmrs, coord1, coord2)
}
\arguments{
    \item{pos}{
The CpG positions. It must be in the structure of a matrix.
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
    \item{dmrs}{
The DMRs output of the hummingbirdPostAdjustment function.
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
hummingbirdGraph(pos=m_pos, normM=m_normM, normUM=m_normUM, abnormM=m_abnormM, abnormUM=m_abnormUM, dmrs=hmmbird2$DMRs, coord1=107991, coord2=108350)

}

