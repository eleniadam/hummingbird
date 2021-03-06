hummingbirdEM <- function(experimentInfoControl, experimentInfoCase,
                        binSize = 40){

    ## Check Input
    if(dim(experimentInfoControl)[2] < 1){
        print("Please enter data from at least one replicate.")
        return(0)
    }

    ## Computation
    hmmbirdEM <- hummingbirdEMinternal(
        assays(experimentInfoControl)[["normM"]],
        assays(experimentInfoControl)[["normUM"]],
        assays(experimentInfoCase)[["abnormM"]],
        assays(experimentInfoCase)[["abnormUM"]],
        matrix(rowRanges(experimentInfoControl)@ranges@pos), binSize)

    ## Output
    outEM <- GRanges(matrix(rowRanges(experimentInfoControl)@seqnames@values),
                    IRanges(hmmbirdEM$obs[,3], hmmbirdEM$obs[,4]),
                    distance = hmmbirdEM$obs[,2],
                    norm = hmmbirdEM$normAbnorm[,1],
                    abnorm = hmmbirdEM$normAbnorm[,2],
                    direction = hmmbirdEM$obs[,1])

    return(outEM)
}


