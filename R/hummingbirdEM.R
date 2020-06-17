hummingbirdEM <- function(experimentInfo, binSize = 40){

    ## Check Input
    if(dim(experimentInfo)[2] < 1){
        print("Please enter data from at least one replicate.")
        return(0)
    }

    ## Computation
    hmmbirdEM <- hummingbirdEMinternal(assays(experimentInfo)[["normM"]], assays(experimentInfo)[["normUM"]], assays(experimentInfo)[["abnormM"]], assays(experimentInfo)[["abnormUM"]], matrix(rowRanges(experimentInfo)@ranges@pos), binSize)

    ## Output
    outEM <- GRanges(matrix(rowRanges(experimentInfo)@seqnames@values), IRanges(hmmbirdEM$obs[,3], hmmbirdEM$obs[,4]), distance = hmmbirdEM$obs[,2], norm = hmmbirdEM$normAbnorm[,1], abnorm = hmmbirdEM$normAbnorm[,2], direction = hmmbirdEM$obs[,1])

    return(outEM)
}


