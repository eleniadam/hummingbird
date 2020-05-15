## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----sampleDataset------------------------------------------------------------
library(hummingbird)
data(exampleHummingbird)

## ----sampleDataset_pos--------------------------------------------------------
m_pos[1:6,1]

## ----sampleDataset_norm-------------------------------------------------------
m_normM[1:6,1:4]
m_normUM[1:6,1:4]

## ----sampleDataset_abnorm-----------------------------------------------------
m_abnormM[1:6,1:4]
m_abnormUM[1:6,1:4]

## ----sampleDataset_hmmbird_em-------------------------------------------------
str(hmmbird1)

## ----sampleDataset_hmmbird_postadj--------------------------------------------
str(hmmbird2)

## ----hummingbird--------------------------------------------------------------
library(hummingbird)
data(exampleHummingbird)

## ----hummingbird_em-----------------------------------------------------------
hmmbird1 <- hummingbirdEM(normM=m_normM, normUM=m_normUM, abnormM=m_abnormM, 
abnormUM=m_abnormUM, pos=m_pos, binSize=40)

## ----hummingbird_postAdjustment-----------------------------------------------
hmmbird2 <- hummingbirdPostAdjustment(em=hmmbird1$obs, pos=m_pos, minCpGs=10, 
minLength=100, maxGap=300)
hmmbird2$DMRs

## ----hummingbird_graph--------------------------------------------------------
hummingbirdGraph(pos=m_pos, normM=m_normM, normUM=m_normUM, abnormM=m_abnormM, 
abnormUM=m_abnormUM, dmrs=hmmbird2$DMRs, coord1=107991, coord2=108350)

