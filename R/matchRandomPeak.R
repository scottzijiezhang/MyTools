#' @title generateRandomPeak
#' @param gtf GTF files
#' @param peakWidth mean width of the random peak
#' @param seed Seed for random number
#' @param intervalPerPeak average interval on transcript between each sampled peak 
#' @export
generateRandomPeak <- function(gtf, peakWidth = 70, seed = 111, intervalPerPeak = 1e3){
  set.seed(seed)
  geneModel <- gtfToGeneModel(gtf )
  
  transcribed <- range( geneModel )
  
  randomPeaks <- lapply(transcribed, function(x){
    tx_length <- width(ranges(x))
    n_peak <- max(1,round(tx_length/intervalPerPeak))
    size_peak <-  sapply( rnorm(n_peak,peakWidth,sd = 25) ,max, 10 )
    onTx_loc <- sample((peakWidth/2):(tx_length-peakWidth/2),n_peak)
    return(GRanges(seqnames = seqnames(x), IRanges( start(x)+round(onTx_loc - size_peak/2), start(x)+round(onTx_loc + size_peak/2) ),strand = strand(x)) )
  })
  randomPeaks <- unlist(GRangesList(randomPeaks ) )
  names(randomPeaks) <- 1:length(randomPeaks) 
  
  randomPeaks_anno <- as.GRanges( ChIPseeker::annotatePeak(peak = randomPeaks,TxDb = txdb, genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron","Promoter", "Downstream", "Intergenic") ) )
  
  
  return(randomPeaks_anno)
}


#' @title matchRandomPeak
#' @param peak.gr The GRange object of peak
#' @param randomPeaks The GRange object of random Peaks
#' @param txdb The TxDb object for the annotation
#' @export
matchRandomPeak <- function(peak.gr, randomPeaks ,txdb, nPeaks = 1e10 ){
  
  
  ## Remove target peaks from random peaks
  randomPeaks_filter <- randomPeaks[-c( subjectHits(findOverlaps(peak.gr,randomPeaks) ) ) ]
  
  randomPeakAnno <- as.data.frame( randomPeaks_filter )
  cat("Annotating input peaks\n")
  peakAnno <- ChIPseeker::annotatePeak(peak = peak.gr,TxDb = txdb, genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron","Promoter", "Downstream", "Intergenic") ,sameStrand = T,verbose = F)
  cat("Genomic annotation of the input peak: \n")
  print(peakAnno@annoStat )
    
  cat("Matching genomic annotion for random peaks\n...")
  
  ## 5'UTR
  fiveUTRfreq <-peakAnno@annoStat[,"Frequency"][grep("5' UTR",peakAnno@annoStat[,"Feature"] )]/100
  n5UTR <- length( grep("5' UTR",randomPeakAnno$annotation) )
  
  ## 3'UTR
  threeUTRfreq <-peakAnno@annoStat[,"Frequency"][grep("3' UTR",peakAnno@annoStat[,"Feature"] )]/100
  n3UTR <- length( grep("3' UTR",randomPeakAnno$annotation) )
  cat("...")
  ## 1st Exon
  firstExonfreq <-peakAnno@annoStat[,"Frequency"][grep("1st Exon",peakAnno@annoStat[,"Feature"] )]/100
  nfirstExon <- length( grep("exon 1 of",randomPeakAnno$annotation) )
  
  ## Other Exon
  otherExonfreq <-peakAnno@annoStat[,"Frequency"][grep("Other Exon",peakAnno@annoStat[,"Feature"] )]/100
  nOtherExon <- length( grep("exon 1 of",randomPeakAnno$annotation[grep("Exon",randomPeakAnno$annotation)],invert = T ) )
  cat("...\n")
  ## 1st Intron
  firstIntronfreq <-peakAnno@annoStat[,"Frequency"][grep("1st Intron",peakAnno@annoStat[,"Feature"] )]/100
  nfirstIntron <- length( grep("intron 1 of",randomPeakAnno$annotation) )
  
  ## Other Intron
  otherIntronfreq <-peakAnno@annoStat[,"Frequency"][grep("Other Intron",peakAnno@annoStat[,"Feature"] )]/100
  nOtherIntron <-  length( grep("intron 1 of",randomPeakAnno$annotation[grep("Intron",randomPeakAnno$annotation)],invert = T ) )
  cat("...\n")
  ## DistalIntergenic
  DistalIntergenicfreq <-peakAnno@annoStat[,"Frequency"][grep("Distal Intergenic",peakAnno@annoStat[,"Feature"] )]/100
  nDistalIntergenic <-  length( grep("Distal",randomPeakAnno$annotation) ) 
  cat("...\n")
  ## Determine rate limiting factor
  nSample <- round( 0.98* min( c(n5UTR/fiveUTRfreq, n3UTR/threeUTRfreq, nfirstExon/firstExonfreq, nfirstIntron/firstIntronfreq, nOtherExon/otherExonfreq, nOtherIntron/otherIntronfreq, nDistalIntergenic/DistalIntergenicfreq) ) )
  nSample <- min(nSample, length(peak.gr),nPeaks)
  
  sample_id <-   c(sample( grep("5' UTR",randomPeakAnno$annotation),size = nSample*fiveUTRfreq ),
                   sample( grep("3' UTR",randomPeakAnno$annotation),size = nSample*threeUTRfreq ),
                   sample( grep("exon 1 of",randomPeakAnno$annotation) ,size = nSample*firstExonfreq ),
                   sample( grep("exon 1 of",randomPeakAnno$annotation[grep("Exon",randomPeakAnno$annotation)],invert = T ), size = nSample*otherExonfreq ),
                   sample( grep("intron 1 of",randomPeakAnno$annotation), size = nSample*firstIntronfreq ),
                   sample( grep("intron 1 of",randomPeakAnno$annotation[grep("Intron",randomPeakAnno$annotation)],invert = T ), size = nSample*otherIntronfreq),
                   sample( grep("Distal",randomPeakAnno$annotation), size = nSample*DistalIntergenicfreq)
                   )
  cat("Done")
  return(randomPeaks_filter[unique(sample_id)] )
}

## set inherited function for class csAnno
setMethod("summary", signature(object="csAnno"),
          function(object) {
            
            if (object@hasGenomicAnnotation) {
              return(object@annoStat)
            }
          }
)

