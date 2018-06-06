#' @title getSNPseq
#' @description Given BSgenome object, get the sequence flanking the SNPs.
#' @param snp.vcf A table contain the SNP information. The format is the first five columns of the VCF format plus a column for strand. (CHROM,POS,ID,REF,ALT,strand).
#' @param genome The BSgenome object
#' @param n_flank The number of nucleotide to get around the SNP. 
#' @export
getSNPseq <- function(snp.vcf, genome, n_flank ){
  
  nRef <- sapply(snp.vcf$REF,function(x){nchar(as.character(x))}) # count number of nt for Ref allele
  frontSeq <- getSeq( genome,snp.vcf$CHROM ,start =( snp.vcf$POS -n_flank),end = ( snp.vcf$POS-1 ), strand = snp.vcf$strand,as.character =T )
  tailSeq <- getSeq( genome,snp.vcf$CHROM ,start =( snp.vcf$POS +nRef),end = ( snp.vcf$POS + nRef + n_flank-1 ), strand = snp.vcf$strand,as.character =T )
  
  ## Get Ref and Alt allele for correspond strand
  REFs <- apply(snp.vcf,1, function(x){
    if(x$strand == "-"){
      return(rev.comp(as.character(x$REF)))
    }else{
      return(as.character(x$REF))
    }
  })
  
  ALTs <- apply(snp.vcf,1, function(x){
    if(x$strand == "-"){
      return(rev.comp(as.character(x$ALT)))
    }else{
      return(as.character(x$ALT))
    }
  })
  
  seqComponent <- cbind("strand"=as.character(snp.vcf$strand),frontSeq,tailSeq,REFs,ALTs)
  
  finalSeq <- t( apply(seqComponent,1,function(x){
    if(x["strand"]=="+"){
      REFseq <- paste0(x["frontSeq"],x["REFs"],x["tailSeq"])
      ALTseq <- paste0(x["frontSeq"],x["ALTs"],x["tailSeq"])
    }else{
      REFseq <- paste0(x["tailSeq"],x["REFs"],x["frontSeq"])
      ALTseq <- paste0(x["tailSeq"],x["ALTs"],x["frontSeq"])
    }
    return(cbind(REFseq,ALTseq))
  }) )
  
  colnames(finalSeq) <- c("RefSequence","AltSequence")
  rownames(finalSeq) <- snp.vcf$ID
  return(finalSeq)
}