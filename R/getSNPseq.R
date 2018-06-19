#' @title getSNPseq
#' @description Given BSgenome object, get the sequence flanking the SNPs.
#' @param snp.vcf A table contain the SNP information. The format is the first five columns of the VCF format plus a column for strand. (CHROM,POS,ID,REF,ALT,strand).
#' @param genome The BSgenome object
#' @param n_flank The number of nucleotide to get around the SNP. 
#' @param IUPAC logcal option. Whether to output IUPAC code for SNP instead of outputing two sequence.
#' @export
getSNPseq <- function(snp.vcf, genome, n_flank , IUPAC = FALSE){
  
  nRef <- sapply(snp.vcf$REF,function(x){nchar(as.character(x))}) # count number of nt for Ref allele
  frontSeq <- getSeq( genome,snp.vcf$CHROM ,start =( snp.vcf$POS -n_flank),end = ( snp.vcf$POS-1 ), strand = snp.vcf$strand,as.character =T )
  tailSeq <- getSeq( genome,snp.vcf$CHROM ,start =( snp.vcf$POS +nRef),end = ( snp.vcf$POS + nRef + n_flank-1 ), strand = snp.vcf$strand,as.character =T )
  
  ## Get Ref and Alt allele for correspond strand
  REFs <- apply(snp.vcf,1, function(x){
    if(x["strand"] == "-"){
      return(revComp(as.character(x["REF"])))
    }else{
      return(as.character(x["REF"]))
    }
  })
  
  ALTs <- apply(snp.vcf,1, function(x){
    if(x["strand"] == "-"){
      return(revComp(as.character(x["ALT"])))
    }else{
      return(as.character(x["ALT"]))
    }
  })
  
  ## handle SNP as IUPAC ambigueous code
  if(IUPAC){
    ## filter out Indels, keep only SNPs. 
    keep.id <- which( nchar(REFs) == 1 & nchar(ALTs) == 1 )
    snp_IUPAC <- .SNPtoIUPAC(REFs[keep.id],ALTs[keep.id])
    
    seqComponent <- cbind("strand"=as.character(snp.vcf$strand[keep.id]),"frontSeq"=frontSeq[keep.id],"tailSeq"=tailSeq[keep.id],snp_IUPAC)
    rownames(seqComponent) <- snp.vcf[keep.id,"ID"]
    finalSeq <- apply(seqComponent,1,function(x){
      if(x["strand"]=="+"){
        tmp_seq <- paste0(x["frontSeq"],x["snp_IUPAC"],x["tailSeq"])
      }else{
        tmp_seq <- paste0(x["tailSeq"],x["snp_IUPAC"],x["frontSeq"])
      }
      return(tmp_seq)
    }) 

  }else{
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
  }
  
  return(finalSeq)
}

.SNPtoIUPAC <- function(x,y){
  
  IUPAC_code <- c("W","S","M","K","R","Y","W","S","M","K","R","Y")
  names(IUPAC_code) <- c("AT","CG","AC","GT","AG","CT","TA","GC","CA","TG","GA","TC")
  return( IUPAC_code[paste0(x,y)] )
}

