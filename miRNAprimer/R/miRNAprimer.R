#' @title Derive primers for miRNA stem-loop qPCR.
#'
#' @description This package provides a configurable interface for generating
#'   specific primers for miRNA stem-loop quantitative reverse transcription PCR (RT-qPCR).
#'   It follows both the design principles of PMID: 21732315 and the MIQE of qPCR doi: 10.1373/clinchem.2008.112797.
#'   It can use any stem-loop sequence, and you can modify the extension sequence for the 5' end of the forward primer.
#'   Melting temperature is automatically calculated using the TmCalculator package.
#'
#' @param x, s, r
#'
#' @return Res
#'
#' @examples
#'
#' @export
#'
miRNAprimer <- function(x, s, r) {
  # Derives qPCR primers for miRNA expression analysis 5' to 3'.
  #
  # Args:
  #   x: Mature miRNA sequence/s whose primers is to be derived.
  #   s: Stem-loop sequence to be used for the cDNA synthesis primer.
  #   r: Universal reverse primer associated with the stem-loop.
  #
  # Returns:
  #   The primers for the target miRNA, including cDNA synthesis primer,
  #    universal reverse primer and miRNA specific foward primers with Tm.
  n <- nchar(x)
  # Error handling
  if (n < 15 || n > 25) {
    stop("Arguments x must be mature miRNA sequence: see http://www.mirbase.org.")
  }
  #
  UtoT <- function(x) {
    # Converts U(Uracil:RNA) to T(thymine:DNA).
    #
    # Args:
    #   x: Nucleic acid sequence/s to be converted.
    #
    # Returns:
    #  The convereted nucleic acid sequence/s of same length as input.
    chartr(old = "U", new = "T", x)
  }
  #
  substrS <- function(x, n, side) {
    # Extract subset of original sequence/s from  either ends.
    #
    # Args:
    #   x: Sequence/s whose side is to be extracted.
    #   n: The number of nucleotides from the right.
    #   side: Which end c("F"=left, "L"=right) of the sequence/s to extract.
    #
    # Returns:
    #   Specefied subset of the input sequence.
    nl <- nchar(x)
    # Option handling
    if (side == "F") {
      substring(x, 1, nl-n)
    } else if (side == "L") {
      substring(x, nl-n+1, nl)
    } else {
      stop("Njehie: Side must be F(left) or L(right).")
    }
  }
  #
  strRev <- function(x) {
    # Reverses string left to right.
    #
    # Args:
    #   x: String whose direction is to be reversed.
    #
    # Returns:
    #   The reversed version of the input string.
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  }
  #
  #Convert Stem-loop to complement
  ToComp <- function(x) {
    # Generate the complement of a nucleic acid Sequence/s.
    #
    # Args:
    #   x: Sequence/s whose complement is to be generated.
    #
    # Returns:
    #   The complement of the input sequence.
    bases <- c("A","C","G","T")
    # Handling sequence at nucleotide level
    x <- unlist(strsplit(toupper(x), NULL))
    paste(unlist(lapply(x, function(x) {
      if (x == "A") compSt <- "T"
      if (x == "C") compSt <- "G"
      if (x == "G") compSt <- "C"
      if (x == "T") compSt <- "A"
      if (!x %in% bases) compSt <- "N"
      return(compSt)
    })), collapse="")
  }
  # Note:
  # Lines are in parenthesis to force print.
  (miRNA <- UtoT(x)) #Convert U to T.
  # Get 6nt for stem-loop primer.
  StemlpUnique <- ToComp(strRev(substrS(miRNA, 6, "L")))
  # Get sequence for forward primer design.
  pFWR <- substrS(miRNA, 6, "F")
  # Sliding window
  Lnt <- seq(13, nchar(pFWR), by=1)
  #
  pFWR1 = NULL
  # Generate Raw windows
  for (i in 1:length(Lnt)){
    pFWR1[i] = substring(pFWR, 1, Lnt[i])
  }
  # Stablise melting temperatures, with "GCTGGCA".
  Ext5end <- function(x,y) {
    # Stablise melting temperatures, by extending the 5' end.
    #  Uses different lengths of "GCTGGCA".
    #
    # Args:
    #   x: Primer sequence/s whose 5' end is to be extended.
    #   y: Sequence to be added to the 5' ends.
    #
    # Returns:
    #   Extended sequences.
  # User can change the extender.
  Exts = NULL
  for (i in 1:length(pFWR1)){
    Exts[i] = paste(y, x[i], sep="")
  }
  return(Exts)
  }
  #
  FWRs <- Ext5end(pFWR1,"GCTGGCA")
  FWRsa <- Ext5end(pFWR1,"TGGCA")
  FWRsb <- Ext5end(pFWR1,"GGCA")
  #
TmNN <- function(x) {
  # Calculate primer melting tempurates using SantaLucia J (2004).
  #
  # Args:
  #   x: Primer sequence/s whose melting temperature/s is to be calculated.
  #
  # Returns:
  #   Temperatures in a vector of equal length to input.
  Tms <- NULL
  for (i in 1:length(x)){
    Tms[i] = TmCalculator::Tm_NN(x[i], nn_table = "DNA_NN4")
  }
  return(Tms)
  }
  #
  pFWR1Tm <- TmNN(pFWR1)
  FWRsTm <- TmNN(FWRs)
  FWRsTma <- TmNN(FWRsa)
  FWRsTmb <- TmNN(FWRsb)
  #
  RevTm <- Tm_NN(r, nn_table = "DNA_NN4")
  # Bind pFWR1 primers and their Tms.
  FWR1 <- cbind(Fowards=c(pFWR1), Tm=c(round(pFWR1Tm)))
  # Bind FWR primers and their Tms.
  FWR <- cbind(FWRs, round(FWRsTm))
  FWRa <- cbind(FWRsa, round(FWRsTma))
  FWRb <- cbind(FWRsb, round(FWRsTmb))
  # Outputs
  # Bind REV primer and its Tm.
  Rev <- c(paste("Reverse", r, sep = ":"),
           paste(round(RevTm), sep = ":"))
  # Bind stemloop with miRNA 6nt unique seq
  Stem <- c("Stem-loop:", paste(s,StemlpUnique, sep="|"))
  #All
  Res <- rbind(FWR1,FWR,FWRa,FWRb, Rev, Stem)
 return(Res)
}
