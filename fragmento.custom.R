#####################################################
# This is a custom script that looks
# at fragment sizes in the hybrid capture
# data per window (specified in GRange)
# and calculates short/long fragment ratio
# Short = <150bp
# Long = >150bp
#
# Written by Sami Ul Haq, Apr 7 2023
#


library(Rsamtools)

buncha.bams <- c( Sys.glob("../trimmed_sequences/*BASELINE*mkd.bam") , Sys.glob("../trimmed_sequences/*RELAPSE*.mkd.bam") )

load("/cluster/home/sulhaq/hybrid_capture_regions/hybrid.capture.regions.GRange.RData")

for(eachBam in buncha.bams) {
  bamFile <- BamFile(eachBam)

  # hybridcap <- read.table("regions.coding_and_tss.bed")
  # colnames(hybridcap) <- c("chr", "start", "end")

  #where <- makeGRangesFromDataFrame(hybridcap)
  #save(where, file="hybrid.capture.regions.GRange.RData")


  what <- c("rname", "pos", "isize")
  param <- ScanBamParam(what=what, which=where)

  bamo <- scanBam(bamFile, param=param)

  fragout <- c()
  for(blah in names(bamo)) {
    fragsize <- abs(bamo[[blah]]$isize)
    short.frags <- fragsize[fragsize < 150]
    long.frags <- fragsize[fragsize > 150]
    sl.ratio <- median(short.frags)/median(long.frags)
    fragout <- c(fragout, sl.ratio)
  }
  names(fragout) <- names(bamo)
  save(fragout, file=paste0("../sami_fragmento_analyzo/", eachBam, "_fragments.RData"))
}
