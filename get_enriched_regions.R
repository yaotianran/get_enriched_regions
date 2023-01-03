#!/usr/bin/env Rscript
# use insert to get enrich region instead of reads/segments
library(Rsamtools)
library(stringr)

get_enriched_regions <- function(bam, chromosomes = 'auto', min_coverage = 1, min_count = 1, max_gap = 3, min_mapq = 1) {
   
   result.list = list()
   
   # get regions
   if (class(bam) == 'character') {
      message('Reading bamfile ...')
      param <- Rsamtools::ScanBamParam(flag = scanBamFlag(isDuplicate = F, isSecondaryAlignment = F, isUnmappedQuery = F), mapqFilter = min_mapq)
      galp <- GenomicAlignments::readGAlignmentPairs(bam, param = param)
      if  (length(chromosomes) == 1 && chromosomes == 'auto') {
         chrom.c = levels(seqnames(galp))
         chrom.c = chrom.c[!str_detect(chrom.c, '_')]
      } else {
         chrom.c = chromosomes
      }
      
      frags <- granges(keepSeqlevels(galp, chrom.c, pruning.mode = "coarse"), on.discordant.seqnames = "drop") # convert aligned segments to cfDNA fragments
      
   } else if (class(bam) == 'GRanges') {
      if  (length(chromosomes) == 1 && chromosomes == 'auto') {
         chrom.c = levels(seqnames(bam))
         chrom.c = chrom.c[!str_detect(chrom.c, '_')]
      } else {
         chrom.c = chromosomes
      }
      frags = bam[seqnames(bam) %in% chrom.c]
   }
   
   strand(frags) <- '*'
   message(sprintf('Read %s fragments', length(frags)))
   result.list$fragments = frags
   
   # pileup
   message('Pileuping ...')
   total.gr = reduce(frags)
   singlebase.gr = unlist(tile(total.gr, width = 1))
   query.c = countOverlaps(singlebase.gr, frags)
   singlebase.gr = singlebase.gr[query.c >= min_coverage]
   
   query.c = countOverlaps(total.gr, singlebase.gr)
   total.gr = total.gr[query.c >= min_count]
   
   # fill gaps
   message(sprintf('Filling gaps no greater than %s bases', max_gap))
   gaps.gr = gaps(total.gr)
   gaps.gr = gaps.gr[width(gaps.gr) <= max_gap]
   results.gr = reduce(c(total.gr, gaps.gr))
   message(sprintf('Enriched %s regions', length(results.gr)))
   
   result.list$regions = results.gr
   return(result.list)
   
}


