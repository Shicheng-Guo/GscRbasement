library(rtracklayer)
library(GenomicRanges)
peak_files <- list.files(pattern = "*.narrowPeak$")
peak_files
peak_granges <- lapply(peak_files, import)
peak_granges
peak_grangeslist <- GRangesList(peak_granges)
peak_grangeslist
peak_coverage <- coverage(peak_grangeslist)
peak_coverage
covered_ranges <- slice(peak_coverage, lower=15, rangesOnly=T)
covered_ranges
covered_granges <- GRanges(covered_ranges)
covered_granges
export(covered_granges, "./consensuspeak/consensus.bed")
reduced_covered_granges<-reduce(covered_granges, min.gapwidth=150)
export(reduced_covered_granges, "./consensuspeak/reduced.consensus.bed")
