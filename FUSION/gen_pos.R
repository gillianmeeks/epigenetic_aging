#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(stringr)

main <- function(args) {
  parser <- ArgumentParser(description = "")
  parser$add_argument("bedfile")
  parser$add_argument("dir")
  parser$add_argument("pop")
  parser$add_argument("-o", "--output", type = "character", default = stdout())
  args <- parser$parse_args(args)

  bed <- fread(args$bedfile, sep = "\t")
  N <- ncol(bed) - 4

  cat(paste(c("DIR", "WGT", "ID", "CHR", "P0", "P1", "N"), collapse = "\t"), "\n", file = args$output)

  for (i in 1:nrow(bed)) {
    chr <- bed[[1, 1]]
    name <- bed[[1, "gene"]]
    p0 <- bed[[1, "start"]]
    p1 <- bed[[1, "end"]]
    wgt <- paste0(name, ".wgt.RDat")
    fname <- file.path(args$dir, paste0(name, ".wgt.RDat"))

    if (!file.exists(fname) || file.size(fname) == 0) {
      next
    }

    items <- c(args$dir, wgt, name, chr, p0, p1, N)
    cat(paste(items, collapse = "\t"), "\n", file = args$output)
  }

  return(0)
}

if (basename(commandArgs(trailingOnly = FALSE)[1]) == "Rscript") {
  sys.exit(main(commandArgs(trailingOnly = TRUE)))
}
