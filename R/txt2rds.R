#' txt file to RDS file
#'
#' @param i chromosome indicator, one of 1, ..., 23, "X"; for print purpose only
#' @param input input txt file
#' @param output output Rds file
txt2rds <- function(i, input, output) {
  message(paste("Cleaning Chromosome",i,"..."))

  ### load genetype data
  message("load genetype data...")
  raw_snp = fread2(file = input)

  p = dim(raw_snp)[2]
  N = dim(raw_snp)[1]

  ### Modify patient id to align with phenotype data
  message("remove SNP with ALT containiing `,` ...")
  bad =  grep(",", raw_snp$ALT)
  message(paste("On chromosome ", i, ": ", length(bad)," SNPs were removed...", sep = "" ))
  genotype_data = dplyr::filter(raw_snp, !(1:N %in%  bad))
  colnames(genotype_data) = gsub('\\.', "-", colnames(genotype_data))
  genotype_data[,6:p] <- sapply(genotype_data[,6:p], as.numeric)

  saveRDS(genotype_data, file = output)

  message("Done!")
}
