#' Gwas for single chromosome;
#'
#' @param chr indictor for which chromosome
#' @param model indictor for which model
#' @param phenotype vector of response
#' @param covariate list of names covariates to include in regression model
#' @param batch_size
#' @param ncores
#' @param data_path
#' @param save_path
single_gwas <- function(chr, model, phenotype, covariate="sex",
                        batch_size = 10, ncores = 1,
                        data_path, save_path) {

  if(!model %in% c("baseline", "log")) {message("  Error: invalid model!"); stop()}

  message(paste("Chromesome",chr,Sys.time()))

  message("  ### load genetype data...")
  genotype_data <- readRDS(data_path)
  id <- intersect(colnames(genotype_data), row.names(phenotype))
  info <- select(genotype_data, one_of(c("CHROM","POS","ID")))
  genotype_data <- genotype_data %>% select(id)
  n.snps = dim(genotype_data)[1]

  message("  ### load phenotype data...")
  p = pull(phenotype[id, ],  pa_age)
  if (model=="log") {p = log(p)}
  cov = as.matrix(select(phenotype[id,],covariate))

  message("  ### batch-wise gwas")
  batch = split(1:n.snps, cut(seq_along(1:n.snps), batch_size, labels = FALSE))
  gwas.list <- lapply(1:batch_size, function(j) {
    print(paste("Chromesome",chr,"; batch", j, model, Sys.time()))
    g = as_FBM(t(slice(genotype_data,batch[[j]])))
    big_univLinReg(g, y.train = p, covar.train = cov, ncores = ncores)
  })
  gwas<- list(info = info, out = do.call(rbind, gwas.list))
  saveRDS(gwas, file = save_path)
}
