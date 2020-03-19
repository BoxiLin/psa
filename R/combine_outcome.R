#' Combine GWAS outcome of each batch
#'
#' Combine GWAS outcome of each batch
#'
#' @param path path until the model name
combine_outcome <- function(path){
  res <- lapply(path, function(p) {
    batch_res <- readRDS(p)
    info = batch_res$info
    gwas = batch_res$out
    return(list(info, gwas))
  })
  info <- do.call(rbind, lapply(1:23, function(b){
    res[[b]][[1]]
  }))
  gwas <- do.call(rbind, lapply(1:23, function(b){
    res[[b]][[2]]
  }))
  return(list(info, gwas))
}
