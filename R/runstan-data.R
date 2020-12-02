

#' \code{rundata} Run stan for RWD
#'
#' @title rundata
#'
#' @param dat Data for stan model
#' @param stancode Type of stan model to run
#' @param keepall Whether to keep all results
#' @param seed Seed for stan.  If inform, seed must be length 2
#' @param iter Number of iterations for MCMC
#' @param chains Number of chains
#' @param keepall Whether to keep all results or only summary
#' @export
rundata <- function(data, stancode, seed, iter = 1000,
                    chains = 1,  keepall = F, ...) {

  # run stan model
  fit <- rstan::stan(file= stancode, data= data, seed = seed, iter = iter, chains = chains, ...)

  # get output
  if(keepall) {
    out <- list(data = data, fit = fit, stantype = stancode)
  } else {
    out <- list(data = data, summary = summary(fit)$summary, stantype = stancode)
  }

  return(out)
}


#' \code{tblchain} Format stan output for 1 chain as tibble
#'
#' @title tblchain
#'
#' @param x data from one chain of stan, i.e. stanfitATsim$samples list #1
#'@param by1 Skips for posterior draws
#' @param start Where to start by
#' @export
tblchain <- function(x, by1 = 1, start = 5000) {
  x <- tibble::as_tibble(x)
  #print(dim(x))
  samps <- seq(start, nrow(x), by = by1)
  x <- x %>% dplyr::mutate(., iters = seq(1, nrow(x))) %>%
    dplyr::filter(., iters %in% samps) %>%
    tidyr::pivot_longer(., names_to = "var", values_to = "value", -iters)
  x
}

#' \code{unlstan} Unlist stan samples
#'
#' @title unlstan
#'
#' @param stanfit results from call to stan
#' @param by1 Skips for posterior draws
#' @param start Where to start by
#' @export
unlstan <- function(stanfit, by1 = 1, start = 5000) {
  samples <- lapply(stanfit@sim$samples, function(x) tblchain(x, by1 = by1, start = start))
  print(dim(samples[[1]]))
  names(samples) <- paste0("chain", length(samples))
  result <- tibble::tibble(samples) %>%
    tibble::rowid_to_column(.,"chain") %>%
    tidyr::unnest(samples) %>%
    dplyr::mutate(., chain = factor(chain))
  result
}

