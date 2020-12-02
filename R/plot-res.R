
#' \code{plotstan} Create plots for assessing stan diagnostic
#'
#' @title plotstan
#' @param typesim Which simulation to run.  Currently, local1 - local4, ambient.
#' @param stanres Output from runstan
#' @param dirname Folder name where to save plots
#' @param filename Base filename or typesim if NULL
#' @param by Skips for traceplot, biasplot
#' @param prof Profile data, defaults to prof.
#' @param meansd Matrix with columns corresponding to source, type, mean, sd, defaults to meansd.
#' @param pdf Whether to save pdfs of energy, rhat, bias
#' @param hten Height in inches of energy pdf
#' @param wden Width in inches of energy pdf
#' @param htrh Height in inches of rhat pdf
#' @param htbi Height in inches of bias pdf
#' @param wden Width in inches of bias pdf
#' @export
plotstan <- function(typesim, stanres, dirname = NULL,
                     # only need defaults
                     filename = NULL, by = 5,
                     prof = prof, meansd = meansd, pdf = F,
                     # change heights for outputs
                     hten = 400, wden = 20, htrh = 100, htbi = 10, wdbi = 10) {

  # Get filename if not provisted
  if(is.null(filename)) {
    filename <- typesim
  }

  # Get necessary input
  sources <- dat$dat$f$source
  mat1 <- mat1fun(stanres, sources)


  # Get simulation-specific mean/SD, profiles
  prof <- dplyr::filter(prof, type == typesim)
  meansd <- dplyr::filter(meansd, type == typesim)


  # Plots
  tr <- mytraceplot(stanres, dirname, filename, by = by1)
  pa <- pairsplot(stanres, dirname, filename, mat1, sources)
  en <- energyplot(stanres, dirname, filename, pdf, ht = hten, wd = wden)
  rh <- rhatplot(stanres, dirname, filename, pdf, ht = htrh)
  bi <- biasplot(stanres, dirname, filename,
        mat1, prof, meansd, typesim, by = by1, pdf, ht = htbi, wd = wdbi)


  # output (save energy, rhat, bias)
  list(free = mat1, sources = sources,
       en = en, rh = rh, bi = bi)
}




#' \code{mat1fun} Get constraints with labels
#'
#' @title mat1fun
#'
#' @param stanres Output from runstan
#' @param sources Names of sources
#' @export
mat1fun <- function(stanres, sources) {
  # get constraint names for nVF
  P <- dat$standat$P
  L <- dat$standat$L

  zeromat <- dat$standat$zeromat
  mat1 <- matrix(paste0(rep(sources, P), "-",
                        rep(colnames(zeromat), each = L)),
                 byrow = F, nrow = L)
  mat1 <- as.vector(mat1)[(as.vector(zeromat)) != 0]

  mat1
}





#' \code{mytraceplot} Create traceplots for stan output
#'
#' @title mytraceplot
#'
#' @param stanres Output from runstan
#' @param dirname Folder name where to save plots
#' @param filename Base filename
#' @param by Skips for traceplot
#' @export
mytraceplot <- function(stanres, dirname, filename, by = 5) {
  start1 <- dim(stanres$fit)[1]
  res1 <- unlstan(stanres$fit, by1 = by, start = start1)

  fn1 <-  here::here(file.path(dirname,
                         paste0(filename, "-traceplot.pdf")))

  cols <- RColorBrewer::brewer.pal(4, "Dark2")
  pdf(fn1)
  vars <- unique(res$var)
  chains <- unique(res$chain)
  iters1 <- unique(res$iters)
  iters1 <- seq(ceiling(max(iters1) / 2), max(iters1), by = by)
  res <- dplyr::filter(res, iters %in% iters1)

  par(mfrow = c(2, 2))
  for(i in 1 : length(vars)) {
    res0 <- dplyr::filter(res, var == vars[i])
    min1 <- min(res0$value, na.rm = T)
    max1 <- max(res0$value, na.rm = T)
    for(j in 1 : length(chains)) {
      res1 <- dplyr::filter(res0, chain == chains[j]) # %>% arrange(iters)
      if(j == 1) {
        plot(res1$iters, res1$value, col = cols[j], type = "l",
             xlab = "Iteration", ylab = vars[i], main = vars[i],
             ylim = c(min1, max1))
      } else {
        lines(res1$iters, res1$value, col = cols[j])

      }
    }
  }
  dev.off()


}









#' \code{pairsplot} Create pairs plots for stan output
#'
#' @title pairsplot
#'
#' @param stanres Output from runstan
#' @param dirname Folder name where to save plots
#' @param filename Base filename
#' @param mat1 Names of free elements of F
#' @param sources Names of sources
#' @export
pairsplot <- function(stanres, dirname, filename, mat1, sources) {

  filename1 <- paste0(filename, "-pairs.pdf")
  pdf(here::here(dirname, filename1))

  # mu G/sigma G
  lab1 <- paste0(rep(c("mug", "sigmag"), each = length(sources)), rep(sources, 2))
  pairs(stanres$fit, pars = c("mug", "sigmag"), labels = lab1, condition = "energy")

  # All vf
  types <- c("mug", "sigmag")
  nF <- names(stanres$fit)[grep("nvF", names(stanres$fit))] %>% length()
  nF1 <- ceiling(nF / 8)
  for(j in 1 : 2) {
    k <- 1
    first5 <- paste0(rep(types[j], length(sources)), sources)
    print(c("first5", first5))
    for(i in 1 : nF1) {
      l <- min(c((k + 7), nF))
      print(c(k, l))
      lab1 <- c(first5, mat1[k : l])
      print(lab1)
      pairs(stanres$fit, labels = lab1,
            pars = c(types[j], paste0("vF[", (k : l), "]")), condition = "energy")

      k <- k + 8
    }
  }

  dev.off()

}






#' \code{energyplot} Create energy plots for stan output
#'
#' @title energyplot
#'
#' @param stanres Output from runstan
#' @param dirname Folder name where to save plots
#' @param filename Base filename
#' @param pdf Whether to save pdf, defaults to F
#' @param ht Height of pdf
#' @param wd Width of pdf
#' @export
energyplot <- function(stanres, dirname, filename, pdf = F,
                       ht = 400, wd = 20) {

  # Which items to save
  start1 <- dim(stanres$fit)[1]
  nmax <- start1 - start1 / 2
  cn <- dim(stanres$fit)[2]
  energy <- rstan::get_sampler_params(stanres$fit)
  names(energy) <- paste0("chain", seq(1, cn))

  # Get energy for later iterations
  energy1 <- lapply(energy, function(x) {
    data.frame(x) %>%
      dplyr::select(., energy__) %>%
      dplyr::mutate(., iters = seq(1, nrow(x))) %>%
      dplyr::filter(., iters > nmax )
  }) %>%
    tibble::as_tibble_col(.) %>%
    tibble::rowid_to_column(., "chain") %>%
    tidyr::unnest(., value)

  # fix names
  params <- rstan::extract(stanres$fit, permuted = F) %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(., "iters") %>%
    tidyr::pivot_longer(., names_to = "var", values_to = "val", -iters)

  # Get chain info (for subsetting)
  params1 <- params %>%
    dplyr::mutate(., chain = as.numeric(substr(var, 7, 7)), var = substring(var, 8))

  # Remove some iterations
  energy2 <- dplyr::mutate(energy1, iters = iters - min(iters)) %>%
    dplyr::full_join(., params1)


  # print highest correlations
  dplyr::group_by(energy2, chain, var) %>%
    dplyr::summarize(., cor1 = cor(energy__, val)) %>%
    dplyr::mutate(., abscor = abs(cor1)) %>%
    dplyr::arrange(., desc(abscor)) %>%
    dplyr::slice(1 : 5) %>% print()

  # which to plot
  # samps <- sample(seq(1, nmax), 500, replace = F)
  samps <- seq(nmax, start1, length = 100) %>%
    round() %>% unique()

  # keep chain = 1
  energy3 <- dplyr::filter(energy2, iters %in% samps, chain == 1) %>%
    dplyr::mutate(., chain = factor(chain))
  # could add lables for var here

  # plot
  g1 <- ggplot2::ggplot(energy3, aes(x = val, y = energy__, colour = chain)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~var, scales = "free_x", ncol = 8)

  if(pdf) {
    filename1 <- paste0(filename, "-energy.pdf")

    ggplot2::ggsave(here(file.path(dirname, filename1)),
           g1, height = ht, width = wd, limitsize = F)
  }

  g1
}







#' \code{rhatplot} Create rhat plots for stan output
#'
#' @title rhatplot
#'
#' @param stanres Output from runstan
#' @param dirname Folder name where to save plots
#' @param filename Base filename
#' @param pdf Whether to save pdf, defaults to F
#' @param ht Height of pdf
#' @export
rhatplot <- function(stanres, dirname, filename, pdf = F, ht = 100) {

  # Get Rhat and plot
  rhats <- bayesplot::rhat(stanres$fit)
  g1 <-  bayesplot::mcmc_rhat(rhats) + bayesplot::yaxis_text(hjust = 1)

  # If save
  if(pdf) {
    filename1 <- paste0(filename, "-rhat.pdf")

    ggplot2::ggsave(here(file.path(dirname, filename1)),
           g1, height = ht, limitsize = F)
  }

  g1

}




#' \code{biasplot} Create energy plots for stan output
#'
#' @title energyplot
#'
#' @param stanres Output from runstan
#' @param dirname Folder name where to save plots
#' @param filename Base filename
#' @param by Skips for parameter posterior draws (boxplot)
#' @param pdf Whether to save pdf, defaults to F
#' @param ht Height of pdf
#' @param wd Width of pdf
#' @export
biasplot <- function(stanres, dirname, filename,
                     mat1,
  prof = prof, meansd = meansd, typesim = typesim,
  by = 5,
  pdf = F, ht = 10, wd = 10) {

  # Get data for boxplot of posterior
  iters <- dim(stanres$fit)[1]
  sel <- seq(1, iters, by = by)
  params <- rstan::extract(stanres$fit, permuted = F) %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(., "iters") %>%
    tidyr::pivot_longer(., names_to = "var", values_to = "val", -iters) %>%
    dplyr::mutate(., chain = as.numeric(substr(var, 7, 7)), var = substring(var, 8)) %>%
    dplyr::filter(chain == 1, iters %in% sel) %>%
    dplyr::select(-c(iters, chain)) %>%
    dplyr::separate(var, c("var1", "rowcol"), "\\[") %>%
    dplyr::separate(rowcol, c("row", "col"), ",") %>%
    dplyr::mutate(var1 = gsub("\\.", "", var1),
           col = gsub("\\]", "", col),
           row = gsub("\\]", "", row),
           row = as.numeric(row),
           col = as.numeric(col)) %>%
    dplyr::filter(!(var1 %in% c("nvF", "lp__")))


  # truth
  truth <- stanres$dat
  g <- dplyr::rename(truth$g, row = id) %>%
    tidyr::pivot_longer(-row, values_to = "truth")  %>%
    dplyr::mutate(col = as.numeric(factor(name)), var1 = "G")
  f <- dplyr::filter(prof, type == typesim, is.na(constraint)) %>%
    dplyr::select(poll, source, scale1) %>%
    dplyr::rename(truth = scale1) %>%
    tidyr::unite(name, c(source, poll), sep = "-") %>%
    dplyr::mutate(row = factor(name, levels = mat1, labels = seq(1, length(mat1))),
           row = as.numeric(row),
           var1 = "vF")
  musigg <- dplyr::filter(meansd, type == typesim) %>%
    dplyr::select(-type) %>%
    dplyr::mutate(row = as.numeric(factor(source))) %>%
    tidyr::pivot_longer(-c(source, row),
                 names_to = "var1", values_to = "truth") %>%
    dplyr::mutate(var1 = factor(var1, levels = c("mean", "sd"),
                         labels = c("mug", "sigmag"))) %>%
    dplyr::rename(name = source)
  P <- unique(profdat$poll) %>% length()
  sigmaeps <- data.frame(row = seq(1, P), var1 = "sigmaeps", truth = truth$sigmaeps)
  truth <- dplyr::full_join(g, f) %>%
    dplyr::full_join(musigg) %>% dplyr::full_join(sigmaeps)

  # Find posterior mean
  means <- dplyr::group_by(params) %>%
    dplyr::group_by(var1, row, col) %>%
    dplyr::summarize(mean = mean(val))

  # combine
  dat <- dplyr::full_join(params, truth) %>%
    dplyr::full_join(means)

  # all but G
  dat1 <- dplyr::filter(dat, var1 != "G")
  g1 <- ggplot2::ggplot(data = dat1) +
    ggplot2::geom_boxplot(aes(x = factor(row), y = val)) +
    ggplot2::geom_point(aes(x = row, y = mean), colour = "red",
               shape = 17) +
    ggplot2::geom_point(aes(x = row, y = truth), colour = "blue",
               shape = 8) +
    ggplot2::facet_wrap(~var1, scales = "free")


  # G
  dat1 <- dplyr::filter(dat, var1 == "G")
  g2 <- ggplot2::ggplot(data = dat1) +
    ggplot2::geom_boxplot(aes(x = factor(row), y = val)) +
    ggplot2::geom_point(aes(x = row, y = mean), colour = "red",
               shape = 17) +
    ggplot2::geom_point(aes(x = row, y = truth), colour = "blue",
               shape = 8) +
    ggplot2::facet_wrap(~col, scales = "free")

  # If save
  if(pdf) {
    filename1 <- paste0(filename, "-biasplot-noG.pdf")

    ggplot2::ggsave(here(file.path(dirname, filename1)),
           g1, height = ht, limitsize = F)

    filename1 <- paste0(filename, "-biasplot-G.pdf")

    ggplot2::ggsave(here(file.path(dirname, filename1)),
           g2, height = ht, limitsize = F)
  }

  list(g1, g2)

}