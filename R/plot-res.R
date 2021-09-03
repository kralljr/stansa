
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
#' @param typeplot Whether to plot local or ambient
#' @param plot Whether to plot only bias
#' @export
plotstan <- function(typesim, stanres, dirname = NULL,
                     # only need defaults
                     filename = NULL, by = 5,
                     prof = prof, meansd = meansd, pdf = F,
                     # change heights for outputs
                     hten = 250, wden = 20, htrh = 100, htbi = 10, wdbi = 20,
                     typeplot = "none", plot = "bias") {

  # Get filename if not provisted
  if(is.null(filename)) {
    filename <- typesim
  }

  # if joint/sparse
  name1 <- names(stanres$dat)
  name2 <- names(stanres$standat)
  if("ya" %in% name1) {
    if(typeplot == "local") {
      keeps <- name1[grep("l", name1)]
      keeps2 <- name2[grep("l", name2)]
      keeps2 <- keeps2[-which(keeps2 %in% c( "LBl2", "LBl1"))]

    } else {
      keeps <- name1[grep("a", name1)]
      keeps <- keeps[-which(keeps %in% c("sigmaepsl", "namesfl"))]

      keeps2 <- name2[grep("a", name2)]
      keeps2 <- keeps2[-which(keeps2 %in% c("zeromatl", "onematl", "matchamb", "LBa"))]
      # use ambient truth
      typesim <- "ambient"
    }


    temp <- stanres$dat$namesfl[which(stanres$standat$matchamb == 0), ]
    stanres$mx <- paste(temp[, 2], temp[, 1], sep = "-")

    stanres$dat <- stanres$dat[keeps]
    nchar1 <- nchar(names(stanres$dat))
    names(stanres$dat) <- substr(names(stanres$dat), 1, nchar1- 1)


    stanres$standat <- stanres$standat[keeps2]
    nchar1 <- nchar(names(stanres$standat))
    names(stanres$standat) <- substr(names(stanres$standat), 1, nchar1- 1)
  }

  # Get necessary input
  sources <- stanres$dat$f$source
  mat1 <- mat1fun(stanres, sources)
  cons <- stanres$dat$cons
  labels <- makelabels(cons, sources, mat1)


  # Get simulation-specific mean/SD, profiles
  prof <- dplyr::filter(prof, type == typesim)
  meansd <- dplyr::filter(meansd, type == typesim)


  # # Plots

  if(plot == "bias") {
    bi <- biasplot(stanres, dirname, filename,
                   mat1, prof, meansd, typesim, by = by, pdf, ht = htbi, wd = wdbi,
                   labels = labels, typeplot)
    out <- list()
    out$bi <- bi
  } else if(typeplot != "ambient") {
      tr <- mytraceplot(stanres, dirname, filename, by = by, labels = labels)
      en <- energyplot(stanres, dirname, filename, pdf, ht = hten, wd = wden,
                       labels = labels, typeplot = typeplot)
      rh <- rhatplot(stanres, dirname, filename, pdf, ht = htrh)
      out <- list(tr = tr, en = en, rh = rh, labels = labels,
                  free = mat1, sources = sources)
  } else {
    out <- list()
  }

  if(plot != "bias") {

  bi <- biasplot(stanres, dirname, filename,
                 mat1, prof, meansd, typesim, by = by, pdf, ht = htbi, wd = wdbi,
                 labels = labels, typeplot)
  pa <- pairsplot(stanres, dirname, filename, mat1, sources, cons, typeplot)

  out$pa <- pa
  out$bi <- bi
  }

  out

}





#' \code{getdat} Get ambient vs. local data
#'
#' @title getdat
#'
#' @param standat Stan data
#' @param typeplot Ambient or local
#' @param L number of sources
#' @export
#'
getdat <- function(stanres, typeplot, L) {

  name2 <- names(stanres$data)

  if("ya" %in% name2) {
    if(typeplot == "local") {
      keeps2 <- name2[grep("l", name2)]
      keeps2 <- keeps2[-which(keeps2 %in% c( "LBl2", "LBl1"))]


    } else {

      keeps2 <- name2[grep("a", name2)]
      keeps2 <- keeps2[-which(keeps2 %in% c("zeromatl", "onematl", "matchamb", "LBa"))]
      # use ambient truth
      typesim <- "ambient"
    }


    stanres$data <- stanres$data[keeps2]
    nchar1 <- nchar(names(stanres$data))
    names(stanres$data) <- substr(names(stanres$data), 1, nchar1- 1)
  }

  P <- ncol(stanres$data$y)

  zeromat <- standat$zeromat
  mat1 <- matrix(paste0(rep(sources, P), "-",
                        rep(colnames(zeromat), each = L)),
                 byrow = F, nrow = L)
  mat1 <- as.vector(mat1)[(as.vector(zeromat)) != 0]

  list(stanres = stanres, mat1 = mat1)
}





#' \code{mat1fun} Get constraints with labels
#'
#' @title mat1fun
#'
#' @param standat Stan data
#' @param sources Names of sources
#' @export
mat1fun <- function(standat, sources) {
  # get constraint names for nVF

  if("standat" %in% names(standat)) {
    standat <- standat$standat
  }
  P <- standat$P
  L <- standat$L

  zeromat <- standat$zeromat
  mat1 <- matrix(paste0(rep(sources, P), "-",
                        rep(colnames(zeromat), each = L)),
                 byrow = F, nrow = L)
  mat1 <- as.vector(mat1)[(as.vector(zeromat)) != 0]

  mat1
}


#' \code{varn} Extract name (vs. position)
#'
#' @title varn
#'
#' @param var Variable from stan output
varn <- function(var) {
  var <- gsub("\\.", "", var)
  wh <- stringr::str_locate(var, "\\[")[, 1]
  wh <- wh - 1
  wh[is.na(wh)] <- 1
  var <- stringr::str_sub(var, 1, wh)
  var
}


#' \code{rown} Extract row
#'
#' @title rown
#'
#' @param var Variable from stan output
rown <- function(var) {
  stringr::str_extract(var, "\\[\\d*\\]") %>%
    stringr::str_sub(., 2, nchar(.)-1) %>%
    as.numeric()
}


#' \code{rown} Extract row
#'
#' @title rown
#'
#' @param var Variable from stan output
rownG <- function(var) {
  stringr::str_extract(var, "\\[\\d*\\,") %>%
    stringr::str_sub(., 2, nchar(.)-1) %>%
    as.numeric()
}


#' \code{coln} Extract column
#'
#' @title coln
#'
#' @param var Variable from stan output
coln <- function(var) {
  stringr::str_extract(var, "\\,\\d*\\]") %>%
    stringr::str_sub(., 2, nchar(.)-1) %>%
    as.numeric()
}

#' \code{makelabels} Create names for variables
#'
#' @title makelabels
#'
#' @param cons Names of pollutants in order
#' @param sources Names of sources in order
#' @param mat1 Names of free elements in F in order
#' @export
makelabels <- function(cons, sources, mat1) {
   sourcesd <- tibble::as_tibble(sources) %>% tibble::rowid_to_column() %>%
     dplyr::rename(source = value)
   consd <- tibble::as_tibble(cons) %>% tibble::rowid_to_column()%>%
     dplyr::rename(cons = value)
   mat1d <- tibble::as_tibble(mat1) %>% tibble::rowid_to_column()%>%
     dplyr::rename(mat1 = value)
   labels <- dplyr::full_join(sourcesd, consd) %>% dplyr::full_join(., mat1d)

  labels

}

#' \code{getnames} Label function
#'
#' @title getnames
#'
#' @param var Variable name from stan
getnames <- function(var, labels = labels) {

  varname <- varn(var)
  var1 <- data.frame(var = var, varname = varname)
  var1 <- dplyr::mutate(var1, varid = ifelse(varname %in% c("G", "Ga", "Gl"), coln(var), rown(var)),
                 varid2 = ifelse(varname %in% c("G", "Ga", "Gl"), rownG(var), rown(var)))

  sources <- labels$source[var1$varid]
  mat1 <- labels$mat1[var1$varid]
  cons <- labels$cons[var1$varid]

  var1 <- dplyr::mutate(var1, sources = sources, mat1 = mat1, cons = cons,
                 out = dplyr::case_when(varname %in% c("sigmag", "mug", "G",
                                                       "sigmaga", "muga", "Ga",
                                                       "sigmagl", "mugl", "Gl") ~ sources,
                                 varname %in% c("nvF", "vF",
                                                "nvFa", "vFa",
                                                "nvFl", "vFl") ~ mat1,
                                 varname %in% c("sigmaeps",
                                                "sigmaepsa","sigmaepsl")~ cons), varname1 = paste0(varname, varid2, "-", out))

  var1$varname1
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
mytraceplot <- function(stanres, dirname, filename, by = 5, labels = NULL) {
  start1 <- dim(stanres$fit)[1]
  res1 <- unlstan(stanres$fit, by1 = by, start = start1)
#chain iters var    value
  # <fct> <int> <chr>  <dbl>
  #   1 1     50000 G[1,1] 0.610
  # 2 1     50000 G[2,1] 1.28
  fn1 <-  here::here(file.path(dirname,
                         paste0(filename, "-traceplot.pdf")))

  cols <- RColorBrewer::brewer.pal(4, "Dark2")
  pdf(fn1)
  vars <- unique(res1$var)
  chains <- unique(res1$chain)
  iters1 <- unique(res1$iters)
  iters1 <- seq(ceiling(max(iters1) / 2), max(iters1), by = by)
  res1 <- dplyr::filter(res1, iters %in% iters1)

  par(mfrow = c(2, 2))
  for(i in 1 : length(vars)) {
    res0 <- dplyr::filter(res1, var == vars[i])
    min1 <- min(res0$value, na.rm = T)
    max1 <- max(res0$value, na.rm = T)
    if(!is.null(labels)) {
      lab <- getnames(vars[i], labels)

    } else {
      lab <- vars[i]
    }
    for(j in 1 : length(chains)) {
      res2 <- dplyr::filter(res0, chain == chains[j]) # %>% arrange(iters)
      if(j == 1) {
        plot(res2$iters, res2$value, col = cols[j], type = "l",
             xlab = "Iteration", ylab = lab, main = lab,
             ylim = c(min1, max1))
      } else {
        lines(res2$iters, res2$value, col = cols[j])

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
pairsplot <- function(stanres, dirname, filename, mat1,
                      sources, cons, typeplot, cond1 = "energy") {

  if(typeplot != "ambient") {
    filename1 <- paste0(filename, "-pairs.pdf")

  }else {
    filename1 <- paste0(filename, "-ambient-pairs.pdf")

  }
  pdf(here::here(dirname, filename1))

  # mu G/sigma G
  musig <- c("mug", "sigmag")
  if(typeplot == "local") {
      musig <- paste0(musig, "l")
      append1 <- "l"
      n1 <- 3
  } else if(typeplot == "ambient") {
      musig <- paste0(musig, "a")
      append1 <- "a"
      n1 <- 3
  }else {
    append1 <- ""
    n1 <- 2
  }
  lab1ms <- paste0(rep(musig, each = length(sources)), rep(sources, 2))
  pairs(stanres$fit, pars = musig, labels = lab1ms, condition = cond1)

  # All vf
  types <- musig
  # only vF exactly
  #nF <- names(stanres$fit)[grep("vF", names(stanres$fit))] %>% length()
  nF <- substr(names(stanres$fit), 1, n1)
  vfn <- paste0("vF", append1)
  nF <- length(which(nF == vfn))
  nF1 <- ceiling(nF / 6)
  for(j in 1 : 2) {
    k <- 1
    first5 <- paste0(rep(types[j], length(sources)), sources)
    #print(c("first5", first5))

    # number of pages/plots
    for(i in 1 : nF1) {
      l <- min(c((k + 5), nF))
      #print(c(k, l))
      lab1 <- c(first5, mat1[k : l])
      #print(lab1)
      pairs(stanres$fit, labels = lab1,
            pars = c(types[j], paste0(vfn, "[", (k : l), "]")), condition = cond1)

      k <- k + 6
    }
  }

  lab1 <- paste0("sigmaeps-", cons)
  se1 <- c(paste0("sigmaeps", append1))
  pairs(stanres$fit, labels = lab1,
        pars = se1, condition = cond1)


  # add pairs of sigmaeps with sigmag
  pairs(stanres$fit, pars = c(musig, se1), labels = c(lag1ms, lab1), condition = cond1)


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
pairsjoint <- function(stanres, dirname, filename) {

  filename1 <- paste0(filename, "-pairs-joint.pdf")

  pdf(here::here(dirname, filename1))

  pairs(stanres$fit, pars = musig, labels = lab1, condition = "energy")

  # All vf
  types <- musig
  # only vF exactly
  #nF <- names(stanres$fit)[grep("vF", names(stanres$fit))] %>% length()
  nF <- substr(names(stanres$fit), 1, n1)
  vfn <- paste0("vF", append1)
  nF <- length(which(nF == vfn))
  nF1 <- ceiling(nF / 6)
  for(j in 1 : 2) {
    k <- 1
    first5 <- paste0(rep(types[j], length(sources)), sources)
    #print(c("first5", first5))

    # number of pages/plots
    for(i in 1 : nF1) {
      l <- min(c((k + 5), nF))
      #print(c(k, l))
      lab1 <- c(first5, mat1[k : l])
      #print(lab1)
      pairs(stanres$fit, labels = lab1,
            pars = c(types[j], paste0(vfn, "[", (k : l), "]")), condition = "energy")

      k <- k + 6
    }
  }

  lab1 <- paste0("sigmaeps-", cons)
  pairs(stanres$fit, labels = lab1,
        pars = c(paste0("sigmaeps", append1)), condition = "energy")

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
#' @param labels label dataset
#' @export
energyplot <- function(stanres, dirname, filename, pdf = F,
                       ht = 150, wd = 20, labels = labels, typeplot = "none") {

  # Which items to save
  start1 <- dim(stanres$fit)[1] *2
  nmax <- start1 - start1 * 3/4
  cn <- dim(stanres$fit)[2]

  if(typeplot != "none") {
    pars1 <- c("Ga", "Gl",
                "vFa", "vFl",
                "sigmaga", "sigmagl",
               "muga", "mugl",
                "sigmaepsa", "sigmaepsl")
  } else {
    pars1 <- c("G",
               "vF",
               "sigmag",
               "mug",
               "sigmaeps")
  }

  energy <- rstan::get_sampler_params(stanres$fit)
  names(energy) <- paste0("chain", seq(1, cn))

  # Get energy for later iterations
  energy1 <- lapply(energy, function(x) {
    data.frame(x) %>%
      dplyr::select(., energy__) %>%
      dplyr::mutate(., iters = seq(1, nrow(x))) %>%
      dplyr::filter(., iters >= nmax )
  }) %>%
    tibble::as_tibble_col(.) %>%
    tibble::rowid_to_column(., "chain") %>%
    tidyr::unnest(., value)

  # fix names
  params <- rstan::extract(stanres$fit, permuted = F, pars = pars1) %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(., "iters") %>%
    tidyr::pivot_longer(., names_to = "var", values_to = "val", -iters)

  # Get chain info (for subsetting)
  params1 <- params %>%
    dplyr::mutate(., chain = as.numeric(substr(var, 7, 7)), var = substring(var, 8))

  # Remove some iterations
  energy2 <- dplyr::mutate(energy1, iters = iters - min(iters)) %>%
    dplyr::full_join(., params1)
  # %>%
  #   dplyr::mutate(varname = getnames(var, labels))


  #  highest correlations
  cors <- dplyr::group_by(energy2, chain, var) %>%
    dplyr::summarize(., cor1 = cor(energy__, val)) %>%
    dplyr::mutate(., abscor = abs(cor1)) %>%
    dplyr::arrange(., desc(abscor)) %>%
    dplyr::slice(1 : 5)

  # which to plot
  # samps <- sample(seq(1, nmax), 500, replace = F)
  samps <- seq(1, nmax, length = 1000) %>%
    round() %>% unique()

  # keep chain = 1
  energy3 <- dplyr::filter(energy2, iters %in% samps, chain == 1) %>%
    dplyr::mutate(., chain = factor(chain))
  # could add lables for var here

  # plot
  g1 <- ggplot2::ggplot(energy3, ggplot2::aes(x = val, y = energy__, colour = chain)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~var, scales = "free_x", ncol = 8)

  if(pdf) {
    filename1 <- paste0(filename, "-energy.pdf")

    ggplot2::ggsave(here::here(file.path(dirname, filename1)),
           g1, height = ht, width = wd, limitsize = F)
  }

  list(g1 = g1, cors = cors)
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
rhatplot <- function(stanres, dirname, filename, pdf = F, ht = 120) {

  # Get Rhat and plot
  rhats <- bayesplot::rhat(stanres$fit)
  g1 <-  bayesplot::mcmc_rhat(rhats) + bayesplot::yaxis_text(hjust = 1)

  # If save
  if(pdf) {
    filename1 <- paste0(filename, "-rhat.pdf")

    ggplot2::ggsave(here::here(file.path(dirname, filename1)),
           g1, height = ht, limitsize = F)
  }

  g1

}




#' \code{biasplot} Create energy plots for stan output
#'
#' @title biasplot
#'
#' @param stanres Output from runstan
#' @param dirname Folder name where to save plots
#' @param filename Base filename
#' @param by Skips for parameter posterior draws (boxplot)
#' @param pdf Whether to save pdf, defaults to F
#' @param ht Height of pdf
#' @param wd Width of pdf
#' @param labels labels from makelabels
#' @export
biasplot <- function(stanres, dirname, filename,
                     mat1,
  prof = prof, meansd = meansd, typesim = typesim,
  by = 5,
  pdf = F, ht = 10, wd = 10, labels = labels, typeplot) {


  if(typeplot == "local") {
      append1 <- "l"
  } else if(typeplot == "ambient"){
      append1 <- "a"
  }else {
    append1 <- ""
  }
  # get scaling info
  dat <- stanres$dat
  mean1 <- as.matrix(dat$g[, -1]) %*% as.matrix(dat$f[, -1])
  ynosd <- mean1 + dat$err
  sdscale <- apply(ynosd, 2, sd)
  sddf <- data.frame(sdscale, poll = colnames(dat$y[, -1]))
  prof1 <- dplyr::filter(prof, constraint == 1) %>% dplyr::select(poll, source)
  sdsource <-   dplyr::left_join(prof1, sddf) %>%
    dplyr::select(source, sdscale) %>% na.omit() %>% dplyr::rename(sd1 = sdscale) %>%unique()


  # Get data for boxplot of posterior
  iters <- dim(stanres$fit)[1]
  sel <- seq(1, iters, by = by)

  pars1 <- paste0(c("mug", "sigmag", "sigmaeps", "G", "vF"), append1)

  params <- rstan::extract(stanres$fit, pars = pars1, permuted = F) %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(., "iters") %>%
    #select(contains("chain:1")) %>% # needed if more than 1 chain?
    dplyr::slice(sel) %>%
    tidyr::pivot_longer(., names_to = "var", values_to = "val", -iters) %>%
    dplyr::mutate(., var = substring(var, 9)) #chain = as.numeric(substr(var, 7, 7)),
    #dplyr::filter(chain == 1, iters %in% sel) %>%


  params1 <- dplyr::select(params, var) %>% unique() %>%
    dplyr::mutate(var1 = var) %>%
    tidyr::separate(var1, c("var1", "rowcol"), "\\[") %>%
    tidyr::separate(rowcol, c("row", "col"), ",") %>%
    dplyr::mutate(col = gsub("\\]", "", col),
           row = gsub("\\]", "", row),
           row = as.numeric(row),
           col = as.numeric(col))
    #dplyr::filter(!(var1 %in% c(paste0("nvF", append1), "lp__")))

  params <- dplyr::full_join(params, params1) %>%
    dplyr::select(-c(iters, var))

  params <- params %>%
    dplyr::group_by(var1, row, col) %>%
    dplyr::summarize(q25 = quantile(val, probs = 0.25),
                     med = median(val),
                     q75 = quantile(val, probs = 0.75),
                     mean = mean(val)) %>%
    dplyr::ungroup()

  # Find posterior mean
  # means <- params
  # means <- dplyr::group_by(params) %>%
  #   dplyr::group_by(var1, row, col) %>%
  #   dplyr::summarize(mean = mean(val))

  # truth
  truth <- stanres$dat
  g <- dplyr::rename(truth$g, row = id) %>%
    tidyr::pivot_longer(-row, values_to = "truth")  %>%
    dplyr::mutate(col = as.numeric(factor(name)), var1 = paste0("G", append1), source = name) %>%
    dplyr::full_join(sdsource) %>%
    dplyr::mutate(truth = truth / sd1) %>% dplyr::select(-c(source, sd1))
  f <- dplyr::filter(prof, type == typesim, is.na(constraint)) %>%
    dplyr::select(poll, source, scale1) %>%
    dplyr::full_join(sddf) %>% dplyr::full_join(sdsource) %>%
    dplyr::rename(truth = scale1) %>%
    tidyr::unite(name, c(source, poll), sep = "-") %>%
    dplyr::mutate(row = factor(name, levels = mat1, labels = seq(1, length(mat1))),
           row = as.numeric(row),
           var1 = paste0("vF", append1),   truth = truth / sdscale * sd1) %>%
    dplyr::select(-sdscale) %>% na.omit()

  if(typeplot == "local" & ("mx" %in% names(stanres))) {
    f <- dplyr::filter(f, name %in% stanres$mx) %>%
      dplyr::mutate(name = factor(name, levels = stanres$mx)) %>%
      dplyr::arrange(name) %>%
      dplyr::mutate(name = as.character(name)) %>%
      dplyr::select(-row) %>%
      tibble::rowid_to_column(., "row")

    mat1 <- stanres$mx
  }

  musigg <- dplyr::filter(meansd, type == typesim) %>%
    dplyr::select(-type) %>%
    dplyr::mutate(row = as.numeric(factor(source))) %>%

    # new
    dplyr::mutate(meang = exp(mean + sd^2/2), sdg = sqrt((exp(sd^2) - 1) * meang^2)) %>%
    dplyr::select(-c(mean, sd)) %>%
    tidyr::pivot_longer(-c(source, row),
                 names_to = "var1", values_to = "truth") %>%

    # old
    # tidyr::pivot_longer(-c(source, row),
    #              names_to = "var1", values_to = "truth") %>%
    # dplyr::mutate(var1 = factor(var1, levels = c("mean", "sd"),
    #                      labels = c("mug", "sigmag"))) %>%
    dplyr::full_join(sdsource) %>%
    dplyr::mutate(truth = truth / sd1) %>%


    # for lognormal? i not work
    # dplyr::mutate(truth = ifelse(var1 == "mug", truth + log(sd1)), truth) %>%
    dplyr::select(-sd1) %>%
    dplyr::rename(name = source) %>%
    dplyr::mutate(var1 = ifelse(var1 %in% c("meang", "sdg"), paste0(var1, append1), var1))
  P <- unique(prof$poll) %>% length()

  se1 <- apply(mean1, 2, sd) / 10
  se1 <- se1 / sdscale

  # sigmaeps is variance
  sigmaeps <- data.frame(row = seq(1, P), var1 = paste0("sigmaeps", append1), truth = (se1)^2)
  truth <- dplyr::full_join(g, f) %>%
    dplyr::full_join(musigg) %>% dplyr::full_join(sigmaeps)

  # Find posterior mean

  gval <- dplyr::filter(params, substr(var1, 1, 1) == "G") %>%
    dplyr::group_by(col) %>%
    dplyr::summarize(sdg= sd(mean), meang = mean(mean)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mug = log(meang^2/sqrt(sdg^2 + meang^2)), sigmag = (log(1 + sdg^2/meang^2))) %>%
    tidyr::pivot_longer(-c(col),
                        names_to = "var1", values_to = "mean") %>%
    dplyr::rename(row = col)

  gval1 <- dplyr::filter(gval, var1 %in% c("meang", "sdg")) %>%
    dplyr::mutate(var1 = ifelse(var1 %in% c("meang", "sdg"), paste0(var1, append1), var1))
  gval2 <- dplyr::filter(gval, !(var1 %in% c("meang", "sdg"))) %>%
    dplyr::rename(truth = mean) %>%
    dplyr::mutate(var1 = ifelse(var1 %in% c("mug", "sigmag"), paste0(var1, append1), var1))




  # combine
  dat <- dplyr::full_join(params, truth) %>%
    dplyr::full_join(gval1) %>%
    dplyr::full_join(gval2) %>%
    #dplyr::full_join(means) %>%
    dplyr::mutate(varid = ifelse(var1 == paste0("G", append1), col, row),
           varid2 = row)

  sources <- labels$source[dat$varid]
  mat1 <- labels$mat1[dat$varid]
  cons <- labels$cons[dat$varid]


  if(typeplot == "local" & ("mx" %in% names(stanres))) {

     mat1 <- stanres$mx[dat$varid]
  }

  dat <- dplyr::mutate(dat, sources = sources, mat1 = mat1, cons = cons,

                 out = dplyr::case_when(var1 %in% c("sigmag", "mug", "G",
                                                       "sigmaga", "muga", "Ga",
                                                       "sigmagl", "mugl", "Gl", "sdg", "meang",
                                                    "sdgl", "sdga", "meangl", "meanga") ~ sources,
                                        var1 %in% c("nvF", "vF", "vFl",
                                                       "nvFa", "vFa", "nvFl") ~ mat1,
                                        var1 %in% c("sigmaeps",
                                                       "sigmaepsa","sigmaepsl")~ cons),
                varname1 = paste0(var1, varid2, "-", out))


  # all but G
  append2 <- ifelse(append1 == "l", "a", "l")
  dat1 <- dplyr::filter(dat, !(var1 %in% c(paste0("G", append1), paste0("mug", append2),
                                           paste0("sigmag", append2),
                                           paste0("sigmaeps", append2)) | substr(var1, 1, 2) %in% c("lG", "nv")))
  g1 <- ggplot2::ggplot(data = dat1) +
    # ggplot2::geom_boxplot(ggplot2::aes(x = varname1, y = val)) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = q25)) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = med)) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = q75)) +

    ggplot2::geom_point(ggplot2::aes(x = varname1, y = mean), colour = "red",
               shape = 17) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = truth), colour = "blue",
               shape = 8) +
    ggplot2::facet_wrap(~var1, scales = "free") +
    ggplot2::ylab("Value") +
    ggplot2::theme(axis.text.x =  ggplot2::element_text(angle = 90))


  # G
  dat1 <- dplyr::filter(dat, var1 == paste0("G", append1))
  g2 <- ggplot2::ggplot(data = dat1) +
    # ggplot2::geom_boxplot(ggplot2::aes(x = varname1, y = val)) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = q25)) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = med)) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = q75)) +

    ggplot2::geom_point(ggplot2::aes(x = varname1, y = mean), colour = "red",
               shape = 17) +
    ggplot2::geom_point(ggplot2::aes(x = varname1, y = truth), colour = "blue",
               shape = 8) +
    ggplot2::facet_wrap(~col, scales = "free") +
    ggplot2::ylab("Value") +
    ggplot2::theme(axis.text.x =  ggplot2::element_text(angle = 90))

  # If save
  if(pdf) {
    if(typeplot != "ambient") {
      filename1 <- paste0(filename, "-biasplot-noG.pdf")
    } else {
      filename1 <- paste0(filename, "-biasplot-ambient-noG.pdf")

    }

    ggplot2::ggsave(here::here(file.path(dirname, filename1)),
           g1, height = ht, width = wd, limitsize = F)

    if(typeplot != "ambient") {
      filename1 <- paste0(filename, "-biasplot-G.pdf")
    } else {
      filename1 <- paste0(filename, "-biasplot-ambient-G.pdf")

    }

    ggplot2::ggsave(here::here(file.path(dirname, filename1)),
           g2, height = ht, width = wd, limitsize = F)
  }

  list(g1, g2)

}



#' \code{rescoverage} Calculate coverage of results
#'
#' @title rescoverage
#' @param typesim Which simulation to run.  Currently, local1 - local4, ambient.
#' @param stanres Output from runstan
#' @param by Skips for traceplot, biasplot
#' @param prof Profile data, defaults to prof.
#' @param meansd Matrix with columns corresponding to source, type, mean, sd, defaults to meansd.
#' @export
rescoverage <- function(typesim, stanres, by = 5,
                     prof = prof, meansd = meansd) {

  sources <- stanres$dat$f$source
  mat1 <- mat1fun(stanres, sources)


  # Get simulation-specific mean/SD, profiles
  prof <- dplyr::filter(prof, type == typesim)
  meansd <- dplyr::filter(meansd, type == typesim)

  # get scaling info
  dat <- stanres$dat
  mean1 <- as.matrix(dat$g[, -1]) %*% as.matrix(dat$f[, -1])
  ynosd <- mean1 + dat$err
  sdscale <- apply(ynosd, 2, sd)
  sddf <- data.frame(sdscale, poll = colnames(dat$y[, -1]))
  prof1 <- dplyr::filter(prof, constraint == 1) %>% dplyr::select(poll, source)
  sdsource <-   dplyr::left_join(prof1, sddf) %>%
    dplyr::select(source, sdscale) %>% na.omit() %>% dplyr::rename(sd1 = sdscale)


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
    tidyr::separate(var, c("var1", "rowcol"), "\\[") %>%
    tidyr::separate(rowcol, c("row", "col"), ",") %>%
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
    dplyr::mutate(col = as.numeric(factor(name)), var1 = "G", source = name) %>%
    dplyr::full_join(sdsource) %>%
    dplyr::mutate(truth = truth / sd1) %>% dplyr::select(-c(source, sd1))
  f <- dplyr::filter(prof, type == typesim, is.na(constraint)) %>%
    dplyr::select(poll, source, scale1) %>%
    dplyr::full_join(sddf) %>% dplyr::full_join(sdsource) %>%
    dplyr::rename(truth = scale1) %>%
    tidyr::unite(name, c(source, poll), sep = "-") %>%
    dplyr::mutate(row = factor(name, levels = mat1, labels = seq(1, length(mat1))),
                  row = as.numeric(row),
                  var1 = "vF",   truth = truth / sdscale * sd1) %>%
    dplyr::select(-sdscale)
  musigg <- dplyr::filter(meansd, type == typesim) %>%
    dplyr::select(-type) %>%
    dplyr::mutate(row = as.numeric(factor(source))) %>%
    tidyr::pivot_longer(-c(source, row),
                        names_to = "var1", values_to = "truth") %>%
    dplyr::mutate(var1 = factor(var1, levels = c("mean", "sd"),
                                labels = c("mug", "sigmag"))) %>%
    dplyr::full_join(sdsource) %>%
    dplyr::mutate(truth = truth / sd1) %>% dplyr::select(-sd1) %>%
    dplyr::rename(name = source)
  P <- unique(prof$poll) %>% length()

  se1 <- apply(mean1, 2, sd) / 10
  se1 <- se1 / sdscale


  sigmaeps <- data.frame(row = seq(1, P), var1 = "sigmaeps", truth = se1)
  truth <- dplyr::full_join(g, f) %>%
    dplyr::full_join(musigg) %>% dplyr::full_join(sigmaeps)

  # Find posterior mean
  means <- dplyr::group_by(params) %>%
    dplyr::group_by(var1, row, col) %>%
    dplyr::summarize(lb = quantile(val, probs = 0.025),
                     ub = quantile(val, probs = 0.975), postmean = mean(val))

  # combine
  dat <- dplyr::full_join(means, truth) %>%
    dplyr::select(-sd1) %>%
    dplyr::mutate(coverage = ifelse(truth >= lb & truth <= ub, 1, 0))
  dat

}
