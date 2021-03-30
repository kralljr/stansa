



#' \code{runstan} Simulate all data
#'
#' @title runstan

#' @param N Number of days to simulate.  If not noninformative stan, vector of 2 for days in ambient, days in local
#' @param typesim Which simulation to run.  Currently, local1 - local4, ambient.
#' @param stancode Type of stan model to run
#' @param stantype Type of stan model (can differ from code: for ifelse statements)
#' @param prof Profile data, defaults to prof.
#' @param meansd Matrix with columns corresponding to source, type, mean, sd, defaults to meansd.
#' @param keepall Whether to output all posterior draws and simulated data.
#' @param rmout Whether to remove largest (outliers)
#' @param sderr Standard deviation of error.  If not noninformative stan, can be vector of 2 for ambient local
#' @param seed Seeds for stan data simulation and running (length 2).  If inform, seed must be length 3
#' @param iter Number of iterations to run
#' @param chains Number of chains
#' @param findamb Results from find ambient for informative
#' @export
runstan <- function(N, typesim, stancode = NULL, stantype = NULL,
                    prof = prof, meansd = meansd, keepall = T,
                    rmout = F, sderr = 0.01, seeds = NULL, iter = 1000,
                    chains = 1, findamb = NULL, notes = NULL, fp = NULL, ...) {

  options(warn = 1)
  con <- file(here("logs/warnings.log"))
  sink(con, append=TRUE)
  sink(con, append=TRUE, type="message")
  cl <- paste0(match.call(), collapse = " ")
  print(cl)
  starttime <- Sys.time()

  if(is.null(seeds)) {
    seeds <- sample(seq(1, 10000), 3)
  }


  set.seed(seeds[1])

  # Type of stancode to run
  stantype <- ifelse(!is.null(stancode) & is.null(stantype),
                     sapply(strsplit(substring(stancode, 9), "\\."), function(x) x[[1]]),
                     ifelse(!is.null(stantype), stantype, "none"))

  # all sources match
  if(typesim != "local4" & stantype == "joint") {
    stancode <- "stan-sa-joint-allmatch.stan"
  }

  # simulate data
  dat <- simdat(stantype, typesim, N, prof, meansd, sderr, rmout)

  if(stantype == "inform") {
    # Note type sim is still for local simulation
    warning("Not tested!")
    res <- findamb(dat, profdat, sderr, rmout, typesim, N, iter = iter,
                   chain = chains, seed = seeds[3], ...)
    # add
    dat$stan <- append(dat$stan, )

  }

  # run stan model
  if(!is.null(stancode)) {

    # starttime1 <- stringr::str_replace_all(starttime, " ", "-") %>%
    #   stringr::str_replace_all(., ":", "")
    starttime1 <- gsub(" ", "-", starttime) %>% gsub(":", "", .)
    starttime1 <-
    loggr::log_file(
      here(paste0("logs/mainstan-", starttime1,".log"))
      , subscriptions=c('message', 'warning','stop'),
      log_muffled = T
    )
    fit <- rstan::stan(file= stancode, data= dat$stan,
                seed = seeds[2], iter = iter, chains = chains,
                ...)
  }

  # get output
  if(keepall) {
    out <- list(dat = dat$true, standat = dat$stan, fit = fit)
  } else if(!is.null(stancode)) {
    out <- list(data = dat$true, standat = dat$stan,
                summary = summary(fit)$summary, stantype = stantype)
  } else {
    out <- dat$true
  }

  endtime <- Sys.time()

  # save relevant metadata
  sources <- dat$true$g %>% colnames()

  meta <- data.frame(starttime = starttime, endtime = endtime, fp = fp, notes = notes, typesim = typesim, stancode = stancode,
                     sderr = sderr, N = N, seeds = paste(seeds, collapse = ", "),
                     iter = iter,  chains = chains, call = cl, comments = "")
  coln <- ifelse(file.exists(here("logs/sim-model-log.csv")), F, T)
  write.table(meta, file = here("logs/sim-model-log.csv"),
              sep = ",", row.names = F, col.names = coln, append = T)


  sink()
  sink(type="message")
  return(out)
}



#' \code{findamb} Find ambient priors
#'
#' @title findamb
#'
#' @param dat Simulated local data
#' @param amb Ambient data and stan results + stan data
#' @param profdat Profile data
#' @param meansd Matrix with columns corresponding to source, type, mean, sd
#' @param typesim Which simulation to run
#' @param N Number of days to simulate.  If not noninformative stan, vector of 2 for days in ambient, days in local
#' @param rmout Whether to remove largest
#' @param sderr Standard deviation of error.  If not noninformative stan, can be vector of 2 for ambient local
#' @param seed Seed for stan.  Must be a vector of length 2
#' @export
findamb <- function(dat, amb) {
  warning("Not tested recently!")
  sfit <- amb$summary


  # Need to get meanf, sdf (L by B matrices)
  res <- data.frame(param = rownames(sfit), sfit) %>%
    dplyr::filter(., substr(param, 1, 2) == "Fs") %>%
    dplyr::select(., param, mean, sd) %>%
    tidyr::separate(., param, c("row", "col"), ",") %>%
    dplyr::mutate(., row = as.numeric(substring(row, 7)), col = as.numeric(gsub("\\]", "", col)))
  mean <- dplyr::select(res, -sd) %>%
    tidyr::pivot_wider(., names_from = "col", values_from = "mean") %>%
    dplyr::arrange(., row) %>%
    dplyr::select(., -row)
  sd <- dplyr::select(res, -mean) %>%
    tidyr::pivot_wider(., names_from = "col", values_from =  "sd") %>%
    dplyr::arrange(., row)%>%
    dplyr::select(., -row)

  # get into dimensions of constituents, sources (not free elements)
  pos <- amb$stan$pos
  onemat <- amb$stan$onemat
  zeromat <- amb$stan$zeromat
  mean <- getfmat(mean, pos, onemat, zeromat)
  sd <- getfmat(sd, pos, onemat, zeromat)

  # Get zeromat for local, sources for local
  zeromat <- dat$stan$zeromat
  L <- dat$stan$L
  B <- dat$stan$B
  pos <- dat$stan$pos
  sourcel <- dat$true$f$source
  sourcea <- amb$true$f$source
  polll <- colnames(dat$true$y)[-1]
  polla <- colnames(amb$true$y)[-1]


  # get final meanc/sdc L X B
  meanc <- matrix(nrow = L, ncol = B)
  sdc <- matrix(nrow = L, ncol = B)
  for(l in 1 : L) {

    # if source in ambient
    if(sourcel[l] %in% sourcea) {
      # match source
      source1 <- which(sourcea == sourcel[l])

      # for each free
      for(b in 1 : B) {

        # get pollutant of interest
        poll1 <- polll[pos[l, b]]

        #  if pollutant in ambient
        if(poll1 %in% polla) {
          poll1 <- which(polla == poll1)
          meanc[l, b]<- mean[source1, poll1]
          sdc[l, b] <- sd[source1, poll1]

        # else not inform
        } else {
          meanc[l, b] <- 0
          sdc[l, b] <- 10000
        }

      }


    # uninform otherwise
    } else {
      meanc[l, ] <- rep(0, B)
      sdc[l, ] <- rep(10000, B)
    }

  }


  colnames(meanc) <- NULL
  colnames(sdc) <- NULL
  list(meanf = meanc, sdf = sdc)
}





#' \code{getfmat} Get F matrix from Stan output of free elements
#'
#' @title getfmat
#'
#' @param ffree Fmatrix of free elements from Stan output
#' @param pos position data
#' @param onemat Onemat constraints of 1
#' @param zeromat Zeromat constraints of 0s and 1s
#' @export
getfmat <- function(ffree, pos, onemat, zeromat) {

  warning("Not tested recently!")
  fmat <- zeromat + onemat
  for(i in 1 : nrow(pos)) {
    for(j in 1 : ncol(pos)) {
      fmat[i, pos[i, j]] <- ffree[i, j]
    }
  }
  fmat
}
