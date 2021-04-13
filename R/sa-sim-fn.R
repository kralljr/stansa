

#' \code{simdat} Simulate all data
#'
#' @title simdat
#'
#' @param stanmodel Type of stan model to run.  Currently only noninformative tested.
#' @param typesim Which simulation to run.  Currently, local1 - local4, ambient.
#' @param N Number of days to simulate.  If not noninformative stan, vector of 2 for days in ambient, days in local.
#' @param prof Profile data, defaults to prof
#' @param meansd Matrix with columns corresponding to source, type, mean, sd, defaults to meansd.
#' @param sderr Standard deviation of error. If not noninformative stan, optional vector of 2 for error in ambient, error in local.
#' @param rmout Whether to remove largest (outliers)
#' @export
#'
simdat <- function(stanmodel, typesim, N, prof = prof, meansd = meansd,
  sderr = NULL, rmout = F) {

  outl <- simdat1(typesim, N[1],  prof, meansd, sderr[1], rmout)

  # If not non-informative, need both ambient and local data
  if(stanmodel %in% c("joint", "penalty")) {
    warning("Joint/Penalty: not yet tested!")
    sderra <- ifelse(length(sderr) == 1, sderr, sderr[2])
    outa <- simdat1("ambient", N[2], prof, meansd, sderra, rmout)

    out <- reorgout(outa, outl)


  # just One dataset
  } else {
    out <- outl
  }

  return(out)
}







#' \code{simdat1} Inner functions to simulate each data (e.g., ambient or local)
#'
#' @title simdat1
#'
#' @param typesim Which simulation to run.  Currently, local1 - local4, ambient.
#' @param N Number of days to simulate.
#' @param prof Profile data, defaults to prof.
#' @param meansd Matrix with columns corresponding to source, type, mean, sd, defaults to meansd.
#' @param sderr Standard deviation of error.
#' @param rmout Whether to remove largest (outliers)
#' @export
simdat1 <- function(typesim = "ambient", N = 100, prof = prof,
                    meansd = meansd, sderr = NULL, rmout = F) {

  ## Format
  # format G means
  meansd <- dplyr::filter(meansd, type == typesim) %>%
    dplyr::select(., -type)


  # Format profiles: Keep 1 scale
  prof1 <- dplyr::filter(prof, type == typesim) %>%
    dplyr::select(., poll, source, scale1, constraint)
  f1 <- dplyr::select(prof1, -constraint) %>%
    tidyr::pivot_wider(., names_from = "poll", values_from = "scale1")
  cn <- colnames(f1)[-which(colnames(f1) == "source")] %>% sort()
  f1 <- f1[, c("source", cn)]
  # alphabetize by source
  f1 <- dplyr::arrange(f1, source)
  constraints <- dplyr::select(prof1, -scale1) %>%
    tidyr::pivot_wider(., names_from = "poll",
        values_from = "constraint") %>%
    dplyr::arrange(source) %>%
    dplyr::select(., -source)
  constraints <- constraints[, cn]
  zeromat <- (is.na(constraints)) * 1
  onemat <- (!is.na(constraints) & constraints == 1) * 1

  # Get dimensions
  L <- nrow(f1)
  P <- ncol(f1) - 1

  # Get G
  g1 <- geng(meansd, N, rmout)
  #alphabetize by source
  cn <- colnames(g1)[-which(colnames(g1) == "id")] %>% sort()
  g1 <- g1[, c("id", cn)]


  # Get pollutant data
  g0 <- dplyr::select(g1, -id) %>% as.matrix()
  f0 <- f1
  rn <- f0$source
  f0 <- dplyr::select(f1, -source) %>% as.matrix()
  rownames(f0) <- rn

  mean <- g0 %*% f0

  # Get error matrix
  if(!is.null(sderr)) {
    err <- matrix(rnorm(N * P, sd = sderr), ncol = P)
  } else {
    err <- matrix(rep(apply(mean, 2, sd) / 100, each = N), byrow = F, ncol = P)
  }

  #y <- exp(log(g0 %*% f0) + (err))
  y <- mean + err

  # Standardize as in Nikolov
  # sd1 <- apply(y, 2, sd)
  # y <- sweep(y, 2, sd1, "/")

  colnames(y) <- colnames(f0)
  y1 <- data.frame(id = g1$id, y)


  ae <- identical(colnames(y1)[-1], colnames(f0)) & identical(rownames(f0), colnames(g0))
  if(!ae) {
    print(colnames(y))
    print(colnames(f0))
    print(rownames(f0))
    print(colnames(g0))
    stop("Names do not match")
  }
  ###########
  # Other setup

  # Vector of ones for getting row sums
  ones <- rep(1, P)


  # Position of free elements
  pos <- which(zeromat == 1, arr.ind = T)
  posc <- pos[, "col"]
  posr <- pos[, "row"]

  # Number of free elements
  # Position of free elements
  pos <- t(apply(zeromat, 1, function(x) which(x == 1)))

  # Free elements
  B <- ncol(pos)
  LB <- B * L

  # return stan data and truth
  out <- list(stan = list(N = N, L = L, P = P, B = B, LB = LB, posr = posr, posc = posc,
                          y = y, ones = ones, zeromat = zeromat, onemat = onemat),
              true = list(y = y1, g = g1, f = f1, sigmaeps = sderr, err = err))
}



#' \code{genn} Generate normal data from dataset
#'
#' @title genn
#'
#' @param dat Matrix with columns corresponding to mean, sd
#' @param N Number of days to simulate
#' @export
genn <- function(dat, N) {
  data.frame(id = seq(1, N), conc = rnorm(N, dat$mean, dat$sd))
}







#' \code{geng} Simulate G matrix
#'
#' @title geng
#'
#' @param meansd Matrix with columns corresponding to source, mean, sd
#' @param N Number of days to simulate
#' @param rmout Whether to remove largest
#' @export
geng <- function(meansd, N, rmout = F) {

  # Add extra if need
  N1 <- ifelse(rmout, ceiling(1.0025 * N), N)

  # Generate lognormal
  d1 <- tidyr::nest(meansd, data = c(mean, sd)) %>%
    dplyr::mutate(rand = purrr::map(data, ~genn(., N1))) %>%
    tidyr::unnest(rand) %>%
    dplyr::select(., -data)

  # Remove largest
  if(rmout) {
    keeptop <- function(dat) {
      ords <- order(dat$conc, decreasing = T)[1 : (N1 - N)]
      dplyr::slice(dat, -ords)
    }
    dl <- tidyr::nest(d1, data = c(id, conc)) %>%
      dplyr::mutate(rand = purrr::map(data, ~keeptop(.))) %>%
      tidyr::unnest(rand) %>%
      dplyr::select(., -data)
  }


  d1 <- tidyr::pivot_wider(d1, names_from = "source",
                    values_from = "conc")
  d1

}










#' \code{checkc2} Check rank condition C2 for Bayesian SA
#'
#' @title checkc2
#'
#' @param dat Data frame of NA (unconstrained), 0s, and 1s, consisting of a column of pollutants called poll and columns corresponding to constraints for each source
#' @export
checkc2 <- function(dat) {

  # Add in random normal for unconstrained and convert to matrix
  mat <- tidyr::pivot_longer(dat, -poll, names_to = "source", values_to = "val") %>%
    dplyr::mutate(., val = ifelse(is.na(val), rnorm(n()), val)) %>%
    tidyr::pivot_wider(., names_from = "poll", values_from = "val") %>%
    dplyr::select(., -source) %>% as.matrix()

  # Rank of each submatrix equal to nsources - 1
  ranks <- vector()
  for(i in 1 : nrow(mat)) {
    ranks[i] <- (as.numeric(Matrix::rankMatrix(mat[-i, which(mat[i, ] == 0)])) == (nrow(mat) - 1)) * 1
  }
  # Return whether condition is met
  ifelse(nrow(mat) - sum(ranks) == 0, "condition met", "fail")
}








#' \code{reorgout} Reorganize local/ambient output
#'
#' @title reorgout
#'
#' @param outa Ambient data output from simdat1
#' @param outl Local data output from simdat1
reorgout <- function(outa, outl) {
  warning("Not tested!")


  # add local/ambient names
  names(outa$stan) <- paste0(names(outa$stan), "a")
  names(outl$stan) <- paste0(names(outl$stan), "l")

  names(outa$true) <- paste0(names(outa$true), "a")
  names(outl$true) <- paste0(names(outl$true), "l")

  # get source names
  sourcea <- outa$true$f$source
  sourcel <- outl$true$f$source

  # get pollutant names
  polla <- select(outa$true$y, -id) %>% colnames()
  polll <- select(outl$true$y, -id) %>% colnames()
  #all pollutants
  Ps <- unique(c(polla, polll))
  Ptot <- length(Ps)

  # which local sources in ambient (use their information) vs. which in local
  LsM <- which(sourcel %in% sourcea)
  LsN <- which(!(sourcel %in% sourcea))
  Lmatch <- length(LsM)
  Lnomatch <- length(LsN)
  LsN <- ifelse(Lnomatch == 0, 0, LsN)

  # get total sources
  Ls <- unique(c(sourcea, sourcel))
  Ltot <- length(Ls)
  Lsa <- match(sourcea, Ls)
  # remove no matching sources (not necessary for Fstartot)
  sourcelmatch <- sourcel[sourcel %in% sourcea]
  Lsl <- match(sourcelmatch, Ls)

  BslnotA <- which(!(polll %in% polla))
  BlnotA <- length(BslnotA)

  # Links to total
  La <- outa$stan$La
  Ba <- outa$stan$Ba
  Bsa <- matrix(nrow = La, ncol = Ba)
  hold <- match(polla, Ps)

  for(l in 1 : La) {
    for(b in 1 : Ba) {

      pos <- outa$stan$posa[l, ]
      Bsa[l, ] <- hold[pos]

      # check
      ae <- identical(polla[pos], Ps[hold[pos]])
      #print(ae)
      if(!ae) {
        stop("Ambient matching off")
      }
    }
  }

  # free except those not in ambient
  BlA <- outl$stan$Bl - BlnotA
  Bsl <- matrix(nrow = Lmatch, ncol = BlA)
  posl2 <- matrix(nrow = Lmatch, ncol = BlA)

  hold <- match(polll, Ps)

  for(l in 1 : Lmatch) {
    for(b in 1 : BlA) {
      pos <- outl$stan$posl[LsM[l], ]
      # Which matching pollutants ambient/local
      pos <- pos[!(pos %in% BslnotA)]
      posl2[l, ] <- pos
      Bsl[l, ] <- hold[pos]

      # check
      ae <- identical(polll[pos], Ps[hold[pos]])
      if(!ae) {
        stop("Local matching off")
      }
    }
  }


  # Fix arrays
  LsM <- as.array(LsM)
  LsN <- as.array(LsN)
  BslnotA <- as.array(BslnotA)
  Lsl <- as.array(Lsl)

  if(Lmatch == 1) {
    stop("Code cannot handle 1 matching source for joint currently")
  }

  # format output
  extra <- list(Lmatch = Lmatch, LsM = LsM,
                Ltot = Ltot, Lsa = Lsa, Lsl = Lsl,
                BlnotA = BlnotA, BslnotA = BslnotA, BlA = BlA,
                posl2 = posl2, Ptot = Ptot, Bsa = Bsa, Bsl = Bsl)

  if(Lnomatch > 0) {

    extra <- append(extra, list(Lnomatch = Lnomatch, LsN = LsN))
  }

  stan <- append(outa$stan, outl$stan) %>% append(., extra)
  true <- append(outa$true, outl$true)
  out <- list(stan = stan, true = true)

}

