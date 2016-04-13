
#' Function to run stan using lambda constraints
#' 
#' \code{SAstan} creates .stan file for Source apportionment with constraints
#'
#' @param wd directory for basic stan code
#' @param dat matrix of days (or commutes) by constituents
#' @param lamcon matrix of sources by constituents with NA for unconstrained, 0's and 1's as in Park et al. 2002
#' @param file1 output file name
#' @export
SAstan <- function(wd, dat, lamcon, file1 = "acesa-out.stan") {
  identifySA(lamcon)

  # Add in transformed data
  trans <- paste0("transformed data {\n", rowcol(lamcon), "}\n")
  # Add in constraints
  model1 <- modelF(lamcon)


  stancode <- readLines(file.path(wd, "acesa.stan"))
  # find where to split
  li <- grep("//hold", stancode)

  # Cat all info together
  cat(stancode[1 : (li[1] - 1)], trans, 
      stancode[(li[1] + 1) : (li[2] - 1)], model1, 
      stancode[(li[2] + 1) : length(stancode)],
      file = file1, sep = "\n")

}


#' Function to create row_mark and col_mark as in http://bit.ly/1T3Vj2W
#' 
#' \code{rowcol} Creates stan inputs based on constraint matrix lamcon
#' 
#' @param lamcon matrix of sources by constituents with NA for unconstrained, 0's and 1's as in Park et al. 2002
#' @return pasted character vector to be placed in .stan file giving information on row and col marks
rowcol <- function(lamcon) {
  whAR <- which(is.na(lamcon), arr.ind = T)
  K <- nrow(whAR)
  L <- nrow(lamcon)
  P <- ncol(lamcon)
  rm1 <- "row_mark["
  cm1 <- "col_mark["

  # Specify integers
  out <- paste0("  int<lower=1,upper=", L, "> row_mark[", K, "];\n")
  out2 <- paste0("  int<lower=1,upper=", P, "> col_mark[", K, "];\n")
  out <- paste0(out, out2)
  # Specify each element
  for(k in 1 : K) {
    rmhold <- paste0(rm1, k, "] <- ", whAR[k, 1], "; ")
    cmhold <- paste0(cm1, k, "] <- ", whAR[k, 2], ";") 
    out <- paste0(out, "  ", rmhold, cmhold, "\n")
  }

  return(out)
}


#' Incorporate fixed constraints in STAN
#' 
#' \code{modelF} writes stan code based on lambda constraints for SA
#' 
#' @param lamcon matrix of sources by constituents with NA for unconstrained, 0's and 1's as in Park et al. 2002
#' @return pasted character vector to be placed in .stan file giving information on constraints
modelF <- function(lamcon) {
  whAR <- which(!is.na(lamcon), arr.ind = T)

  # Add in model info
  start0 <- "  F["
  mid0 <- "] <- "
  model <- ""
  for(i in 1 : nrow(whAR)) {
    keep1 <- whAR[i, ]
    index1 <- paste(keep1, collapse = ", ")
    val1 <- lamcon[keep1[1], keep1[2]]
    model <- paste0(model, start0, index1, mid0, val1, ";\n") 
  }
  model
}

