#' Profile values with constraints
#' 
#' A dataset of PM2.5 constituent profiles for sources
#' in US EPA Region 3 from Ivey et al (2017).  These are
#' matched to four simulation study types.
#'  The variables are:
#' 
#' \itemize{
#' 		\item poll. PM2.5 chemical constituent
#' 		\item source. Source name
#' 		\item val. Concentration of pollutant in source in mug/m^3
#' 		\item type. Simulation type (ambient, local1 - local4)
#' 		\item constraint. Constraint (0 or 1) for source apportionment
#' 		\item scale1. Profile value scaled to the constraint = 1 pollutant
#'		...
#'}
#'
#' @docType data
#' @keywords datasets
#' @format A tibble with 281 rows and 6 variables
#' @name prof
NULL