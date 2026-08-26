#' @description A climate4R CDM grid containing daily NCEP Reanalysis data for the Iberian Peninsula (1991 and 1992, complete monthly series).
#' @details This R data object was obtained using the function \code{loadGridData} from \pkg{loadeR}, from the User Data Gateway,
#'  using the following parameters:
#'
#'  \code{loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/ncepReanalysis1/ncepReanalysis1_4xDaily.ncml",
#'                     var = var,
#'                     years = 1991:1992,
#'                     season = 1:12,
#'                     lonLim = c(-10, 5),
#'                     latLim <- c(35, 43),
#'                     time = "DD")
#'       }
#'
#' , with the appropriate \code{aggr.m} argument value depending on the cell method applied to \code{var} for daily aggregation.
#' @references
#' \url{http://www.meteo.unican.es/climate4r}
