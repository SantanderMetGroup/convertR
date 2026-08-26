##     hurs2w.R Mixing ratio from realtive humidity
##
##     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
##
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @title Mixing (mass) ratio of water vapor
#' @description Calculate the mass ratio of water vapor in air from relative humidity, pressure and temperature
#' @param ps Surface pressure
#' @param tas Near-surface air temperature
#' @param hurs Relative humidity grid
#' @return A climate4R CDM grid
#' @author J. Bedia, S. Herrera
#' @keywords internal
#' @template templateUnits
#' @export
#' @template templateRefHumidity
#' @details
#' The mixing ratio usually refers to the mass ratio (\eqn{\zeta_i}), which is defined as the mass of a constituent \eqn{m_i}
#' divided by the total mass of all other constituents in a mixture:
#'
#' \deqn{\zeta_{i}=\frac{m_i}{m_{tot}-m_{i}}}
#'
#' The mass ratio of water vapor in air can be used to describe humidity (\url{https://en.wikipedia.org/wiki/Mixing_ratiohttps://en.wikipedia.org/wiki/Mixing_ratio}).
#'
#' @import transformeR
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @importFrom utils packageVersion
#' @family derivation
#' @family humidity
#' @template templateRefHumidity
#' @examples
#' data("ps.iberia")
#' data("tas.iberia")
#' data("huss.iberia")
#' # First, hurs is derived from huss:
#' hurs <- huss2hurs(huss.iberia, ps.iberia, tas.iberia)
#' w <- hurs2w(ps.iberia, tas.iberia, hurs)
#' \dontrun{
#' require(visualizeR)
#' spatialPlot(climatology(w), rev.colors = TRUE, backdrop.theme = "coastline")
#' }

hurs2w <- function(ps, tas, hurs) {
    # Consistency checks:
    if (isMultigrid(hurs)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(hurs))
    # Redim to have members:
    hurs %<>% redim(member = TRUE)
    message("[", Sys.time(), "] Calculating saturation pressure...")
    ws <- suppressMessages(tas2ws(tas = tas, ps = ps)) %>% redim(member = TRUE)
    tas <- NULL
    # Check dim
    suppressMessages(checkDim(hurs, ws))
    # Check season
    checkSeason(hurs, ws)
    # Check units
    u1 <- getGridUnits(hurs)
    if (u1 != "%") {
        if (!ud.are.convertible(u1, "%")) {
            stop("Non compliant ps units (should be convertible to \'%')")
        }
        message("[", Sys.time(), "] Converting units ...")
        hurs %<>% udConvertGrid(new.units = "%") %>% redim(member = TRUE)
    }
    coords <- getCoordinates(hurs)
    n.mem <- getShape(hurs, "member")
    l <- lapply(1:n.mem, function(x) {
        h <- subsetGrid(hurs, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        w <- subsetGrid(ws, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        w <- w * h * .01
        ps$Data <- mat2Dto3Darray(w, x = coords$x, y = coords$y)
        return(ps)
    })
    w <- suppressWarnings(bindGrid(l, dimension = "member"))
    hurs <- ps <- ws <- NULL
    w$Variable$varName <- "w"
    w$Variable$level <- NULL
    attr(w$Variable, "units") <- "1"
    attr(w$Variable, "longname") <- "humidity_mixing_ratio"
    attr(w$Variable, "description") <- "Humidity mixing ratio of a parcel of moist air is the ratio of the mass of water vapor to the mass of dry air"
    attr(w, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(w, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(w)
}
