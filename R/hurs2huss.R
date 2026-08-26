##     hurs2huss.R Relative humidity to specific humidity conversion
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


#' @title Relative humidity to specific humidity conversion
#' @description Relative humidity to specific humidity conversion using surface pressure and temperature
#' @param hurs Relative humidity grid
#' @param ps Surface pressure
#' @param tas Near-surface air temperature
#' @return A climate4R CDM grid of specific humidity (in kg/kg)
#' @author J. Bedia, S. Herrera
#' @template templateUnits
#' @export
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
#' hurs <- huss2hurs(huss = huss.iberia, ps = ps.iberia, tas = tas.iberia)
#' huss <- hurs2huss(hurs = hurs, ps = ps.iberia, tas = tas.iberia)
#' identical(huss$Data, huss.iberia$Data)
#' # Not identical due to rounding errors:
#' range(huss$Data - huss.iberia$Data)
#' \dontrun{
#' require(visualizeR)
#' spatialPlot(climatology(hurs), rev.colors = TRUE, backdrop.theme = "coastline")
#' }

hurs2huss <- function(hurs, ps, tas) {
    w <- hurs2w(hurs = hurs, ps = ps, tas = tas) %>% redim(member = TRUE)
    ps <- tas <- NULL
    coords <- getCoordinates(w)
    n.mem <- getShape(w, "member")
    l <- lapply(1:n.mem, function(x) {
        a <- subsetGrid(w, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        a <- a / (1 + a)
        hurs$Data <- mat2Dto3Darray(a, x = coords$x, y = coords$y)
        return(hurs)
    })
    hurs <- w <- NULL
    huss <- suppressWarnings(bindGrid(l, dimension = "member"))
    huss$Variable$varName <- "huss"
    huss$Variable$level <- NULL
    attr(huss$Variable, "units") <- "kg/kg"
    attr(huss$Variable, "longname") <- "specific_humidity"
    attr(huss$Variable, "description") <- "The ratio of the mass of water vapor to the total mass of the moist air parcel"
    attr(huss, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(huss, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(huss)
}
