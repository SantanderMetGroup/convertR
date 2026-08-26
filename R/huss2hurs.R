##     huss2hurs.R Specific humidity to relative humidity conversion
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


#' @title Specific humidity to relative humidity conversion
#' @description Specific humidity to relative humidity conversion using surface pressure and temperature
#' @param huss Specific humidity grid
#' @param ps Surface pressure
#' @param tas Near-surface air temperature
#' @return A climate4R grid
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
#' hurs <- huss2hurs(huss.iberia, ps.iberia, tas.iberia)
#' \dontrun{
#' require(visualizeR)
#' spatialPlot(climatology(hurs), rev.colors = TRUE, backdrop.theme = "coastline")
#' }

huss2hurs <- function(huss, ps, tas) {
    # Consistency checks:
    if (isMultigrid(huss) | isMultigrid(ps) | isMultigrid(tas)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(huss) | isGrid(ps) | isGrid(tas))
    # Redim to have members:
    huss  %<>%  redim(member = TRUE)
    ps %<>%  redim(member = TRUE)
    tas  %<>%  redim(member = TRUE)
    # Check dim
    suppressMessages(checkDim(huss, ps, tas))
    # Check season
    checkSeason(huss, ps, tas)
    # Check units
    u.huss <- getGridUnits(huss)
    if (u.huss != "kg/kg") {
        if (!ud.are.convertible(u.huss, "kg/kg")) {
            stop("Non compliant huss units (should be convertible to \'kg/kg\')")
        }
        message("[", Sys.time(), "] Converting units ...")
        huss <- udConvertGrid(huss, new.units = "kg/kg")
    }
    huss %<>% redim(member = TRUE)
    message("[", Sys.time(), "] Calculating Relative humidity ...")
    ws <- suppressMessages(tas2ws(tas, ps)) %>% redim(member = TRUE)
    tas <- NULL
    coords <- getCoordinates(huss)
    n.mem <- getShape(huss, "member")
    l <- lapply(1:n.mem, function(x) {
        a <- subsetGrid(huss, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        b <- subsetGrid(ws, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        w <- a/(1 - a) # w = mixing (mass) ratio
        a <- 100 * w/b
        a[a > 100] <- 100
        a[a < 0] <- 0
        ps$Data <- mat2Dto3Darray(a, x = coords$x, y = coords$y)
        a <- b <- NULL
        return(ps)
    })
    huss <- ws <- ps <- NULL
    hurs <- suppressWarnings(bindGrid(l, dimension = "member"))
    hurs$Variable$varName <- "hurs"
    hurs$Variable$level <- NULL
    attr(hurs$Variable, "units") <- "%"
    attr(hurs$Variable, "longname") <- "Surface_air_relative_humidity"
    attr(hurs$Variable, "description") <- "Estimated relative humidity from saturation pressure and specific humidity"
    attr(hurs, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(hurs, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(hurs)
}
