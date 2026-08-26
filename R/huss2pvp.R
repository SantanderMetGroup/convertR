##     huss2pvp.R Specific humidity to partial vapor pressure
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


#' @title Specific humidity to partial vapor pressure conversion
#' @description Specific humidity to partial vapor pressure conversion using surface pressure and specific humidity
#' @param huss Specific humidity grid
#' @param ps Surface pressure
#' @return A climate4R CDM grid of partial vapor pressure (in Pascals)
#' @author J. Bedia, S. Herrera
#' @template templateUnits
#' @export
#' @import transformeR
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @importFrom utils packageVersion
#' @family derivation
#' @family humidity
#' @family pressure
#' @template templateRefHumidity
#' @details The partial pressure of a gaseous constituent of air (water vapor in this case)
#'  is the pressure which it alone would exert with unchanged temperature and number of moles per unit volume.
#' @examples
#' data("ps.iberia")
#' data("huss.iberia")
#' pvp <- huss2pvp(huss = huss.iberia, ps = ps.iberia)
#' \dontrun{
#' require(visualizeR)
#' spatialPlot(climatology(pvp), rev.colors = TRUE, backdrop.theme = "coastline")
#' }

huss2pvp <- function(huss, ps) {
    # Consistency checks:
    if (isMultigrid(huss) | isMultigrid(ps)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(huss) | isGrid(ps))
    # Redim to have members:
    huss  %<>%  redim(member = TRUE)
    ps %<>%  redim(member = TRUE)
    # Check dim
    suppressMessages(checkDim(huss, ps))
    # Check season
    checkSeason(huss, ps)
    # Check units
    u.huss <- getGridUnits(huss)
    u.ps <- getGridUnits(ps)
    if (u.huss != "kg/kg") {
        if (!ud.are.convertible(u.huss, "kg/kg")) {
            stop("Non compliant huss units (should be convertible to \'kg/kg\')")
        }
        message("[", Sys.time(), "] Converting units ...")
        huss <- udConvertGrid(huss, new.units = "kg/kg")
    }
    if (u.ps != "Pa") {
        if (!ud.are.convertible(u.ps, "Pa")) {
            stop("Non compliant huss units (should be convertible to \'Pa\')")
        }
        message("[", Sys.time(), "] Converting units ...")
        ps <- udConvertGrid(ps, new.units = "Pa")
    }
    huss %<>% redim(member = TRUE)
    ps %<>% redim(member = TRUE)
    message("[", Sys.time(), "] Calculating partial vapour pressure ...")
    source(file.path(find.package(package = "convertR"), "constants.R"), local = TRUE)
    pvp <- ps
    coords <- getCoordinates(huss)
    n.mem <- getShape(huss, "member")
    l <- lapply(1:n.mem, function(x) {
        h <- subsetGrid(huss, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        p <- subsetGrid(ps, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        p <- h * p / (epsilon * (1 - h) + h)
        pvp$Data <- mat2Dto3Darray(p, x = coords$x, y = coords$y)
        return(pvp)
    })
    huss <- ps <- NULL
    pvp <- suppressWarnings(bindGrid(l, dimension = "member"))
    pvp$Variable$varName <- "pvp"
    pvp$Variable$level <- NULL
    attr(pvp$Variable, "units") <- "Pa"
    attr(pvp$Variable, "longname") <- "water_vapor_partial_pressure_in_air"
    attr(pvp$Variable, "description") <- "Estimated partial vapour pressure"
    attr(pvp, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(pvp, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(pvp)
}

