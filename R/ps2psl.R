##     ps2psl.R Surface pressure reduction to sea-level
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


#' @title Surface pressure reduction to sea-level
#' @description Estimate the sea-level pressure from surface pressure, using temperature and surface geopotential height
#' @param ps Surface pressure
#' @param tas Near-surface air temperature
#' @param zgs surface geopotential height. This is a static variable (i.e., constant in time).
#'  Thus, member and time dimensions must be either singletons or just dropped.
#' @return A climate4R CDM grid of estimated sea-level pressure (in Pascals)
#' @author J. Bedia, S. Herrera
#' @template templateUnits
#' @export
#' @template templateRefPressure
#' @import transformeR
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @importFrom utils packageVersion
#' @family derivation
#' @family pressure
#' @examples
#' data("ps.iberia")
#' data("tas.iberia")
#' data("zgs.iberia")
#' psl.hat <- ps2psl(ps = ps.iberia, tas = tas.iberia, zgs = zgs.iberia)
#' require(transformeR)
#' # Error map:
#' psl.err <- gridArithmetics(psl.hat, psl.iberia, operator = "-")
#' \dontrun{
#' require(visualizeR)
#' err <- climatology(psl.err)
#' spatialPlot(err, set.min = -10000,
#'             at = seq(-10000, 10000, 500),
#'             color.theme = "Spectral", rev.colors = TRUE,
#'             backdrop.theme = "coastline",
#'             main = "Mean estimate error (Pascals)")
#' }

ps2psl <- function(ps, tas, zgs) {
    # Consistency checks:
    if (isMultigrid(ps) | isMultigrid(tas) | isMultigrid(zgs)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(ps) | isGrid(tas) | isGrid(zgs))
    # Redim to have members:
    ps %<>% redim(member = TRUE)
    tas %<>% redim(member = TRUE)
    # Check dim
    suppressMessages(checkDim(ps, tas))
    suppressMessages(checkDim(ps, zgs, dimensions = c("lat", "lon")))
    # Check season
    checkSeason(tas, ps)
    # Check units
    ps.u <- getGridUnits(ps)
    if (ps.u != "Pa") {
        if (!ud.are.convertible(ps.u, "Pa")) {
            stop("Non compliant ps units (should be convertible to Pascals)")
        }
        message("[", Sys.time(), "] Converting pressure units ...")
        ps %<>% udConvertGrid(new.units = "Pa") %>% redim(member = TRUE)
    }
    tas.u <- getGridUnits(tas)
    if (tas.u != "K") {
        if (!ud.are.convertible(tas.u, "K")) {
            stop("Non compliant tas units (should be convertible to Kelvin)")
        }
        message("[", Sys.time(), "] Converting temperature units ...")
        tas %<>% udConvertGrid(new.units = "K") %>% redim(member = TRUE)
    }
    zgs.u <- getGridUnits(zgs)
    if (zgs.u != "m") {
        if (!ud.are.convertible(zgs.u, "m")) {
            stop("Non compliant zgs units (should be convertible to meters)")
        }
        message("[", Sys.time(), "] Converting surface geopotential height units ...")
        zgs %<>% udConvertGrid(new.units = "m") %>% redim(drop = TRUE)
    }
    coords <- getCoordinates(ps)
    npoints <- do.call("expand.grid", coords) %>% nrow()
    ntimes <- getShape(ps, "time")
    zgs %<>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat() %>% as.vector() %>% matrix(nrow = ntimes, ncol = npoints, byrow = TRUE)
    n.mem <- getShape(ps, "member")
    message("[", Sys.time(), "] Reducing pressure level ...")
    source(file.path(find.package(package = "convertR"), "constants.R"), local = TRUE)
    psl <- ps
    l <- lapply(1:n.mem, function(x) {
        p <- subsetGrid(ps, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        t <- subsetGrid(tas, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        aux <- reducePressure(p = p, t = t, zgs = zgs, direction = "ps2psl")
        psl$Data <- mat2Dto3Darray(aux, coords$x, coords$y)
        return(psl)
    })
    psl <- suppressWarnings(bindGrid(l, dimension = "member"))
    tas <- ps <- zgs <- NULL
    psl$Variable$varName <- "psl"
    psl$Variable$level <- NULL
    attr(psl$Variable, "units") <- "Pa"
    attr(psl$Variable, "longname") <- "air_pressure_at_mean_sea_level"
    attr(psl$Variable, "description") <- "Estimated sea-level pressure from pressure at surface"
    attr(psl, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(psl, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(psl)
}

