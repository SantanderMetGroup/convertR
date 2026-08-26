##     psl2ps.R Sea-level pressure to surface pressure
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


#' @title Sea-level pressure to surface pressure
#' @description Estimate the surface pressure from sea-level pressure, using temperature and surface geopotential height
#' @param psl Sea-surface pressure
#' @param tas Near-surface air temperature
#' @param zgs surface geopotential height. This is a static variable (i.e., constant in time).
#'  Thus, member and time dimensions must be either singletons or just dropped.
#' @return A climate4R CDM grid of estimated surface pressure (in Pascals)
#' @author J. Bedia, S. Herrera
#' @template templateUnits
#' @export
#' @template templateRefPressure
#' @import transformeR
#' @details Surface pressure is often needed for other variable derivations, but in many models is often not available, while mean sea-level
#' pressure is usually provided.
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom udunits2 ud.are.convertible
#' @importFrom utils packageVersion
#' @family derivation
#' @family pressure
#' @examples
#' data("psl.iberia")
#' data("tas.iberia")
#' data("zgs.iberia")
#' ps.hat <- psl2ps(psl = psl.iberia, tas = tas.iberia, zgs = zgs.iberia)
#' require(transformeR)
#' # Error map:
#' data("ps.iberia")
#' ps.err <- gridArithmetics(ps.hat, ps.iberia, operator = "-")
#' \dontrun{
#' require(visualizeR)
#' err <- climatology(ps.err)
#' spatialPlot(err, set.max = 10000,
#'             at = seq(-10000, 10000, 500),
#'             color.theme = "Spectral", rev.colors = TRUE,
#'             backdrop.theme = "coastline",
#'             main = "Mean estimate error (Pascals)")
#' }

psl2ps <- function(psl, tas, zgs) {
    # Consistency checks:
    if (isMultigrid(psl) | isMultigrid(tas) | isMultigrid(zgs)) stop("Multigrids are not an allowed input")
    stopifnot(isGrid(psl) | isGrid(tas) | isGrid(zgs))
    # Redim to have members:
    psl %<>% redim(member = TRUE)
    tas %<>% redim(member = TRUE)
    # Check dim
    suppressMessages(checkDim(psl, tas))
    suppressMessages(checkDim(psl, zgs, dimensions = c("lat", "lon")))
    # Check season
    checkSeason(tas, psl)
    # Check units
    psl.u <- getGridUnits(psl)
    if (psl.u != "Pa") {
        if (!ud.are.convertible(psl.u, "Pa")) {
            stop("Non compliant psl units (should be convertible to Pascals)")
        }
        message("[", Sys.time(), "] Converting pressure units ...")
        psl %<>% udConvertGrid(new.units = "Pa") %>% redim(member = TRUE)
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
    coords <- getCoordinates(psl)
    npoints <- do.call("expand.grid", coords) %>% nrow()
    ntimes <- getShape(psl, "time")
    zgs %<>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat() %>% as.vector() %>% matrix(nrow = ntimes, ncol = npoints, byrow = TRUE)
    n.mem <- getShape(psl, "member")
    message("[", Sys.time(), "] Reducing pressure level ...")
    source(file.path(find.package(package = "convertR"), "constants.R"), local = TRUE)
    ps <- psl
    l <- lapply(1:n.mem, function(x) {
        p <- subsetGrid(psl, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        t <- subsetGrid(tas, members = x, drop = TRUE) %>% redim(member = FALSE) %>% extract2("Data") %>% array3Dto2Dmat()
        aux <- reducePressure(p = p, t = t, zgs = zgs, direction = "psl2ps")
        ps$Data <- mat2Dto3Darray(aux, coords$x, coords$y)
        return(ps)
    })
    ps <- suppressWarnings(bindGrid(l, dimension = "member"))
    tas <- psl <- zgs <- NULL
    ps$Variable$varName <- "ps"
    ps$Variable$level <- NULL
    attr(ps$Variable, "units") <- "Pa"
    attr(ps$Variable, "longname") <- "surface_pressure"
    attr(ps$Variable, "description") <- "Estimated surface pressure from sea-level pressure"
    attr(ps, "origin") <- paste0("Calculated with R package 'convertR' v", packageVersion("convertR"))
    attr(ps, "URL") <- "https://github.com/SantanderMetGroup/convertR"
    message("[", Sys.time(), "] Done.")
    invisible(ps)
}
