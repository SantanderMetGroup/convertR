##     udConvertGrid.R Grid unit conversion
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

#' @title Unit conversion
#' @description Unit conversion for climate4R objects
#' @param grid Input grid (either regular, rotated or station data).
#' @param new.units Character string vector indicating a udunits-parseable unit definition. See Details.
#' @return Returns (invisibly) a grid similar to the input, but with new units
#' @details
#'
#' \strong{Unit conversion and consistency checks}
#'
#' This function is a wrapper of the \code{\link[udunits2]{ud.convert}} function from the UCAR's UDUNITS
#' R bindings, tailored to the climate4R CDM. The function performs a number of checks in order to ensure that: (i) the
#' units of the input grid and the new units are recognised by the UDUNITS database and (ii), that
#' the conversion between input and output units is possible. This is internally achieved by the functions
#' (i) \code{\link[udunits2]{ud.is.parseable}} and (ii) \code{\link[udunits2]{ud.are.convertible}}
#'  from package \pkg{udunits2}.
#'
#' \strong{Multigrid support}
#'
#' The function supports multigrids as input. Note that in this case, the length of the \code{new.units} vector should match the
#' number of variables within the multigrid (shape of \code{'var'} dimension), and preserve its ordering (otherwise failing in
#' the above-mentioned consistency checks).
#'
#' @import transformeR
#' @importFrom udunits2 ud.is.parseable ud.are.convertible ud.convert
#' @importFrom magrittr %<>%
#' @seealso \code{\link[transformeR]{getGridUnits}}, for accessing the \code{"units"} attribute of a grid and
#' \code{\link[transformeR]{setGridUnits}}, to manually modify it.
#' @export
#' @author J Bedia
#' @references
#'
#' \url{https://www.unidata.ucar.edu/software/udunits}
#'
#' @examples
#' library(transformeR)
#' require(climate4R.datasets)
#' data("NCEP_Iberia_ta850")
#' # Converting from Kelvin to Celsius:
#' getGridUnits(NCEP_Iberia_ta850)
#' range(NCEP_Iberia_ta850$Data, na.rm = TRUE)
#' ta850_degC <- udConvertGrid(grid = NCEP_Iberia_ta850, new.units = "celsius")
#' getGridUnits(ta850_degC)
#' range(ta850_degC$Data, na.rm = TRUE)
#' # Converting from Pascals to mmHg:
#' data("NCEP_Iberia_psl")
#' getGridUnits(NCEP_Iberia_psl)
#' range(NCEP_Iberia_psl$Data, na.rm = TRUE)
#' psl_mmHg <- udConvertGrid(NCEP_Iberia_psl, new.units = "mmHg")
#' getGridUnits(psl_mmHg)
#' range(psl_mmHg$Data, na.rm = TRUE)
#' # Dealing with multigrids:
#' mg <- makeMultiGrid(NCEP_Iberia_psl, NCEP_Iberia_ta850)
#' getVarNames(mg)
#' getGridUnits(mg)
#' # Convert all variables:
#' mg.conv <- udConvertGrid(mg, new.units = c("mmHg", "celsius"))
#' getGridUnits(mg.conv)

udConvertGrid <- function(grid, new.units) {
    stopifnot(isGrid(grid))
    u1 <- suppressMessages(getGridUnits(grid))
    grid %<>% redim(var = TRUE)
    nvar <- getShape(grid, "var")
    varnames <- getVarNames(grid)
    if (nvar != length(new.units)) stop("Inconsistent number of variables an \'new.units\' vector length")
    l <- lapply(1:nvar, function(x) {
        gr <- suppressMessages(subsetGrid(grid, var = varnames[x]))
        if (!ud.is.parseable(u1[x])) stop("Non-parseable unit string in input grid (\'ud.is.parseable\' returned FALSE)")
        if (!ud.is.parseable(new.units[x])) stop("Non-parseable \'new.units\' string (\'ud.is.parseable\' returned FALSE)")
        if (!ud.are.convertible(u1[x], new.units[x])) stop("Non-convertible units \'ud.are.convertible\' returned FALSE)")
        gr$Data %<>%  ud.convert(u1 = u1[x], u2 = new.units[x])
        gr %<>% setGridUnits(unit.string = new.units[x])
    })
    grid <- suppressWarnings(makeMultiGrid(l))
    invisible(grid)
}
