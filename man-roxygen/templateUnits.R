#' @section Input units:
#' In principle, the derivation functions will work regardless of the input units, as all units are internally
#'  converted as necessary according to the specific formulae definitions.
#'   However, the units string (see \code{\link[transformeR]{getGridUnits}}) must be parseable. In the same vein,
#'   the input units need to be convertible to the required ones. Unit consistency is internally achieved
#'   by the function \code{\link{udConvertGrid}}.
