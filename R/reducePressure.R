#' @title Reduce pressure level
#' @description Internal function with formula for presssure reduction
#' @param p A 2D matrix of (sea-level/surface) pressure values (in Pascals)
#' @param t A 2D matrix of near-surface temperature (in K)
#' @param direction Direction of the conversion, either \code{"ps2psl"} or \code{"psl2ps"}.
#' @param zgs A matrix of geopotential height (static variable, all rows are equal, in meters)
#' @return A 2D matrix with reduced pressure data. Need to be reconverted to a grid using \code{\link[transformeR]{mat2Dto3Darray}} afterwards.
#' @details This function is a internal for the derivation functions \code{\link{ps2psl}} and \code{\link{psl2ps}}.
#' @template templateRefPressure
#' @keywords internal
#' @seealso \code{\link{psl2ps}} and \code{\link{ps2psl}}, the exported functions using this one

reducePressure <- function(p, t, zgs, direction) {
    direction <- match.arg(direction, choices = c("psl2ps", "ps2psl"))
    source(file.path(find.package(package = "convertR"), "constants.R"), local = TRUE)
    ind <- which(abs(zgs) >= .001)
    auxGamma <- NA * p
    To <- t + GammaST * zgs / g
    ind1 <- intersect(which(To > 290.5 & t <= 290.5), ind)
    auxGamma[ind1] <- g * (290.5 - t[ind1]) / zgs[ind1]
    ind <- setdiff(ind, ind1)
    ind1 <- intersect(which(To > 290.5 & t > 290.5), ind)
    auxGamma[ind1] <- 0
    t[ind1] <- 0.5 * (255 + t[ind1])
    ind <- setdiff(ind,ind1)
    ind1 <- intersect(which(t < 255), ind)
    auxGamma[ind1] <- GammaST
    t[ind1] <- 0.5 * (255 + t[ind1])
    ind <- setdiff(ind, ind1)
    auxGamma[ind] <- GammaST
    ind <- which(abs(zgs) >= .001)
    term <- switch(direction,
                   "psl2ps" = -zgs[ind],
                   "ps2psl" = zgs[ind])
    paux <- p
    paux[ind] <- p[ind] * exp((term / (Rd * t[ind])) * (1 - 0.5 * (auxGamma[ind] * zgs[ind]) / (g * t[ind]) + (1/3) * ((auxGamma[ind]*zgs[ind])/(g*t[ind]))^2))
    return(paux)
}
