########## .display_round ##########

#' @title Round for human-facing display, without silently showing exact zeros
#' @description Rounds \code{x} to \code{digits} decimal places, as \code{round()} does, except
#' for elements that would otherwise display as \code{0} (i.e.\ round to zero while the original
#' value is non-zero); those elements fall back to \code{signif(x, digits)} so that a small but
#' genuinely non-zero estimate is not misread as exactly zero in printed summaries
#' @param x Numeric vector, matrix, or array to round for display
#' @param digits Integer indicating the number of decimal places (or significant digits, for the
#' zero-rounding fallback). Default \code{3}
#' @keywords internal
.display_round <- function(x, digits = 3)
{
  r <- round(x, digits)
  zeroed <- !is.na(x) & x != 0 & r == 0
  if (any(zeroed)) r[zeroed] <- signif(x[zeroed], digits)
  r
}
