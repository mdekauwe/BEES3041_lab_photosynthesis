####
##  General funcs that the photosynthesis model depends on
##
##  author: Martin De Kauwe
##  date: 13th March, 2019
##  email: mdekauwe@gmail.com
####

is_close <- function(a, b, rel_tol=1e-09, abs_tol=0.0) {
  #
  #   Check if two values are close ...
  #
  #   Args:
  #   -----
  #   a : float
  #     value to compare
  #   b : float
  #     value to compare
  #   rel_tol : float
  #     relative tolerance, max allowed diff btw a and b
  #   abs_tol : float
  #     minimum absolute tolerance, must be at least zero
  #   c : float
  #     co-efficient
  #
  #   Returns:
  #   -------
  #   boolean : logical
  #     Return True if the values a and b are close to each other and False
  #     otherwise.

  return ( abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol) )

}

quadratic <- function(a, b, c, large=FALSE) {
  #
  #   Minimilist quadratic solution as root for J solution should always
  #   be positive, so I have excluded other quadratic solution steps. I am
  #   only returning the smallest of the two roots
  #
  #   Args:
  #   -----
  #   a : float
  #     co-efficient
  #   b : float
  #     co-efficient
  #   c : float
  #     co-efficient
  #
  #   Returns:
  #   -------
  #   val : float
  #     positive root
  #

  # discriminant
  d <- b**2.0 - 4.0 * a * c

  if ( any(d < 0.0) ) {
    stop("imaginary root found")
  }

  if (large) {

    if ( any((is_close(a, 0.0)) & (b > 0.0)) ) {
      root <- -c / b
    } else if ( any((is_close(a, 0.0)) & (is_close(b, 0.0))) ) {
        root <- 0.0
        if (c != 0.0) {
          stop("Cant solve quadratic")
        }
    } else {
        root <- (-b + sqrt(d)) / (2.0 * a)
    }

  } else {

    if ( any((is_close(a, 0.0)) & (b > 0.0)) ) {
      root <- -c / b
    } else if ( any((is_close(a, 0.0)) & (is_close(b, 0.0))) ) {
      root <- 0.0
      if (c != 0.0) {
        stop('Cant solve quadratic')
      }
    } else {
      root <- (-b - sqrt(d)) / (2.0 * a)
    }

  }

  return ( root )

}
