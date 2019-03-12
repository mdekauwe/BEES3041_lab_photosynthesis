
is_close <- function(a, b, rel_tol=1e-09, abs_tol=0.0) {
  #
  # Check if two values are close ...
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
