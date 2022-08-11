#' Accept family arguments in the same way as glm(), convert to standardized
#' form, and throw comprehensible errors.
#'
#' @importFrom rlang caller_arg
#' @keywords internal
normalize_family <- function(family, arg = caller_arg(family)) {
  msg <- "{.arg {arg}} must be a family object, function returning a family, or the name of a family"
  call <- caller_env()

  # based on the first few lines of glm()
  tryCatch({
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame(n = 2))
    }
  }, error = function(e) {
    cli_abort(
      c(msg, "x" = "family {.val {family}} is not a known family"),
      call = call
    )
  })

  tryCatch({
    if (is.function(family)) {
      family <- family()
    }
  }, error = function(e) {
    cli_abort(c(msg, "x" = "value provided is a function that does not return a family"),
              call = call)
  })

  if (!inherits(family, "family")) {
    cli_abort(c(msg, "x" = "family provided has class {.cls {class(family)}}"),
              call = call)
  }

  return(family)
}
