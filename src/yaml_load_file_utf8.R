# This function is lifted from rstudio/rmarkdown R/output_format.R
# the yaml UTF-8 bug has been fixed https://github.com/viking/r-yaml/issues/6
# but yaml >= 2.1.14 Win/Mac binaries are not available for R < 3.2.0, so we
# still need the mark_utf8 trick
#' @importFrom utils packageVersion
yaml_load_utf8 <- function(string, ...) {
  string <- paste(string, collapse = '\n')
  if (packageVersion('yaml') >= '2.1.14') {
    yaml::yaml.load(string, ...)
  } else {
    mark_utf8(yaml::yaml.load(enc2utf8(string), ...))
  }
}

yaml_load_file_utf8 <- function(input, ...) {
  yaml_load_utf8(readLines(input, encoding = 'UTF-8'), ...)
}
