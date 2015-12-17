#' SSDM package Global User Interface
#'
#' User interface of the SSDM package.
#'
#' @return Open a window with a shiny app to use the SSDM package with an
#'   user-friendly interface.
#'
#' @details If your environmental variables have an important size, you should
#'   gave enough memory to the interface with the (\code{maxmem} parameter).
#'
#' @examples
#' \dontrun{
#' gui()
#' }
#'
#' @export
gui = function () {
  appDir <- system.file("shiny", "gui", package = "SSDM")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
