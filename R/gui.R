#' SSDM package Graphic User Interface
#'
#' User interface of the SSDM package.
#'
#' @param port char. The TCP port that the application should listen on (see
#'   \code{\link[shiny]{runApp}} for more details).
#' @param host char. The IPv4 address that the application should listen on (see
#'   \code{\link[shiny]{runApp}} for more details).
#' @param working.directory char. Directory in which the application will run.
#'
#' @return Open a window with a shiny app to use the SSDM package with an
#'   user-friendly interface.
#'
#' @details If your environmental variables have an important size, you should
#'   give enough memory to the interface with the (\code{maxmem} parameter).
#'   Note that only one instance of gui can be run at a time.
#'
#' @examples
#' \dontrun{
#' gui()
#' }
#'
#' @export
gui <- function (port = getOption("shiny.port"),
                 host = getOption("shiny.host", "127.0.0.1"),
                 working.directory = getwd()) {

  appDir <- system.file("shiny", "gui", package = "SSDM")
  if (appDir == "") {
    stop("Could not find shiny directory. Try re-installing `SSDM`.", call. = FALSE)
  }

  ui <- source(file.path(appDir, 'ui.R'))
  serverWD <- source(file.path(appDir, 'server.R'))
  shiny::runApp(shinyApp(ui = ui, server = serverWD(working.directory)),
                display.mode = "normal", port = port, host = host)
  rm(ui, serverWD, envir = .GlobalEnv)
}
