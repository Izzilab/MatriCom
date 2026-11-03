#' Launch MatriCom locally
#'
#' Starts the MatriCom Shiny application on this computer. Files of any size
#' can be uploaded; analysis runs entirely locally.
#'
#' @param launch.browser Open the app in the default web browser.
#' @param port Port to listen on (default: chosen by Shiny).
#' @param ... Further arguments passed to [shiny::runApp()].
#' @return Runs until the app is closed.
#' @export
matricom.app <- function(launch.browser = TRUE, port = getOption("shiny.port"), ...) {
  app_dir <- system.file("webApp", package = "MatriCom")
  if (!nzchar(app_dir)) stop("MatriCom app files not found; reinstall the package.")
  shiny::runApp(app_dir, launch.browser = launch.browser, port = port, ...)
}
