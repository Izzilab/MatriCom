.onAttach <- function(libname, pkgname) {
  packageStartupMessage("MatriCom ", utils::packageVersion("MatriCom"), " - run matricom.app() to start")
}
