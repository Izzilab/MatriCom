\name{matricom.app}
\alias{matricom.app}
\title{Launch MatriCom locally}
\usage{
matricom.app(launch.browser = TRUE, port = getOption("shiny.port"), ...)
}
\arguments{
\item{launch.browser}{Open the app in the default web browser.}
\item{port}{Port to listen on (default: chosen by Shiny).}
\item{...}{Further arguments passed to \code{\link[shiny]{runApp}}.}
}
\value{Runs until the app is closed.}
\description{
Starts the MatriCom Shiny application on this computer. Files of any size can
be uploaded; analysis runs entirely locally.
}
