bundle_hint <- paste0(
  "ggmsaview was installed without its JavaScript. Install it with ",
  "remotes::install_github(\"GMOD/ggmsaview\"), or in a checkout run ",
  "`node scripts/build-js.mjs` beside a built JBrowseMSA checkout."
)

has_bundle <- function(path) {
  nzchar(system.file(path, package = "ggmsaview"))
}

.onAttach <- function(libname, pkgname) {
  if (!has_bundle("htmlwidgets/lib/react-msaview.umd.js")) {
    packageStartupMessage(bundle_hint)
  }
}
