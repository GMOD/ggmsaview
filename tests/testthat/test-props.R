# msaview() wraps MSAViewer, and a prop it cannot reach is a gap a reader finds
# by not finding it. scripts/build-js.mjs writes the props MSAViewer takes into
# inst/msaviewer-props.json, beside the bundle built from the same commit.

viewer <- function() {
  file <- system.file("msaviewer-props.json", package = "ggmsaview")
  skip_if(!nzchar(file), "scripts/build-js.mjs has not run")
  jsonlite::fromJSON(file)
}

# deliberate absences, each with its reason; a prop in neither list fails
not_exposed <- c(
  onCellHover = "a message per pointer move floods a Shiny websocket; the widget reports clicks and the viewport instead",
  onModel = "hands back a live MST model, which no Shiny input can carry"
)
# callbacks the widget turns into Shiny inputs rather than props
as_events <- c(onCellClick = "_click", onViewportChange = "_viewport",
               onSelectionChange = "_selection")

# the props$x <- y keys msaview() builds, plus height, which an htmlwidget
# sizes itself with
reached <- function() {
  code <- paste(deparse(msaview), collapse = "\n")
  built <- regmatches(code, gregexpr("props\\$\\w+ <-", code))[[1]]
  c(sub("props\\$(\\w+) <-", "\\1", built), "height")
}

test_that("the prop list is readable, and covers all of MSAViewerProps", {
  props <- viewer()$props
  expect_true(all(c("msa", "residueEncoding") %in% props))
  expect_gt(length(props), 20)
})

test_that("every MSAViewer prop reaches msaview()", {
  props <- viewer()$props
  missing <- setdiff(props, c(reached(), names(not_exposed), names(as_events)))
  expect_equal(missing, character(0),
    info = "add the argument to msaview() in R/msaview.R, assign it in the props list, document it with @param, and rerun roxygen")
})

test_that("msaview() passes no prop MSAViewer does not have", {
  expect_equal(setdiff(reached(), viewer()$props), character(0))
})

test_that("the callbacks the widget turns into Shiny inputs are wired", {
  js <- paste(readLines(system.file("htmlwidgets", "msaview.js",
                                    package = "ggmsaview")), collapse = "\n")
  expect_true(all(vapply(as_events, grepl, logical(1), x = js, fixed = TRUE)))
})

test_that("a prop left out of msaview() is still one MSAViewer has", {
  expect_equal(setdiff(names(not_exposed), viewer()$props), character(0))
})
