msa <- ">Homo sapiens\nMKAANSE\n>Mouse\nMKA-NSE"

test_that("the props JSON keeps arrays of one and integer columns", {
  w <- msaview(msa = msa) +
    geom_msa_highlight(start = 2, end = 4, label = "motif") +
    geom_msa_highlight(columns = 7)
  json <- as.character(props_json(w$x$props))
  expect_match(json, '"highlights":[{', fixed = TRUE)
  expect_match(json, '"highlightColumns":[7]', fixed = TRUE)
})

test_that("a file type the CLI cannot write is refused before it runs", {
  expect_error(msaview_save(msaview(msa = msa), "fig.jpg", cli = "false"),
               "must end in")
  expect_error(msaview_save(list(), "fig.svg"), "msaview()")
})

test_that("a failing CLI stops with its output", {
  script <- tempfile(fileext = ".R")
  writeLines('cat("broken\\n", file = stderr()); quit(status = 3)', script)
  rscript <- file.path(R.home("bin"), "Rscript")
  expect_error(
    msaview_save(msaview(msa = msa), tempfile(fileext = ".svg"),
                 cli = c(rscript, shQuote(script))),
    "broken"
  )
})

# MSAVIEW_CLI names the built CLI, e.g. packages/cli/dist/index.js
test_that("the saved figure draws the widget's layers", {
  cli <- Sys.getenv("MSAVIEW_CLI")
  skip_if(cli == "", "MSAVIEW_CLI is unset")
  w <- msaview(msa = msa) +
    geom_msa_highlight(start = 2, end = 4, label = "motif")
  svg <- tempfile(fileext = ".svg")
  msaview_save(w, svg, cli = c("node", cli))
  text <- paste(readLines(svg, warn = FALSE), collapse = "\n")
  expect_match(text, "motif", fixed = TRUE)
  expect_match(text, ">Mouse<", fixed = TRUE)
})

test_that("selection_only passes --selection, and needs a selection", {
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "args <- commandArgs(TRUE)",
    "out <- args[which(args == \"-o\") + 1]",
    "writeLines(paste(args, collapse = \" \"), out)"
  ), script)
  rscript <- file.path(R.home("bin"), "Rscript")
  w <- msaview(msa = msa) + geom_msa_selection(2, 4)
  svg <- tempfile(fileext = ".svg")
  msaview_save(w, svg, selection_only = TRUE, cli = c(rscript, shQuote(script)))
  expect_match(readLines(svg), "--selection", fixed = TRUE)
  msaview_save(w, svg, cli = c(rscript, shQuote(script)))
  expect_no_match(readLines(svg), "--selection", fixed = TRUE)
  expect_error(
    msaview_save(msaview(msa = msa), svg, selection_only = TRUE, cli = "false"),
    "no selection"
  )
})

test_that("a selection-only figure draws the selected rows alone", {
  cli <- Sys.getenv("MSAVIEW_CLI")
  skip_if(cli == "", "MSAVIEW_CLI is unset")
  w <- msaview(msa = msa) + geom_msa_selection(2, 4, rows = "Mouse")
  svg <- tempfile(fileext = ".svg")
  msaview_save(w, svg, selection_only = TRUE, cli = c("node", cli))
  text <- paste(readLines(svg, warn = FALSE), collapse = "\n")
  expect_match(text, ">Mouse<", fixed = TRUE)
  expect_no_match(text, ">Homo_sapiens<", fixed = TRUE)
})
