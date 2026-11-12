msa <- ">Homo sapiens\nACGTACGT\n>Mouse\nACGTACGA"

test_that("a layer and the matching argument build the same props", {
  variants <- data.frame(
    row = "Homo sapiens", start = 2, end = 4, label = "R248Q",
    stringsAsFactors = FALSE
  )
  composed <- msaview(msa = msa) + geom_msa_highlight(variants)
  passed <- msaview(msa = msa, highlights = list(
    list(row = "Homo sapiens", start = 2, end = 4, label = "R248Q")
  ))
  # each surface keeps the field order it was handed, and a JSON object carries
  # no order, so compare the fields themselves
  by_name <- function(hs) lapply(hs, function(h) h[order(names(h))])
  expect_equal(by_name(composed$x$props$highlights),
               by_name(passed$x$props$highlights))
})

test_that("a highlight data frame gives one highlight per row", {
  variants <- data.frame(start = c(2, 5), end = c(3, 6), label = c("a", "b"))
  w <- msaview(msa = msa) + geom_msa_highlight(variants)
  expect_length(w$x$props$highlights, 2)
  expect_equal(w$x$props$highlights[[2]]$label, "b")
})

test_that("a highlight data frame carries only the columns it has", {
  w <- msaview(msa = msa) +
    geom_msa_highlight(data.frame(start = 2, end = 4))
  expect_equal(names(w$x$props$highlights[[1]]), c("start", "end"))
})

# a named vector is written out as FASTA, which sanitizes the headers, so the
# highlight's row has to take the same treatment to still name a row
test_that("a layer sanitizes a row name the way the alignment does", {
  w <- msaview(msa = c("Homo sapiens" = "ACGTACGT", Mouse = "ACGTACGA")) +
    geom_msa_highlight(data.frame(
      row = "Homo sapiens", start = 2, end = 4, stringsAsFactors = FALSE
    ))
  expect_equal(w$x$props$highlights[[1]]$row, "Homo_sapiens")
  expect_match(w$x$props$msa, ">Homo_sapiens", fixed = TRUE)
})

test_that("a highlight carries its mark from the argument or a column", {
  w <- msaview(msa = msa) +
    geom_msa_highlight(start = 2, end = 4, label = "RBD", mark = "bracket") +
    geom_msa_highlight(data.frame(
      start = 5, end = 6, mark = "bracket", stringsAsFactors = FALSE
    ))
  expect_equal(
    vapply(w$x$props$highlights, function(h) h$mark, ""),
    c("bracket", "bracket")
  )
})

test_that("highlight layers accumulate", {
  w <- msaview(msa = msa) +
    geom_msa_highlight(start = 1, end = 2) +
    geom_msa_highlight(start = 5, end = 6)
  expect_length(w$x$props$highlights, 2)
  expect_equal(w$x$props$highlights[[1]]$start, 1)
  expect_equal(w$x$props$highlights[[2]]$start, 5)
})

test_that("a highlight needs a data frame, a span, or rows", {
  expect_error(geom_msa_highlight(), "needs a data frame")
  expect_error(
    geom_msa_highlight(data.frame(first = 1)),
    "must have 'start' and 'end' columns"
  )
  expect_error(geom_msa_highlight(1:3), "takes a data frame")
})

test_that("whole-row highlights keep rows an array", {
  w <- msaview(msa = msa) + geom_msa_highlight(rows = "Mouse", label = "clade")
  json <- as.character(jsonlite::toJSON(
    w$x$props$highlights[[1]], auto_unbox = TRUE
  ))
  expect_match(json, '"rows":["Mouse"]', fixed = TRUE)
})

test_that("a track reads its values from a data frame column", {
  w <- msaview(msa = msa) +
    geom_msa_track(data.frame(dnds = c(1, 2, 3)), value = "dnds",
                   name = "dN/dS", domain_max = 2)
  track <- w$x$props$columnTracks[[1]]
  expect_equal(as.numeric(track$values), c(1, 2, 3))
  expect_equal(track$name, "dN/dS")
  expect_equal(track$encoding$y$scale$domainMax, 2)
})

test_that("a text track's range is the scale of its letters", {
  w <- msaview(msa = msa) +
    geom_msa_track(text = "12", name = "frame", kind = "text",
                   range = c("1" = "#aaaaaa"))
  expect_match(
    as.character(jsonlite::toJSON(w$x$props$columnTracks[[1]]$encoding,
                                  auto_unbox = TRUE)),
    '{"color":{"scale":{"domain":["1"],"range":["#aaaaaa"]}}}', fixed = TRUE
  )
})

test_that("a bar track with no range or domain_max carries no encoding", {
  w <- msaview(msa = msa) + geom_msa_track(c(1, 2), name = "score")
  expect_null(w$x$props$columnTracks[[1]]$encoding)
})

test_that("a track takes a numeric vector directly", {
  w <- msaview(msa = msa) + geom_msa_track(c(1, 2), name = "score")
  expect_equal(as.numeric(w$x$props$columnTracks[[1]]$values), c(1, 2))
})

test_that("a track id comes from its name", {
  w <- msaview(msa = msa) + geom_msa_track(c(1), name = "dN/dS ratio")
  expect_equal(w$x$props$columnTracks[[1]]$id, "dn-ds-ratio")
})

test_that("a track without a name or a values column is an error", {
  expect_error(geom_msa_track(c(1, 2)), "needs a name")
  expect_error(
    geom_msa_track(data.frame(other = 1), name = "x"),
    "no 'value' column"
  )
})

test_that("an arc track takes its pairs from a data frame", {
  w <- msaview(msa = msa) +
    geom_msa_track(name = "contacts", kind = "arc",
                   arcs = data.frame(start = c(1, 3), end = c(5, 7)))
  expect_length(w$x$props$columnTracks[[1]]$arcs, 2)
})

test_that("track layers accumulate", {
  w <- msaview(msa = msa) +
    geom_msa_track(c(1), name = "one") +
    geom_msa_track(c(2), name = "two")
  expect_equal(
    vapply(w$x$props$columnTracks, function(t) t$name, character(1)),
    c("one", "two")
  )
})

test_that("domains land as gff text, and a URL as a filehandle", {
  gff <- "##gff-version 3\nHomo_sapiens\tsrc\tDomain\t1\t4\t.\t.\t.\tName=IPR1"
  w <- msaview(msa = msa) + geom_msa_domains(gff)
  expect_match(w$x$props$gff, "IPR1", fixed = TRUE)

  remote <- msaview(msa = msa) +
    geom_msa_domains("https://example.com/domains.gff")
  expect_equal(remote$x$props$gffFilehandle$uri,
               "https://example.com/domains.gff")
  expect_false("gff" %in% names(remote$x$props))
})

test_that("the residue scale sets the scheme and the channel", {
  w <- msaview(msa = msa) + scale_residue_color("clustal", encoding = "color")
  expect_equal(w$x$props$colorScheme, "clustal")
  expect_equal(w$x$props$residueEncoding, "color")
})

test_that("the residue scale takes a color per letter", {
  w <- msaview(msa = msa) + scale_residue_color(c(K = "#1f77b4", D = "#d62728"))
  expect_equal(
    as.character(jsonlite::toJSON(w$x$props$colorScheme, auto_unbox = TRUE)),
    '{"domain":["K","D"],"range":["#1f77b4","#d62728"]}'
  )
})

test_that("an invalid channel is an error from either surface", {
  expect_error(scale_residue_color(encoding = "background"),
               "residue_encoding must be 'fill' or 'color'")
  expect_error(msaview(msa = msa, residue_encoding = "background"),
               "residue_encoding must be 'fill' or 'color'")
})

test_that("theme_msa sets the chrome under its camelCase names", {
  w <- msaview(msa = msa) +
    theme_msa("dark", col_width = 14, draw_tree = FALSE, hide_header = TRUE)
  expect_equal(w$x$props$theme, "dark")
  expect_equal(w$x$props$colWidth, 14)
  expect_false(w$x$props$drawTree)
  expect_true(w$x$props$hideHeader)
})

test_that("theme_msa leaves out what it is not given", {
  w <- msaview(msa = msa) + theme_msa(col_width = 14)
  for (name in c("theme", "rowHeight", "drawTree", "hideHeader")) {
    expect_false(name %in% names(w$x$props))
  }
})

test_that("a later layer replaces an earlier setting", {
  w <- msaview(msa = msa) + theme_msa(col_width = 9) + theme_msa(col_width = 14)
  expect_equal(w$x$props$colWidth, 14)
})

test_that("coord_msa builds the region, sanitizing its row", {
  w <- msaview(msa = msa) + coord_msa(2, 6, row = "Homo sapiens")
  expect_equal(w$x$props$region, list(start = 2L, end = 6L, row = "Homo_sapiens"))
})

test_that("coord_msa without a row spans alignment columns", {
  w <- msaview(msa = msa) + coord_msa(2, 6)
  expect_equal(w$x$props$region, list(start = 2L, end = 6L))
})

test_that("adding something that is not a layer is an error", {
  expect_error(msaview(msa = msa) + 1, "must be an ggmsaview layer")
  expect_error(msaview(msa = msa) + list(a = 1), "must be an ggmsaview layer")
})

test_that("a composed viewer is still an htmlwidget", {
  w <- msaview(msa = msa) + theme_msa("dark")
  expect_s3_class(w, "htmlwidget")
  expect_s3_class(w, "msaview")
})

test_that("a layer prints the props it carries", {
  expect_output(print(theme_msa("dark")), "ggmsaview layer")
  expect_output(print(geom_msa_highlight(start = 1, end = 2)), "highlights")
})

test_that("coord_tree sets the order and roots on an outgroup", {
  w <- msaview(msa = msa) +
    coord_tree(order = "ladderize", root = c("Homo sapiens"))
  expect_equal(w$x$props$treeOrder, "ladderize")
  expect_equal(as.character(w$x$props$treeRoot$outgroup), "Homo_sapiens")
})

test_that("coord_tree passes the midpoint root as a string", {
  w <- msaview(msa = msa) + coord_tree(root = "midpoint")
  expect_equal(w$x$props$treeRoot, "midpoint")
})

test_that("a tree order outside the four is refused", {
  expect_error(coord_tree(order = "sideways"), "tree_order")
  expect_error(coord_tree(order = c("input", "ladderize")), "tree_order")
})

test_that("a rotate clade reaches the props", {
  w <- msaview(msa = msa) +
    geom_msa_clade(c("Homo sapiens", "Mouse"), tips = 2, mark = "rotate")
  expect_equal(w$x$props$clades[[1]]$mark, "rotate")
})

test_that("a features track takes its spans from a data frame", {
  spans <- data.frame(
    start = c(1, 5), end = c(4, 7), strand = c("+", "."),
    type = c("helix", "strand"), color = c(NA, "#9467bd")
  )
  w <- msaview(msa = msa) +
    geom_msa_track(spans, name = "Secondary structure", kind = "features",
                   row = "Homo sapiens", fill = "type", label = "type",
                   position = "strandpile", range = "set1", height = 30)
  track <- w$x$props$columnTracks[[1]]
  expect_equal(track$kind, "features")
  expect_equal(track$row, "Homo_sapiens")
  expect_equal(track$position, "strandpile")
  expect_equal(track$height, 30)
  expect_equal(track$encoding$fill$field, "type")
  expect_equal(track$encoding$fill$scale$scheme, "set1")
  expect_equal(track$encoding$text$field, "type")
  json <- as.character(props_json(w$x$props))
  expect_match(
    json,
    paste0('"features":[{"start":1,"end":4,"strand":"+","type":"helix"},',
           '{"start":5,"end":7,"type":"strand","color":"#9467bd"}]'),
    fixed = TRUE
  )
})

test_that("a features track with no fill or label carries no encoding", {
  w <- msaview(msa = msa) +
    geom_msa_track(data.frame(start = 1, end = 2), name = "spans",
                   kind = "features", color = "#cccccc")
  track <- w$x$props$columnTracks[[1]]
  expect_null(track$encoding)
  expect_equal(track$color, "#cccccc")
  expect_match(as.character(props_json(w$x$props)),
               '"features":[{"start":1,"end":2}]', fixed = TRUE)
})

test_that("a features track's domain_max ends its fill ramp, a bar's its y", {
  features <- msaview(msa = msa) +
    geom_msa_track(data.frame(start = 1, end = 2, score = 3), name = "f",
                   kind = "features", fill = "score", type = "linear",
                   domain_max = 5)
  fill <- features$x$props$columnTracks[[1]]$encoding$fill
  expect_equal(fill$scale$domainMax, 5)
  expect_null(features$x$props$columnTracks[[1]]$encoding$y)

  bar <- msaview(msa = msa) + geom_msa_track(c(1, 2), name = "b", domain_max = 5)
  expect_equal(bar$x$props$columnTracks[[1]]$encoding$y$scale$domainMax, 5)
})

test_that("a bar track's scale arguments color its bars, its ends set its height", {
  w <- msaview(msa = msa) +
    geom_msa_track(c(0.4, 1.8), name = "dN/dS", domain_min = 0.1,
                   domain_max = 2, domain_quantile = 0.9, log = TRUE,
                   type = "threshold", domain = 1,
                   range = c("#2c7fb8", "#d7301f"), title = "ratio")
  encoding <- w$x$props$columnTracks[[1]]$encoding
  expect_equal(encoding$y$scale,
               list(type = "log", domainMin = 0.1, domainMax = 2,
                    domainQuantile = 0.9))
  expect_null(encoding$color$field)
  expect_equal(encoding$color$legend$title, "ratio")
  expect_match(
    as.character(jsonlite::toJSON(encoding$color, auto_unbox = TRUE)),
    '"scale":{"type":"threshold","domain":[1],"range":["#2c7fb8","#d7301f"]}',
    fixed = TRUE
  )
})

test_that("a bar track's color alone is the one color, with no encoding", {
  w <- msaview(msa = msa) +
    geom_msa_track(c(1, 2), name = "score", color = "#c0392b")
  track <- w$x$props$columnTracks[[1]]
  expect_equal(track$color, "#c0392b")
  expect_null(track$encoding)
})

test_that("a strip track's ends set its color ramp", {
  w <- msaview(msa = msa) +
    geom_msa_track(c(30, 90), name = "pLDDT", kind = "strip",
                   type = "linear", range = "viridis", domain_min = 0,
                   domain_max = 100)
  track <- w$x$props$columnTracks[[1]]
  expect_equal(track$kind, "strip")
  expect_null(track$encoding$y)
  expect_equal(track$encoding$color$scale,
               list(type = "linear", scheme = "viridis", domainMin = 0,
                    domainMax = 100))
})

test_that("log is a bar track's alone", {
  expect_error(geom_msa_track(c(1, 2), name = "s", kind = "strip", log = TRUE),
               "bar track's height")
})

test_that("an arc track's scale reads the arcs' value column", {
  w <- msaview(msa = msa) +
    geom_msa_track(data.frame(start = 1, end = 5, value = 0.7),
                   name = "contacts", kind = "arc", type = "linear")
  track <- w$x$props$columnTracks[[1]]
  expect_equal(track$arcs[[1]]$value, 0.7)
  expect_equal(track$encoding$color$scale$type, "linear")
})

test_that("a features track needs spans with a start and an end", {
  expect_error(geom_msa_track(name = "f", kind = "features"), "needs 'features'")
  expect_error(
    geom_msa_track(data.frame(from = 1), name = "f", kind = "features"),
    "needs 'start' and 'end'"
  )
  expect_error(
    geom_msa_track(data.frame(start = 1, end = 2), name = "f",
                   kind = "features", position = "stack"),
    "identity"
  )
})

test_that("a continuous scale writes its ends, its ramp and its legend", {
  w <- msaview(msa = msa) +
    geom_msa_strip("year", type = "linear", range = "viridis",
                   domain_min = 2000, domain_max = 2024, domain_mid = 2012,
                   domain_quantile = 0.99, reverse = TRUE,
                   missing_label = "no date")
  json <- as.character(props_json(w$x$props$rowPanels[[1]]$encoding$fill))
  expect_equal(json, paste0(
    '{"field":"year","scale":{"type":"linear","scheme":"viridis",',
    '"domainMin":2000,"domainMax":2024,"domainMid":2012,',
    '"domainQuantile":0.99,"reverse":true},',
    '"legend":{"missingLabel":"no date"}}'
  ))
})

test_that("a threshold scale keeps its cuts numbers and its labels an array", {
  w <- msaview(msa = msa) +
    scale_row_color("plddt", channel = "featureFill", type = "threshold",
                    domain = c(50, 70, 90),
                    range = c("#ff7d45", "#ffdb13", "#65cbf3", "#0053d6"),
                    labels = c("Very low", "Low", "Confident", "Very high"))
  json <- as.character(props_json(w$x$props$encodings[[1]]))
  expect_match(json, '"domain":[50,70,90]', fixed = TRUE)
  expect_match(json, '"labels":["Very low","Low","Confident","Very high"]',
               fixed = TRUE)
})

test_that("a label array of one stays an array", {
  w <- msaview(msa = msa) +
    scale_row_color("clade", domain = "a", range = "#e41a1c", labels = "A")
  expect_match(as.character(props_json(w$x$props$encodings[[1]])),
               '"legend":{"labels":["A"]}', fixed = TRUE)
})

test_that("an unknown scale type is refused", {
  expect_error(scale_row_color("year", type = "continuous"), "'linear'")
  expect_error(geom_msa_strip("year", type = "sqrt"), "'threshold'")
})

test_that("a features panel takes a position", {
  w <- msaview(msa = msa) + geom_msa_features(position = "strandpile")
  expect_equal(w$x$props$rowPanels[[1]]$position, "strandpile")
  expect_error(geom_msa_features(position = "dodge"), "strandpile")
})

test_that("a selection trims the view only when asked", {
  plain <- msaview(msa = msa) + geom_msa_selection(2, 5)
  expect_false("trimToSelection" %in% names(plain$x$props))
  trimmed <- msaview(msa = msa) + geom_msa_selection(2, 5, trim = TRUE)
  expect_true(trimmed$x$props$trimToSelection)
  expect_true(msaview(msa = msa, trim_to_selection = TRUE)$x$props$trimToSelection)
})

test_that("the row layers merge into one row transform", {
  w <- msaview(msa = msa) +
    filter_msa_rows("host", values = c("duck", "goose")) +
    filter_msa_rows("year", range = c(2010, NA)) +
    arrange_msa_rows("host", desc("year")) +
    facet_msa_rows("host", sort = "wild")
  expect_equal(
    as.character(props_json(w$x$props$rowTransform)),
    paste0(
      '{"filter":[{"field":"host","oneOf":["duck","goose"]},',
      '{"field":"year","range":[2010,null]}],',
      '"sort":[{"field":"host"},{"field":"year","order":"descending"}],',
      '"facet":{"field":"host","sort":["wild"]}}'
    )
  )
})

test_that("a later sort or facet replaces an earlier one", {
  w <- msaview(msa = msa) +
    arrange_msa_rows("host") + arrange_msa_rows("year") +
    facet_msa_rows("host") + facet_msa_rows("clade")
  expect_length(w$x$props$rowTransform$sort, 1)
  expect_equal(w$x$props$rowTransform$sort[[1]]$field, "year")
  expect_equal(w$x$props$rowTransform$facet$field, "clade")
})

test_that("the row layers add to a row transform msaview() was given", {
  w <- msaview(msa = msa, row_transform = list(
    filter = list(list(field = "host", oneOf = "duck"))
  )) + filter_msa_rows("year", range = c(NA, 2020))
  filters <- w$x$props$rowTransform$filter
  expect_length(filters, 2)
  expect_match(as.character(props_json(filters)),
               '[{"field":"host","oneOf":["duck"]},{"field":"year","range":[null,2020]}]',
               fixed = TRUE)
})

test_that("arrange_msa_rows reads desc(), variables and order", {
  field <- "year"
  keys <- function(layer) layer$into$rowTransform$set$sort
  expect_equal(keys(arrange_msa_rows(desc(field))),
               list(list(field = "year", order = "descending")))
  # dplyr's desc, written as the call reads it, which R CMD check would
  # otherwise take for an undeclared dependency
  expect_equal(keys(eval(str2lang('arrange_msa_rows(dplyr::desc("year"))'))),
               list(list(field = "year", order = "descending")))
  expect_equal(keys(arrange_msa_rows(c("host", field))),
               list(list(field = "host"), list(field = "year")))
  expect_equal(
    keys(arrange_msa_rows("host", "year", order = c("descending", "ascending"))),
    list(list(field = "host", order = "descending"),
         list(field = "year", order = "ascending"))
  )
  expect_equal(keys(arrange_msa_rows(desc("year"), order = "ascending")),
               list(list(field = "year", order = "ascending")))
})

test_that("the row layers refuse what the viewer would not read", {
  expect_error(arrange_msa_rows(), "needs a field")
  expect_error(arrange_msa_rows(2010), "as strings")
  expect_error(arrange_msa_rows("a", "b", order = c("ascending", "x", "y")),
               "one per field")
  expect_error(arrange_msa_rows("a", order = "down"), "'descending'")
  expect_error(filter_msa_rows("year"), "'range'")
  expect_error(filter_msa_rows("year", values = 1, range = c(1, 2)), "'range'")
  expect_error(filter_msa_rows("year", range = 2010), "c\\(min, max\\)")
  expect_error(msaview(msa = msa, row_transform = list(group = "x")), "facet")
})

test_that("a row transform layer prints the prop it fills", {
  expect_output(print(facet_msa_rows("host")), "rowTransform")
})
