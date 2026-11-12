skip_if_not_installed("QuickJSR")
skip_if_not(has_bundle("js/msaview-core.js"), "scripts/build-js.mjs has not run")

msa <- ">human\nMKAANSE\n>mouse\nMKA-NSE\n>fish\nMRA-NSQ"
tree <- "((fish:1,mouse:1):1,human:2);"

test_that("rows come out in the tree's order with their visible columns", {
  lay <- msaview_layout(msaview(msa = msa, tree = tree))
  expect_equal(lay$rows$name, c("fish", "mouse", "human"))
  expect_equal(lay$rows$seq[3], "MKAANSE")
  expect_equal(lay$num_columns, 7)
  expect_true(lay$tree$by_length)
  expect_equal(sort(unique(c(lay$tree$segments$y, lay$tree$segments$yend))),
               c(1, 1.5, 2, 3))
})

test_that("a static scheme colors by letter, a dynamic one by column", {
  static <- msaview_layout(msaview(msa = msa, color_scheme = "clustal"))
  expect_true(all(is.na(static$colors$col)))
  expect_match(static$colors$fill, "^#[0-9a-f]{6,8}$")
  dynamic <- msaview_layout(
    msaview(msa = msa, color_scheme = "clustalx_protein_dynamic")
  )
  expect_false(anyNA(dynamic$colors$col))
})

test_that("a collapsed clade and hidden gaps reach the tables", {
  lay <- msaview_layout(
    msaview(msa = ">a\nMK-A\n>b\nMK-A\n>c\nMKWA", tree = "((a,b),c);",
            allowed_gappyness = 50)
  )
  expect_equal(lay$num_columns, 3)
})

test_that("domains and highlights land in visible columns", {
  gff <- "##gff-version 3\nhuman\tpfam\tdomain\t2\t4\t.\t.\t.\tName=Dom"
  w <- msaview(msa = msa, gff = gff) +
    geom_msa_highlight(row = "mouse", start = 4, end = 5, label = "site")
  lay <- msaview_layout(w)
  expect_equal(lay$domains$row, 1)
  expect_equal(c(lay$domains$start, lay$domains$end), c(2, 4))
  expect_equal(lay$domains$label, "Dom")
  expect_equal(c(lay$highlights$start, lay$highlights$end), c(5, 6))
  expect_equal(lay$highlights$label, "site")
})

test_that("R reads a URL the props name, since QuickJS fetches nothing", {
  fasta <- tempfile(fileext = ".fa")
  writeLines(msa, fasta)
  props <- inline_remote(list(
    msaUrl = paste0("file://", fasta),
    treeFilehandle = list(uri = paste0("file://", fasta))
  ))
  expect_equal(props$msa, msa)
  expect_null(props$msaUrl)
  expect_null(props$treeFilehandle)
})

test_that("highlight_columns draws as 1-based column highlights", {
  lay <- msaview_layout(msaview(msa = msa, highlight_columns = c(2, 3, 6)))
  expect_equal(lay$highlights$start, c(2, 6))
  expect_equal(lay$highlights$end, c(3, 6))
})

test_that("the plot builds, with the tree beside the alignment", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa, tree = tree, gff =
    "##gff-version 3\nhuman\tpfam\tdomain\t2\t4\t.\t.\t.\tName=Dom")
  p <- msaview_plot(w, tracks = "sequence-logo")
  expect_s3_class(p, "patchwork")
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, p, width = 6, height = 2)
  expect_gt(file.size(file), 1000)

  alone <- msaview_plot(msaview(msa = msa), tracks = FALSE)
  expect_s3_class(alone, "ggplot")
  built <- ggplot2::ggplot_build(alone)
  expect_true(any(vapply(built$data, function(d) "label" %in% names(d),
                         logical(1))))
})

test_that("bracket highlights draw a strip and outline their columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa) +
    geom_msa_highlight(start = 2, end = 6, label = "outer", mark = "bracket") +
    geom_msa_highlight(start = 3, end = 3, label = "inner", mark = "bracket") +
    geom_msa_highlight(start = 7, end = 7, label = "band")
  lay <- msaview_layout(w)
  expect_equal(lay$highlights$mark, c("bracket", "bracket", "band"))
  strip <- Filter(function(t) t$kind == "highlights", lay$tracks)[[1]]
  expect_equal(strip$id, "highlight-brackets")
  expect_equal(strip$bands$label, c("outer", "inner"))
  expect_equal(strip$bands$lane, c(0, 1))

  built <- ggplot2::ggplot_build(highlight_track(strip))
  expect_true(any(vapply(built$data, function(d)
    all(c("outer", "inner") %in% d$label), logical(1))))

  with_strip <- ggplot2::ggplot_build(alignment_plot(lay, FALSE, 2, TRUE))
  rects <- Filter(function(d) all(c("xmin", "fill") %in% names(d)),
                  with_strip$data)
  outlines <- do.call(rbind, lapply(rects, function(d) d[is.na(d$fill), ]))
  expect_equal(sort(outlines$xmin), c(1.5, 2.5))

  p <- msaview_plot(w, tracks = c("conservation", "highlight-brackets"))
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, p, width = 6, height = 3)
  expect_gt(file.size(file), 1000)
})

test_that("the layout lists every track, flagging the ones the viewer shows", {
  lay <- msaview_layout(msaview(msa = msa))
  ids <- vapply(lay$tracks, function(t) t$id, character(1))
  on <- vapply(lay$tracks, function(t) t$on, logical(1))
  expect_equal(ids, c("conservation", "property-conservation", "sequence-logo"))
  expect_equal(on, c(TRUE, TRUE, FALSE))
  expect_length(lay$tracks[[1]]$values, 7)
  logo <- lay$tracks[[3]]$letters
  expect_true(all(logo$ymax > logo$ymin))
  expect_lte(max(logo$ymax), lay$tracks[[3]]$max_bits)
})

test_that("a column track of your own reaches the figure", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa) +
    geom_msa_track(c(0.1, 0.5, 1, 0, 0, 0.2, 0.9), name = "score") +
    geom_msa_track(name = "pairs", kind = "arc",
                   arcs = data.frame(start = 1, end = 5))
  lay <- msaview_layout(w)
  ids <- vapply(lay$tracks, function(t) t$id, character(1))
  expect_true(all(c("score", "pairs") %in% ids))
  arcs <- lay$tracks[[which(ids == "pairs")]]$arcs
  expect_equal(c(arcs$start, arcs$end), c(1, 5))
  p <- msaview_plot(w, tracks = c("score", "pairs", "sequence-logo"))
  # three track rows over the alignment, beside the row names; patchwork holds
  # the last panel outside its patches
  expect_length(p$patches$plots, 7)
  expect_error(msaview_plot(w, tracks = "nope"), "no track named 'nope'")
})

test_that("features and row panels come out as the mark the viewer draws", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  gff <- paste(
    "##gff-version 3",
    "human\tdemo\tgene\t1\t4\t.\t+\t.\tName=genA",
    "mouse\tdemo\tgene\t2\t6\t.\t.\t.\tName=genB",
    sep = "\n"
  )
  w <- msaview(msa = msa, gff = gff) +
    geom_msa_rowdata(data.frame(name = c("human", "mouse", "fish"),
                                clade = c("a", "a", "b"))) +
    geom_msa_strip("clade") +
    geom_msa_features(x = "column", label = "Name")
  lay <- msaview_layout(w)
  # an arrow where the feature has a strand, a box where it has none
  expect_equal(sort(as.vector(table(lay$features$shapes$group))), c(4, 5))
  expect_false(lay$fill_cells)
  expect_equal(vapply(lay$row_panels, function(p) p$kind, ""),
               c("strip", "features"))
  expect_equal(lay$row_panels[[1]]$cells$row, 1:3)
  expect_equal(lay$row_panels[[1]]$title, "clade")
  expect_equal(lay$row_panels[[1]]$cells$key, c("a", "a", "b"))
  expect_setequal(unique(lay$features$shapes$key), c("genA", "genB"))
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, msaview_plot(w, tracks = FALSE), width = 6, height = 2)
  expect_gt(file.size(file), 1000)
})

test_that("a bracket clade and the tree's scale bar reach the figure", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa, tree = tree) +
    geom_msa_clade(mrca = c("fish", "mouse"), tips = 2, mark = "bracket",
                   label = "pair", color = "#1565c0")
  lay <- msaview_layout(w)
  gutter <- lay$clade_gutter
  expect_equal(c(gutter$bars$first, gutter$bars$last), c(1, 2))
  expect_equal(gutter$labels$label, "pair")
  expect_gt(gutter$width, 0)
  expect_false(is.null(lay$tree$scale_bar$step))
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, msaview_plot(w, tracks = FALSE), width = 6, height = 2)
  expect_gt(file.size(file), 1000)
})

test_that("a row encoding colors the names, the rows and the branches", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa, tree = tree) +
    geom_msa_rowdata(data.frame(name = c("human", "mouse", "fish"),
                                group = c("land", "land", "sea"))) +
    scale_row_color("group", channel = "tipLabel") +
    scale_row_color("group", channel = "rowTint") +
    scale_row_color("group", channel = "branch")
  lay <- msaview_layout(w)
  expect_false(anyNA(lay$rows$label_color))
  expect_false(anyNA(lay$rows$tint))
  expect_gt(sum(!is.na(lay$tree$segments$color)), 0)
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, msaview_plot(w, tracks = FALSE), width = 6, height = 2)
  expect_gt(file.size(file), 1000)
})

test_that("an NA in a layer's data frame takes the default", {
  lay <- msaview_layout(msaview(msa = msa) + geom_msa_highlight(
    data.frame(start = c(1, 3), end = c(2, 4), color = c(NA, "#123456"))))
  expect_equal(lay$highlights$color, c(NA, "#123456"))
})

test_that("features = FALSE leaves the overlay out and fills the cells", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  gff <- "##gff-version 3\nhuman\tpfam\tdomain\t2\t4\t.\t.\t.\tName=Dom"
  w <- msaview(msa = msa, gff = gff)
  expect_false(msaview_layout(w)$fill_cells)
  bare <- without_overlay(msaview_layout(w))
  expect_null(bare$features)
  expect_true(bare$fill_cells)
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, msaview_plot(w, features = FALSE), width = 6, height = 2)
  expect_gt(file.size(file), 1000)
})

years <- data.frame(name = c("human", "mouse", "fish"),
                    year = c(2001, 2010, 2020), host = c("land", "land", "sea"))

test_that("a ramp's legend comes out as its stops and its two ends", {
  w <- msaview(msa = msa) + geom_msa_rowdata(years) +
    geom_msa_strip("year", type = "linear", range = "viridis",
                   domain_max = 2015)
  legend <- msaview_layout(w)$legends[[1]]
  expect_equal(legend$title, "year")
  expect_length(legend$value, 0)
  g <- legend$gradient
  expect_length(g$offset, 11)
  expect_equal(g$color[1], "#440154")
  expect_equal(c(g$min, g$max), c(2001, 2015))
  expect_equal(g$max_label, "≥2015")
  expect_false(g$log)
})

test_that("a ramp strip draws through a gradient scale over its ends", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa) + geom_msa_rowdata(years) +
    geom_msa_strip("year", type = "log", domain_min = 1000, domain_max = 3000)
  lay <- msaview_layout(w)
  expect_true(lay$legends[[1]]$gradient$log)
  p <- row_panel_plot(lay$row_panels[[1]], lay, 3)
  scale <- p$scales$get_scales("fill")
  expect_s3_class(scale, "ScaleContinuous")
  expect_equal(scale$get_limits(), log10(c(1000, 3000)))
  expect_equal(scale$get_labels(), c("1000", "3000"))
  built <- ggplot2::ggplot_build(p)
  tiles <- Filter(function(d) "fill" %in% names(d), built$data)[[1]]
  # ggplot2 interpolates between the 11 stops in Lab, the viewer along the
  # ramp's own table, so the two agree to a few units per channel
  drift <- abs(grDevices::col2rgb(tiles$fill) -
                 grDevices::col2rgb(lay$row_panels[[1]]$cells$color))
  expect_lte(max(drift), 12)
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, msaview_plot(w, tracks = FALSE), width = 6, height = 2)
  expect_gt(file.size(file), 1000)
})

test_that("a threshold strip keys each interval under its own label", {
  skip_if_not_installed("ggplot2")
  w <- msaview(msa = msa) + geom_msa_rowdata(years) +
    geom_msa_strip("year", type = "threshold", domain = c(2005, 2015),
                   labels = c("early", "mid", "late"))
  lay <- msaview_layout(w)
  legend <- lay$legends[[1]]
  expect_equal(legend$label, c("early", "mid", "late"))
  expect_equal(legend$value, c("< 2005", "2005 – 2015", "≥ 2015"))
  p <- row_panel_plot(lay$row_panels[[1]], lay, 3)
  scale <- p$scales$get_scales("fill")
  expect_equal(as.character(scale$get_breaks()), legend$value)
  expect_equal(unname(scale$get_labels()), legend$label)
})

test_that("a legend's note reaches the layout and prints when plotting", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa) + geom_msa_rowdata(years) +
    geom_msa_strip("year", type = "linear", domain_min = 2020,
                   domain_max = 2000)
  lay <- msaview_layout(w)
  expect_equal(lay$scale_problems$record, "rowPanels[0].encoding.fill")
  expect_equal(lay$scale_problems$slot, "scale.domainMax")
  expect_match(lay$legends[[1]]$note, "domainMax is below domainMin")
  expect_message(msaview_plot(w), "year: domainMax is below domainMin")
})

test_that("a features track's spans reach the figure as the span mark", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  spans <- data.frame(start = c(1, 3, 6), end = c(4, 5, 7),
                      strand = c("+", NA, "-"),
                      type = c("helix", "strand", "helix"))
  w <- msaview(msa = msa) +
    geom_msa_track(spans, name = "SS", kind = "features", fill = "type",
                   label = "type", range = c(helix = "#e15759",
                                             strand = "#edc948"))
  lay <- msaview_layout(w)
  track <- Filter(function(t) t$kind == "features", lay$tracks)[[1]]
  expect_true(track$on)
  expect_equal(track$title, "type")
  # two strand arrows, five points each, and one box
  expect_equal(sort(as.vector(table(track$shapes$group))), c(4, 5, 5))
  expect_setequal(unique(track$shapes$key), c("helix", "strand"))
  expect_true(all(track$shapes$y >= 0 & track$shapes$y <= track$height))
  expect_equal(lay$legends[[1]]$value, c("helix", "strand"))

  built <- ggplot2::ggplot_build(track_plot(track, lay))
  polygons <- Filter(function(d) "group" %in% names(d) && "fill" %in% names(d),
                     built$data)[[1]]
  expect_setequal(unique(polygons$fill), c("#e15759", "#edc948"))
  p <- msaview_plot(w, tracks = "ss")
  file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(file, p, width = 6, height = 2)
  expect_gt(file.size(file), 1000)
})

test_that("a trimmed view draws the selected columns under their own numbers", {
  skip_if_not_installed("ggplot2")
  lay <- msaview_layout(msaview(msa = msa) + geom_msa_selection(3, 6, trim = TRUE))
  expect_equal(lay$num_columns, 4)
  expect_equal(lay$file_columns, 3:6)
  expect_equal(lay$rows$seq[1], "AANS")
  ruler <- column_ruler(lay$file_columns)
  expect_equal(ruler$labels, lay$file_columns[ruler$breaks])
  expect_equal(column_ruler(1:100)$labels, c(20, 40, 60, 80))
})

test_that("facet bands draw as brackets in the clade gutter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  w <- msaview(msa = msa) + geom_msa_rowdata(years) +
    facet_msa_rows("host", sort = "sea")
  lay <- msaview_layout(w)
  expect_equal(lay$rows$name, c("fish", "human", "mouse"))
  gutter <- lay$clade_gutter
  expect_equal(gutter$labels$label, c("sea", "land"))
  expect_equal(c(gutter$bars$first, gutter$bars$last), c(1, 2, 1, 3))
  built <- ggplot2::ggplot_build(gutter_plot(gutter, lay, 3))
  expect_true(any(vapply(built$data, function(d)
    all(c("sea", "land") %in% d$label), logical(1))))
})

test_that("a filter and a sort reorder the rows the figure draws", {
  w <- msaview(msa = msa) + geom_msa_rowdata(years) +
    filter_msa_rows("year", range = c(2005, NA)) +
    arrange_msa_rows(desc("year"))
  expect_equal(msaview_layout(w)$rows$name, c("fish", "mouse"))
})
