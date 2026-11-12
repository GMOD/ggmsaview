the <- new.env(parent = emptyenv())

utils::globalVariables(".data")

need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(what, " needs the ", pkg, " package: install.packages(\"", pkg, "\")",
         call. = FALSE)
  }
}

js_context <- function() {
  if (is.null(the$ctx)) {
    need("QuickJSR", "A static figure")
    if (!has_bundle("js/msaview-core.js")) stop(bundle_hint, call. = FALSE)
    ctx <- QuickJSR::JSContext$new(stack_size = 64 * 1024 * 1024)
    ctx$source(file = system.file("js", "msaview-core.js", package = "ggmsaview"))
    the$ctx <- ctx
  }
  the$ctx
}

as_frame <- function(x) {
  if (length(x) == 0 || length(x[[1]]) == 0) {
    return(NULL)
  }
  as.data.frame(lapply(x, function(col) {
    values <- unlist(col)
    type <- if (is.character(values)) character(1) else
      if (is.logical(values) || is.null(values)) logical(1) else numeric(1)
    na <- type[NA_integer_]
    vapply(col, function(v) if (is.null(v)) na else v, type)
  }), stringsAsFactors = FALSE)
}

#' Lay out a viewer as tables
#'
#' Runs the viewer's own model, in an embedded QuickJS engine, over the
#' widget's props and returns what it would draw as data frames: the rows in
#' tree order, the color of each letter, the tree's segments, the domains,
#' highlights and clades in visible columns, the legends and every column
#' track the viewer can draw above the alignment. Collapsed clades, a focus,
#' hidden gappy columns and the color scheme are all applied, as on screen.
#' [msaview_plot()] draws the result.
#'
#' @param widget A viewer from [msaview()], with any layers added.
#' @return An `msaview_layout`: a list with `rows`, `colors`, `tree`,
#'   `domains`, `features`, `row_panels`, `highlights`, `clades`, `legends`,
#'   `scale_problems`, `tracks`, `clade_gutter`, `num_columns`,
#'   `file_columns`, `row_height`, `col_width`, `font_size`, `reference` and
#'   `fill_cells`.
#'
#'   `features` and each features panel hold `shapes`, one polygon per span
#'   (an arrow where the feature has a strand) in pixels across and rows down,
#'   with the `key` its fill scale reads; the `labels` that fit; and the
#'   legend `title`. A bar or point panel holds `marks`, each row's mark in
#'   the panel's pixels from `base` to `x`, a text panel `texts`, the text
#'   each row prints, and both the `key` and `color` of each mark a color
#'   field reads. A panel with an axis, a bar, a point or a features panel
#'   over positions, holds its ticks in `axis`, their pixel positions `at`
#'   and their `label`s. `rows` holds each row's `name` and the `label` its
#'   tip prints.
#'
#'   Each track is a list with its `id`, `name`, `kind`, pixel `height`,
#'   whether the viewer shows it (`on`), and its content: `bars` for a bar,
#'   each a `value` from 0 to 1 of the track height, `letters` for a logo,
#'   `arcs` for an arc, `cells` for a strip or a text track, and `shapes`,
#'   `labels` and `title` for a features track, its shapes in pixels down
#'   from the track's top. A bar, cell or arc a scale over the track's own
#'   values colors carries the value it read as `key`, which the legend the
#'   track's `title` names lists. The tree's `nodes` are the internal nodes
#'   the `nodePoint` channel colors, keyed by their support.
#'
#'   Each legend is a list with its `title`, a `note` the viewer prints under
#'   it, its swatches as the parallel vectors `value`, `label` and `color`,
#'   and for a linear or log scale a `gradient`: the ramp's stops (`offset`
#'   and `color`), its ends as numbers (`min`, `max`) and as the legend prints
#'   them (`min_label`, `max_label`), and whether the scale is `log`.
#'   `scale_problems` lists each contradiction in a color object as its
#'   `record`, `rule`, `slot` and `message`. `file_columns` gives the file's
#'   column number of each drawn column, which hidden gappy columns and
#'   `trim_to_selection` make differ from its position.
#' @importFrom utils globalVariables
#' @export
#' @examples
#' \dontrun{
#' msaview_layout(msaview(msa = ">a\nMKAANSE\n>b\nMKA-NSE"))
#' }
msaview_layout <- function(widget) {
  if (!inherits(widget, "msaview")) {
    stop("widget must be a viewer made by msaview()")
  }
  props <- inline_remote(widget$x$props)
  out <- js_context()$call("msaviewFigure", as.character(props_json(props)))
  raw <- jsonlite::fromJSON(out, simplifyVector = FALSE)
  tree <- raw$tree
  structure(
    list(
      rows = data.frame(
        name = vapply(raw$rows$name, as.character, character(1)),
        label = vapply(raw$rows$label, as.character, character(1)),
        seq = vapply(raw$rows$seq, as.character, character(1)),
        label_color = vapply(raw$rows$labelColor, function(v)
          if (is.null(v)) NA_character_ else v, character(1)),
        tint = vapply(raw$rows$tint, function(v)
          if (is.null(v)) NA_character_ else v, character(1)),
        stringsAsFactors = FALSE
      ),
      colors = as_frame(raw$colors),
      tree = if (!is.null(tree)) list(
        segments = as_frame(tree$segments),
        collapsed = as_frame(tree$collapsed),
        nodes = as_frame(tree$nodes),
        node_title = tree$nodeTitle,
        by_length = isTRUE(tree$byLength),
        scale_bar = tree$scaleBar
      ),
      clade_gutter = if (!is.null(raw$cladeGutter)) list(
        width = raw$cladeGutter$width,
        bar_x = raw$cladeGutter$barX,
        bars = as_frame(raw$cladeGutter$bars),
        labels = as_frame(raw$cladeGutter$labels)
      ),
      domains = as_frame(raw$domains),
      highlights = as_frame(raw$highlights),
      clades = as_frame(raw$clades),
      legends = lapply(raw$legends, parse_legend),
      scale_problems = records_frame(raw$scaleProblems,
                                     c("record", "rule", "slot", "message")),
      features = parse_spans(raw$features),
      row_panels = lapply(raw$rowPanels, parse_row_panel),
      tracks = lapply(raw$tracks, parse_track),
      num_columns = raw$numColumns,
      file_columns = vapply(raw$fileColumns, as.numeric, numeric(1)),
      row_height = raw$rowHeight,
      col_width = raw$colWidth,
      font_size = raw$fontSize,
      reference = raw$reference,
      fill_cells = isTRUE(raw$fillCells),
      letter_mode = isTRUE(raw$letterMode),
      letters_over_domains = isTRUE(raw$lettersOverDomains),
      tree_area_width = raw$treeAreaWidth
    ),
    class = "msaview_layout"
  )
}

parse_spans <- function(x) {
  if (is.null(x)) return(NULL)
  list(title = x$title, shapes = as_frame(x$shapes),
       labels = as_frame(x$labels))
}

parse_row_panel <- function(p) {
  c(
    list(kind = p$kind, header = p$header, width = p$width, title = p$title,
         size = p$size, padding = p$padding, cells = as_frame(p$cells),
         marks = as_frame(p$marks), texts = as_frame(p$texts),
         axis = if (!is.null(p$axis)) list(
           at = vapply(p$axis$at, as.numeric, numeric(1)),
           label = strings(p$axis$label)
         )),
    if (identical(p$kind, "features")) parse_spans(p)
  )
}

parse_track <- function(t) {
  c(
    list(
      id = t$id, name = t$name, kind = t$kind, on = isTRUE(t$on),
      height = t$height, title = t$title, max_bits = t$maxBits,
      bars = as_frame(t$bars),
      letters = as_frame(t$letters),
      arcs = as_frame(t$arcs),
      cells = as_frame(t$cells),
      bands = as_frame(t$bands)
    ),
    if (identical(t$kind, "features")) parse_spans(t)
  )
}

strings <- function(x) vapply(x, as.character, character(1))

parse_legend <- function(l) {
  g <- l$gradient
  list(
    id = l$id, title = l$title, note = l$note,
    value = strings(l$value), label = strings(l$label),
    color = strings(l$color),
    gradient = if (!is.null(g)) list(
      offset = vapply(g$offset, as.numeric, numeric(1)),
      color = strings(g$color),
      min = g$min, max = g$max,
      min_label = g$minLabel, max_label = g$maxLabel,
      log = isTRUE(g$log)
    )
  )
}

# a JSON array of objects as a data frame of the named fields, NULL when empty
records_frame <- function(records, fields) {
  if (length(records) == 0) return(NULL)
  as.data.frame(lapply(stats::setNames(fields, fields), function(field) {
    vapply(records, function(r) as.character(r[[field]] %||% NA), character(1))
  }), stringsAsFactors = FALSE)
}

# QuickJS fetches nothing, so R reads each URL the props name into the text
# prop the viewer would have filled from it
inline_remote <- function(props) {
  for (layer in c("msa", "tree", "gff")) {
    url_prop <- paste0(layer, "Url")
    handle_prop <- paste0(layer, "Filehandle")
    url <- props[[handle_prop]]$uri %||% props[[url_prop]]
    if (!is.null(url) && is.null(props[[layer]])) {
      props[[layer]] <- paste(readLines(url, warn = FALSE), collapse = "\n")
    }
    props[[url_prop]] <- NULL
    props[[handle_prop]] <- NULL
  }
  props
}

#' @export
print.msaview_layout <- function(x, ...) {
  cat("<msaview_layout>", nrow(x$rows), "rows x", x$num_columns, "columns",
      if (!is.null(x$tree)) "with a tree", "\n")
  invisible(x)
}

# one letter per cell, rows top to bottom, as a character matrix
letter_matrix <- function(layout) {
  n <- layout$num_columns
  chars <- strsplit(layout$rows$seq, "", fixed = TRUE)
  m <- matrix(NA_character_, nrow = length(chars), ncol = n)
  for (i in seq_along(chars)) {
    k <- min(length(chars[[i]]), n)
    if (k > 0) m[i, seq_len(k)] <- chars[[i]][seq_len(k)]
  }
  m
}

cell_colors <- function(layout, letters, field) {
  colors <- layout$colors
  out <- matrix(NA_character_, nrow(letters), ncol(letters))
  if (is.null(colors)) return(out)
  if (all(is.na(colors$col))) {
    out[] <- colors[[field]][match(letters, colors$letter)]
  } else {
    cols <- col(letters)
    out[] <- colors[[field]][match(paste(cols, letters),
                                   paste(colors$col, colors$letter))]
  }
  out
}

reference_matches <- function(layout, letters) {
  ref <- match(layout$reference, layout$rows$name)
  if (is.null(layout$reference) || is.na(ref)) {
    return(matrix(FALSE, nrow(letters), ncol(letters)))
  }
  same <- sweep(letters, 2, letters[ref, ], "==")
  same[is.na(same)] <- FALSE
  same[ref, ] <- FALSE
  same
}

#' Draw a viewer as a ggplot
#'
#' Draws an [msaview_layout()] with ggplot2: the alignment as a raster of
#' colored cells with its letters on top, the domains and highlights over it,
#' the row names down its side, and the tree beside it. The figure is ordinary
#' ggplot2, so `ggsave()`, themes and patchwork apply, and knitr renders it in
#' PDF and Word output.
#'
#' Strips, gene arrows, node points and the bars, strip cells and text cells
#' of a column track map `fill` through the scale the viewer resolved, and
#' points, texts and arcs map `colour` through it: a `scale_fill_manual()` of
#' its colors for a palette, a threshold or an identity scale, and a
#' `scale_fill_gradientn()` over its ramp's stops for a linear or log scale.
#' A legend the viewer gives a note, such as a contradiction in a color
#' object's scale, prints the note as a message.
#'
#' The row panels draw as ggtree's `geom_facet()` facets beside the tree: a
#' bar panel as `geom_rect()` bars from 0, a point panel as `geom_point()`
#' and a text panel as `geom_text()`, each over an axis at the panel's top
#' with the ticks and labels the viewer's header band draws. A features panel
#' over positions draws its position axis the same way. The tip labels print
#' the field a `rowLabel` encoding names.
#'
#' @param x A viewer from [msaview()], or its [msaview_layout()].
#' @param letters Draw each residue's letter. By default, when the alignment
#'   has at most 4000 cells and its letters come out at least 1.2 mm.
#' @param letter_size Letter size in mm. By default as large as a column of a
#'   figure about 8 inches wide holds, at most 3.
#' @param row_labels Name each row down the alignment's side. By default, when
#'   there are at most 60 rows.
#' @param tree_width The tree's width as a fraction of the alignment's.
#' @param tracks The column tracks to draw above the alignment: `TRUE` for
#'   the ones the viewer shows, `FALSE` for none, or the ids or names of the
#'   ones to draw, such as `"sequence-logo"`, which the viewer hides by
#'   default. `msaview_layout(x)$tracks` lists them all.
#' @param features Draw the GFF's spans over the alignment. `FALSE` leaves the
#'   overlay out and colors the cells as with no GFF, for a GFF that feeds a
#'   row panel, or whose coordinates are not the alignment's. The panels still
#'   draw it.
#' @return A ggplot, or a patchwork of the tree, the alignment and the tracks
#'   when there is more than one panel.
#' @export
#' @examples
#' \dontrun{
#' w <- msaview(msa = ">human\nMKAANSE\n>mouse\nMKA-NSE",
#'              tree = "(human:0.1,mouse:0.2);") +
#'   geom_msa_highlight(start = 2, end = 4, label = "motif")
#' msaview_plot(w)
#' ggplot2::ggsave("alignment.pdf", msaview_plot(w), width = 6, height = 2)
#' }
msaview_plot <- function(x, letters = NULL, letter_size = NULL, row_labels = NULL,
                         tree_width = 0.3, tracks = TRUE, features = TRUE) {
  need("ggplot2", "msaview_plot()")
  layout <- if (inherits(x, "msaview_layout")) x else msaview_layout(x)
  if (!features) {
    layout <- without_overlay(layout)
  }
  layout <- merge_feature_keys(layout)
  for (legend in layout$legends) {
    if (!is.null(legend$note)) message(legend$title, ": ", legend$note)
  }
  letter_size <- letter_size %||% min(3, 200 / max(layout$num_columns, 1))
  letters <- letters %||% (nrow(layout$rows) * layout$num_columns <= 4000 &&
                             letter_size >= 1.2)
  row_labels <- row_labels %||% (nrow(layout$rows) <= 60)
  shown <- if (layout$num_columns > 0) pick_tracks(layout$tracks, tracks)
  strip <- any(vapply(shown, function(t) identical(t$kind, "highlights"),
                      logical(1)))
  aln <- if (layout$num_columns > 0) {
    alignment_plot(layout, letters, letter_size, strip)
  }
  if (length(shown) == 0 && is.null(layout$tree) && !row_labels &&
      is.null(layout$clade_gutter) && length(layout$row_panels) == 0) {
    return(aln)
  }
  need("patchwork", "A figure with a tree, row names, a row panel or a track")
  # one column per panel left to right, and the track rows above the main one
  columns <- list()
  if (!is.null(layout$tree)) {
    columns$tree <- list(plot = tree_plot(layout),
                         width = grid::unit(tree_width, "null"))
  }
  if (row_labels) {
    columns$labels <- list(
      plot = labels_plot(layout),
      # names are text at a fixed size, so the column is measured in it
      width = grid::unit(0.3 + 0.15 * max(nchar(layout$rows$label), 1), "cm")
    )
  }
  gutter <- layout$clade_gutter
  if (!is.null(gutter)) {
    columns$gutter <- list(
      plot = gutter_plot(gutter, layout, 3),
      # the labels are text at a fixed size, so the column is measured in it
      width = grid::unit(0.6 + 0.2 * max(nchar(gutter$labels$label), 0), "cm")
    )
  }
  for (i in seq_along(layout$row_panels)) {
    panel <- layout$row_panels[[i]]
    columns[[paste0("panel", i)]] <- list(
      plot = row_panel_plot(panel, layout, 3),
      width = grid::unit(tree_width * panel$width / layout$tree_area_width,
                         "null")
    )
  }
  if (!is.null(aln)) {
    columns$aln <- list(plot = aln, width = grid::unit(1, "null"))
  }
  track_rows <- lapply(shown, function(track) {
    lapply(names(columns), function(name) {
      if (name == "aln") track_plot(track, layout, letter_size)
      else patchwork::plot_spacer()
    })
  })
  plots <- c(unlist(track_rows, recursive = FALSE),
             lapply(columns, function(col) col$plot))
  # the alignment's pixel height, capped at a screenful, so tracks over a
  # tall alignment keep a readable share of the figure
  heights <- c(vapply(shown, function(t) t$height, numeric(1)),
               min(max(nrow(layout$rows), 1) * layout$row_height, 800))
  patchwork::wrap_plots(
    plots, ncol = length(columns), byrow = TRUE, heights = heights,
    widths = do.call(grid::unit.c, lapply(columns, function(col) col$width)),
    guides = "collect"
  )
}

# the layout as the viewer draws it with its domains switched off: no spans,
# and the cells take the scheme's fill again unless the letters carry it
without_overlay <- function(layout) {
  layout$features <- NULL
  layout$domains <- NULL
  layout$letters_over_domains <- FALSE
  layout$fill_cells <- !layout$letter_mode
  layout$legends <- Filter(function(l) !identical(l$title, "Domains"),
                           layout$legends)
  layout
}

# The overlay and a features panel often color the same genes the same way,
# under the default title and the panel's field. Titled alike, patchwork
# collects the two keys as one.
merge_feature_keys <- function(layout) {
  overlay <- layout$features$shapes
  if (is.null(overlay)) return(layout)
  pairs <- function(shapes) unique(paste(shapes$key, shapes$fill))
  for (panel in layout$row_panels) {
    if (!is.null(panel$shapes) &&
        setequal(pairs(overlay), pairs(panel$shapes))) {
      layout$features$title <- panel$title
      break
    }
  }
  layout
}

pick_tracks <- function(all, tracks) {
  if (isTRUE(tracks)) {
    return(Filter(function(t) t$on, all))
  }
  if (isFALSE(tracks) || is.null(tracks)) {
    return(list())
  }
  ids <- vapply(all, function(t) t$id, character(1))
  names <- vapply(all, function(t) t$name, character(1))
  unknown <- setdiff(tracks, c(ids, names))
  if (length(unknown) > 0) {
    stop("no track named ", paste0("'", unknown, "'", collapse = ", "),
         "; the viewer has ", paste0("'", ids, "'", collapse = ", "),
         call. = FALSE)
  }
  Filter(function(t) t$id %in% tracks || t$name %in% tracks, all)
}

#' @rdname msaview_plot
#' @param object A viewer from [msaview()].
#' @param ... Passed to [msaview_plot()].
#' @exportS3Method ggplot2::autoplot
autoplot.msaview <- function(object, ...) {
  msaview_plot(object, ...)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# pretty file column numbers at the drawn columns that carry them, staying off
# the last column, whose label the panel edge cuts
column_ruler <- function(file_columns) {
  numbers <- pretty(range(file_columns))
  at <- match(numbers, file_columns)
  keep <- !is.na(at) & at <= length(file_columns) * 0.96
  list(breaks = at[keep], labels = numbers[keep])
}

#' @exportS3Method knitr::knit_print
knit_print.msaview <- function(x, ..., options = NULL) {
  if (knitr::is_html_output() ||
      !requireNamespace("QuickJSR", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    return(NextMethod())
  }
  knitr::knit_print(msaview_plot(x), ..., options = options)
}

row_scale <- function(layout, labels) {
  n <- nrow(layout$rows)
  ggplot2::scale_y_continuous(
    limits = c(-n - 0.5, -0.5), expand = c(0, 0),
    breaks = if (labels) -seq_len(n) else NULL,
    labels = if (labels) layout$rows$name else NULL
  )
}

bare_theme <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )
}

alignment_plot <- function(layout, letters, letter_size, strip = FALSE) {
  n <- nrow(layout$rows)
  m <- layout$num_columns
  aes <- ggplot2::aes
  letter_mat <- letter_matrix(layout)
  fills <- cell_colors(layout, letter_mat, "fill")
  matches <- reference_matches(layout, letter_mat)
  draw_letters <- letters

  ruler <- column_ruler(layout$file_columns %||% seq_len(m))
  p <- ggplot2::ggplot() + bare_theme() + row_scale(layout, FALSE) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(limits = c(0.5, m + 0.5), expand = c(0, 0),
                                breaks = ruler$breaks, labels = ruler$labels,
                                position = "top")
  if (layout$fill_cells && n > 0 && m > 0) {
    raster <- fills
    raster[matches] <- "#0000000a"
    raster[is.na(raster)] <- "transparent"
    p <- p + ggplot2::annotation_raster(
      grDevices::as.raster(raster), xmin = 0.5, xmax = m + 0.5,
      ymin = -n - 0.5, ymax = -0.5, interpolate = FALSE
    )
  }

  p <- p + tint_layer(layout)
  domains <- layout$domains
  p <- p + span_layers(layout$features, 3, layout,
                       x = function(x) x / layout$col_width + 0.5)

  if (draw_letters && n > 0 && m > 0) {
    p <- p + ggplot2::geom_text(
      data = letter_cells(layout, letter_mat, matches),
      aes(x = .data$x, y = .data$y, label = .data$letter, colour = .data$color),
      size = letter_size, family = "sans"
    )
  }

  p <- p + highlight_layers(layout, strip) + clade_layers(layout, m)
  p
}

# the rowTint channel's wash across a row, its alpha already applied
tint_layer <- function(layout) {
  rows <- which(!is.na(layout$rows$tint))
  if (length(rows) == 0) return(list())
  ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -rows - 0.5,
                    ymax = -rows + 0.5, fill = layout$rows$tint[rows])
}

# the tip labels, the row names or the rowLabel channel's field,
# right-aligned against the panels beside them, each in the color the tipLabel
# channel gives it, over its row's tint
labels_plot <- function(layout) {
  n <- nrow(layout$rows)
  d <- data.frame(
    y = -seq_len(n), name = layout$rows$label,
    color = ifelse(is.na(layout$rows$label_color), "grey20",
                   layout$rows$label_color)
  )
  ggplot2::ggplot() + bare_theme() + row_scale(layout, FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(axis.text = ggplot2::element_blank()) +
    tint_layer(layout) +
    ggplot2::geom_text(
      data = d, ggplot2::aes(x = 1, y = .data$y, label = .data$name),
      colour = d$color, hjust = 1, size = min(2.8, 120 / max(n, 1)),
      family = "sans"
    )
}

# each letter's color as renderMSABlock picks it: against the domain box under
# it, the reference wash, its own colored cell, or the plain background, and
# in letter-color mode the scheme's color darkened until it reads
letter_cells <- function(layout, letter_mat, matches) {
  primary <- "#000000de"
  color <- if (layout$fill_cells) {
    cell_colors(layout, letter_mat, "text")
  } else if (layout$letter_mode) {
    cell_colors(layout, letter_mat, "letterColor")
  } else {
    matrix(primary, nrow(letter_mat), ncol(letter_mat))
  }
  color[is.na(color) | matches] <- primary
  domains <- layout$domains
  if (layout$letters_over_domains && !is.null(domains)) {
    for (i in seq_len(nrow(domains))) {
      cols <- seq(domains$start[i], min(domains$end[i], ncol(color)))
      color[domains$row[i], cols] <- domains$text[i]
    }
  }
  cells <- data.frame(
    x = as.vector(col(letter_mat)), y = -as.vector(row(letter_mat)),
    letter = as.vector(letter_mat), color = as.vector(color),
    match = as.vector(matches), stringsAsFactors = FALSE
  )
  cells$letter[cells$match] <- "."
  cells[!is.na(cells$letter), ]
}

# the span mark as the bridge laid it out, x and y through the caller's
# mappings from pixels or rows, filled by the value its scale reads
span_layers <- function(spans, letter_size, layout, x,
                        y = function(y) -y - 0.5) {
  if (is.null(spans$shapes)) return(list())
  shapes <- spans$shapes
  shapes$x <- x(shapes$x)
  shapes$y <- y(shapes$y)
  fill <- fill_mapping(spans$title, shapes$key, shapes$fill, layout)
  shapes$key <- fill$key
  layers <- list(
    ggplot2::geom_polygon(
      data = shapes,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group,
                   fill = .data$key),
      colour = "grey25", linewidth = 0.25, show.legend = TRUE
    ),
    fill$scale
  )
  labels <- spans$labels
  if (!is.null(labels)) {
    labels$x <- x(labels$x)
    labels$y <- y(labels$y)
    size <- letter_size * labels$size * layout$row_height / layout$font_size
    layers <- c(layers, list(ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      colour = labels$color, size = size, family = "sans"
    )))
  }
  layers
}

# The fill scale a mark's keys map through, from the legend titled alike, and
# the keys as that scale reads them: numbers for a ramp, strings otherwise.
# `aesthetics = "colour"` gives a point, a text or an arc the same scale.
fill_mapping <- function(title, key, fill, layout, aesthetics = "fill") {
  legend <- Find(function(l) identical(l$title, title), layout$legends)
  if (!is.null(legend$gradient)) {
    list(key = suppressWarnings(as.numeric(key)),
         scale = gradient_scale(legend, aesthetics))
  } else {
    list(key = key, scale = fill_scale(title, key, fill, legend, aesthetics))
  }
}

# A linear or log scale as the viewer's ramp: its stops between its two ends,
# a value past an end in the end color, and the ends labeled as the viewer
# labels them, with a <= or >= where the field's values run past them
gradient_scale <- function(legend, aesthetics = "fill") {
  g <- legend$gradient
  ggplot2::scale_fill_gradientn(
    aesthetics = aesthetics,
    name = legend$title, colours = g$color, values = g$offset,
    limits = c(g$min, g$max), breaks = c(g$min, g$max),
    labels = c(g$min_label, g$max_label), oob = squish,
    transform = if (g$log) "log10" else "identity",
    na.value = missing_fill(legend)
  )
}

squish <- function(x, range, ...) pmin(pmax(x, range[1]), range[2])

# the grey a record with no value draws in, which a legend lists under the
# empty value when it names that row
missing_fill <- function(legend) {
  grey <- legend$color[legend$value == ""]
  if (length(grey) > 0) grey[[1]] else "grey70"
}

# a manual scale the viewer's colors seed, listed in the order of its legend
# and then by name, exon-2 before exon-10, so a scale_fill_*() added after it
# replaces it as it would any other. The breaks are the values the legend
# keys, labeled with its labels. Where keys share a color, as the alternating
# shades of a gene's exons do, a key would not tell them apart, so the scale
# draws none.
fill_scale <- function(title, key, fill, legend, aesthetics = "fill") {
  first <- !duplicated(key) & !is.na(key)
  values <- stats::setNames(fill[first], key[first])
  # every entry of the legend, so two panels sharing it build the same guide
  # and patchwork collects them as one
  if (!is.null(legend)) {
    values <- c(stats::setNames(legend$color, legend$value),
                values[setdiff(names(values), legend$value)])
  }
  # a legend lists what the viewer's key lists, so a panel whose scale maps a
  # few values keeps the rest out of the key
  breaks <- if (!is.null(legend)) legend$value else
    names(values)[natural_order(names(values))]
  # limits as well as breaks, since a manual scale drops the entries a panel's
  # own data lacks, and then two panels' guides differ and never merge; the
  # limits keep the data's own keys too, which color but stay out of the key
  ggplot2::scale_fill_manual(
    aesthetics = aesthetics,
    name = title, values = values, breaks = breaks,
    labels = if (!is.null(legend)) legend$label else ggplot2::waiver(),
    limits = union(breaks, names(values)),
    na.value = if (!is.null(legend)) missing_fill(legend) else "grey70",
    guide = if (anyDuplicated(values[breaks])) "none" else "legend"
  )
}

natural_order <- function(x) {
  number <- suppressWarnings(as.numeric(sub("^\\D*(\\d+).*$", "\\1", x)))
  order(sub("\\d.*$", "", x), number, x)
}

# A panel beside the tree, x in its own pixels. A panel with an axis, a bar,
# a point or a features panel over positions, labels it above the panel with
# the header across over it, as the viewer's header band draws it.
row_panel_plot <- function(panel, layout, letter_size) {
  axis <- panel$axis
  p <- ggplot2::ggplot() + bare_theme() + row_scale(layout, FALSE) +
    ggplot2::scale_x_continuous(
      limits = c(0, panel$width), expand = c(0, 0), position = "top",
      breaks = axis$at,
      labels = if (!is.null(axis)) axis$label else ggplot2::waiver(),
      name = panel$header
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.text.x.top = if (!is.null(axis)) {
        ggplot2::element_text(size = 6.5)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x.top = if (!is.null(axis)) {
        ggplot2::element_line(colour = "grey30", linewidth = 0.3)
      } else {
        ggplot2::element_blank()
      },
      axis.line.x.top = if (!is.null(axis)) {
        ggplot2::element_line(colour = "grey30", linewidth = 0.3)
      } else {
        ggplot2::element_blank()
      },
      # a header wider than its strip turns on its side, as the viewer's does
      axis.title.x.top = if (is.null(axis) &&
                             nchar(panel$header) * 6 > panel$width) {
        ggplot2::element_text(size = 9, angle = 90, hjust = 0, vjust = 0.5)
      } else {
        ggplot2::element_text(size = 9)
      }
    )
  cells <- panel$cells
  if (!is.null(cells)) {
    fill <- fill_mapping(panel$title, cells$key, cells$color, layout)
    cells$key <- fill$key
    p <- p +
      ggplot2::geom_tile(
        data = cells,
        ggplot2::aes(x = panel$width / 2, y = -.data$row, fill = .data$key),
        width = panel$width, height = 1, show.legend = TRUE
      ) +
      fill$scale
  }
  p + value_marks(panel, layout, letter_size) +
    span_layers(panel, letter_size, layout, x = identity)
}

# A bar or point panel's marks, ggtree's geom_facet with geom_col or
# geom_point: a bar from the panel's base to its value, 0.8 of the row high
# as the viewer leaves a gap between rows, and a dot at its value, `size`
# pixels across. A text panel prints each row's text from its left edge.
# Marks a color field reads map their key through its scale, and the rest
# draw grey, or a text in the theme's text color.
value_marks <- function(panel, layout, letter_size) {
  marks <- panel$marks %||% panel$texts
  if (is.null(marks)) return(list())
  aesthetic <- if (panel$kind == "bar") "fill" else "colour"
  m <- keyed_marks(marks, panel$title, layout, aesthetic)
  own <- m$own
  own$color[is.na(own$color)] <- "#000000de"
  geom <- switch(panel$kind,
    bar = list(
      fn = ggplot2::geom_rect,
      aes = ggplot2::aes(xmin = pmin(.data$base, .data$x),
                         xmax = pmax(.data$base, .data$x),
                         ymin = -.data$row - 0.4, ymax = -.data$row + 0.4),
      args = list()
    ),
    point = list(
      fn = ggplot2::geom_point,
      aes = ggplot2::aes(x = .data$x, y = -.data$row),
      args = list(size = panel$size * letter_size / layout$font_size)
    ),
    text = list(
      fn = ggplot2::geom_text,
      aes = ggplot2::aes(x = panel$padding, y = -.data$row,
                         label = .data$text),
      args = list(hjust = 0, family = "sans",
                  size = min(2.8, 120 / max(nrow(layout$rows), 1)))
    )
  )
  keyed_aes <- geom$aes
  keyed_aes[[aesthetic]] <- quote(.data$key)
  list(
    if (nrow(own) > 0) {
      do.call(geom$fn, c(list(data = own, mapping = geom$aes),
                         stats::setNames(list(own$color), aesthetic),
                         geom$args))
    },
    if (nrow(m$mapped) > 0) {
      do.call(geom$fn, c(list(data = m$mapped, mapping = keyed_aes,
                              show.legend = TRUE),
                         geom$args))
    },
    m$scale
  )
}

# Marks split by whether a scale colored them: `mapped` carries each key as
# the scale the legend titled `title` reads it, through `scale`, and `own`
# draws in its own color
keyed_marks <- function(marks, title, layout, aesthetic = "fill") {
  keyed <- !is.na(marks$key)
  mapped <- marks[keyed, ]
  scale <- if (nrow(mapped) > 0) {
    fill_mapping(title, mapped$key, mapped$color, layout, aesthetic)
  }
  if (!is.null(scale)) mapped$key <- scale$key
  list(own = marks[!keyed, ], mapped = mapped, scale = scale$scale)
}

# the clade bracket mark: a bar beside each clade's rows and its label, read
# across the rows where they are taller than the font and up the bar where not
gutter_plot <- function(gutter, layout, letter_size) {
  p <- ggplot2::ggplot() + bare_theme() + row_scale(layout, FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, gutter$width), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(axis.text = ggplot2::element_blank())
  bars <- gutter$bars
  if (!is.null(bars)) {
    bars$color[is.na(bars$color)] <- "#000000de"
    p <- p + ggplot2::geom_segment(
      data = bars,
      ggplot2::aes(x = gutter$bar_x, xend = gutter$bar_x,
                   y = -.data$first + 0.45, yend = -.data$last - 0.45),
      colour = bars$color, linewidth = 1, lineend = "butt"
    )
  }
  labels <- gutter$labels
  if (!is.null(labels)) {
    labels$color[is.na(labels$color)] <- "#000000de"
    across <- labels[labels$horizontal, ]
    up <- labels[!labels$horizontal, ]
    if (nrow(across) > 0) {
      p <- p + ggplot2::geom_text(
        data = across,
        ggplot2::aes(x = .data$left, y = -0.5 - .data$top - .data$height / 2,
                     label = .data$label),
        colour = across$color, hjust = 0, size = letter_size, family = "sans"
      )
    }
    if (nrow(up) > 0) {
      p <- p + ggplot2::geom_text(
        data = up,
        ggplot2::aes(x = .data$left, y = -0.5 - .data$top - .data$height,
                     label = .data$label),
        colour = up$color, hjust = 0, vjust = 1, angle = 90, size = letter_size,
        family = "sans"
      )
    }
  }
  p
}

# With the highlight strip shown, a bracket-mark span is an outline in its own
# color, opaque, and its label is in the strip; every other highlight is a band
highlight_layers <- function(layout, strip = FALSE) {
  h <- layout$highlights
  if (is.null(h)) return(list())
  bracketed <- strip & !is.na(h$start) & h$mark %in% "bracket"
  outlined <- h[bracketed, ]
  h <- h[!bracketed, ]
  outlines <- if (nrow(outlined) > 0) highlight_outlines(outlined, layout)
  if (nrow(h) == 0) return(list(outlines))
  n <- nrow(layout$rows)
  m <- layout$num_columns
  h$xmin <- ifelse(is.na(h$start), 0.5, h$start - 0.5)
  h$xmax <- ifelse(is.na(h$end), m + 0.5, h$end + 0.5)
  h$ymin <- ifelse(is.na(h$row), -n - 0.5, -h$row - 0.5)
  h$ymax <- ifelse(is.na(h$row), -0.5, -h$row + 0.5)
  # a color carrying its own alpha draws at it, one without at 0.28
  own <- ifelse(is.na(h$color), "#ff8c00", h$color)
  fill <- ifelse(grDevices::col2rgb(own, alpha = TRUE)["alpha", ] < 255, own,
                 grDevices::adjustcolor(own, alpha.f = 0.28))
  stroke <- ifelse(is.na(h$color), "#d25a00", h$color)
  labelled <- h[!is.na(h$label), ]
  list(
    ggplot2::geom_rect(
      data = h,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = fill, colour = stroke, linewidth = 0.6
    ),
    if (nrow(labelled) > 0) highlight_labels(labelled),
    outlines
  )
}

highlight_outlines <- function(h, layout) {
  n <- nrow(layout$rows)
  h$xmin <- h$start - 0.5
  h$xmax <- h$end + 0.5
  stroke <- ifelse(is.na(h$color), "grey45", opaque(h$color))
  ggplot2::geom_rect(
    data = h,
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
    ymin = -n - 0.5, ymax = -0.5, fill = NA, colour = stroke, linewidth = 0.3
  )
}

opaque <- function(color) {
  grDevices::rgb(t(grDevices::col2rgb(color)), maxColorValue = 255)
}

# a label at the top right of its band; with ggrepel, labels of neighbouring
# bands push each other apart instead of printing over one another
highlight_labels <- function(labelled) {
  aes <- ggplot2::aes(x = .data$xmax, y = .data$ymax, label = .data$label)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_label_repel(
      data = labelled, mapping = aes, size = 2.8, linewidth = 0,
      fill = grDevices::adjustcolor("white", 0.85), direction = "y",
      hjust = 0, min.segment.length = 0, max.overlaps = Inf, seed = 1
    )
  } else {
    ggplot2::geom_label(
      data = labelled, mapping = aes, hjust = 0, vjust = 1, size = 2.8,
      label.size = 0, fill = "white", alpha = 0.85
    )
  }
}

clade_layers <- function(layout, m) {
  cl <- layout$clades
  if (is.null(cl)) return(list())
  cl <- cl[cl$mark == "highlight", ]
  if (nrow(cl) == 0) return(list())
  list(ggplot2::geom_rect(
    data = cl,
    ggplot2::aes(xmin = 0.5, xmax = m + 0.5,
                 ymin = -.data$last - 0.5, ymax = -.data$first + 0.5),
    fill = cl$color, inherit.aes = FALSE
  ))
}

tree_plot <- function(layout) {
  t <- layout$tree
  n <- nrow(layout$rows)
  # the rows' limits go on the coord, which clips nothing, so the scale bar can
  # sit above the first row
  p <- ggplot2::ggplot() + bare_theme() +
    ggplot2::scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(-n - 0.5, -0.5), clip = "off") +
    ggplot2::theme(axis.text = ggplot2::element_blank())
  cl <- layout$clades
  if (!is.null(cl)) {
    marked <- cl[cl$mark == "highlight", ]
    if (nrow(marked) > 0) {
      p <- p + ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf,
        ymin = -marked$last - 0.5, ymax = -marked$first + 0.5,
        fill = marked$color
      )
    }
  }
  p <- p + tint_layer(layout)
  if (!is.null(t$segments)) {
    # each edge in the branch channel's color, black where it gives none
    color <- t$segments$color %||% rep(NA_character_, nrow(t$segments))
    p <- p + ggplot2::geom_segment(
      data = t$segments,
      ggplot2::aes(x = .data$x, y = -.data$y, xend = .data$xend,
                   yend = -.data$yend),
      colour = ifelse(is.na(color), "black", color),
      linewidth = 0.4, lineend = "square"
    )
  }
  if (!is.null(t$collapsed)) {
    p <- p + ggplot2::geom_polygon(
      data = t$collapsed,
      ggplot2::aes(x = .data$x, y = -.data$y, group = .data$group),
      fill = "grey70", colour = "grey30", linewidth = 0.3
    )
  }
  p <- p + node_points(t, layout)
  bar <- t$scale_bar
  if (!is.null(bar)) {
    p <- p +
      ggplot2::annotate("segment", x = 0, xend = bar$step, y = 0, yend = 0,
                        linewidth = 0.4) +
      ggplot2::annotate("text", x = bar$step / 2, y = 0.15, label = bar$label,
                        vjust = 0, size = 2.5, family = "sans")
  }
  p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02)))
}

# the bubble on each internal node the nodePoint channel colors, filled by
# its support through the channel's scale, ggtree's geom_nodepoint
node_points <- function(tree, layout) {
  nodes <- tree$nodes
  if (is.null(nodes)) return(list())
  fill <- fill_mapping(tree$node_title, nodes$key, nodes$color, layout)
  nodes$key <- fill$key
  list(
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .data$x, y = -.data$y, fill = .data$key),
      shape = 21, size = 2, colour = "grey20", stroke = 0.3,
      show.legend = TRUE
    ),
    fill$scale
  )
}

track_plot <- function(track, layout, letter_size = 2.5) {
  m <- layout$num_columns
  p <- switch(track$kind,
    bar = bar_track(track, layout),
    strip = strip_track(track, layout),
    logo = logo_track(track),
    arc = arc_track(track, layout),
    text = text_track(track, layout, letter_size),
    highlights = highlight_track(track),
    features = feature_track(track, layout),
    ggplot2::ggplot()
  )
  p +
    ggplot2::scale_x_continuous(limits = c(0.5, m + 0.5), expand = c(0, 0)) +
    bare_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 8, hjust = 0,
                                         margin = ggplot2::margin(b = 2))
    ) +
    ggplot2::labs(title = track$name)
}

# each bar's height a fraction of the track's, as the viewer's value scale
# maps it, in the color a scale over the track's values gives it
bar_track <- function(track, layout) {
  p <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  if (is.null(track$bars)) return(p)
  m <- keyed_marks(track$bars, track$title, layout)
  p + list(
    if (nrow(m$own) > 0) {
      ggplot2::geom_col(data = m$own, ggplot2::aes(x = .data$col, y = .data$value),
                        fill = m$own$color, width = 1)
    },
    if (nrow(m$mapped) > 0) {
      ggplot2::geom_col(
        data = m$mapped,
        ggplot2::aes(x = .data$col, y = .data$value, fill = .data$key),
        width = 1, show.legend = TRUE, key_glyph = "rect"
      )
    },
    m$scale
  )
}

# a cell per column the full track height, a heatmap row of the track's values
strip_track <- function(track, layout) {
  p <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  if (is.null(track$cells)) return(p)
  m <- keyed_marks(track$cells, track$title, layout)
  p + list(
    if (nrow(m$own) > 0) {
      ggplot2::geom_tile(data = m$own, ggplot2::aes(x = .data$col, y = 0.5),
                         width = 1, height = 1, fill = m$own$color)
    },
    if (nrow(m$mapped) > 0) {
      ggplot2::geom_tile(
        data = m$mapped, ggplot2::aes(x = .data$col, y = 0.5, fill = .data$key),
        width = 1, height = 1, show.legend = TRUE, key_glyph = "rect"
      )
    },
    m$scale
  )
}

# each bracket-mark highlight as the strip draws it: a bar with end ticks over
# its columns, or a triangle over a single residue, one lane per nesting
# level, the label above; a label that would overprint another is dropped
highlight_track <- function(track) {
  b <- track$bands
  if (is.null(b)) return(ggplot2::ggplot())
  lanes <- max(b$lane) + 1
  b$color <- ifelse(is.na(b$color), "grey45", b$color)
  b$y <- lanes - b$lane - 0.85
  b$x0 <- b$start - 0.5
  b$x1 <- b$end + 0.5
  spans <- b[b$start < b$end, ]
  singles <- b[b$start == b$end, ]
  ticks <- rbind(
    data.frame(x = spans$x0, y = spans$y, color = spans$color),
    data.frame(x = spans$x1, y = spans$y, color = spans$color)
  )
  labelled <- b[!is.na(b$label), ]
  labelled$x <- ifelse(labelled$start == labelled$end, labelled$start + 0.8,
                       labelled$x0 + 0.3)
  ggplot2::ggplot() +
    ggplot2::scale_y_continuous(limits = c(0, lanes), expand = c(0, 0)) +
    (if (nrow(spans) > 0) list(
      ggplot2::geom_segment(
        data = spans,
        ggplot2::aes(x = .data$x0, xend = .data$x1, y = .data$y,
                     yend = .data$y),
        colour = spans$color, linewidth = 0.4
      ),
      ggplot2::geom_segment(
        data = ticks,
        ggplot2::aes(x = .data$x, xend = .data$x, y = .data$y,
                     yend = .data$y + 0.3),
        colour = ticks$color, linewidth = 0.4
      )
    )) +
    (if (nrow(singles) > 0) ggplot2::geom_point(
      data = singles, ggplot2::aes(x = .data$start, y = .data$y + 0.15),
      shape = 25, size = 1.6, colour = singles$color, fill = singles$color
    )) +
    (if (nrow(labelled) > 0) ggplot2::geom_text(
      data = labelled, ggplot2::aes(x = .data$x, y = .data$y + 0.34,
                                    label = .data$label),
      colour = labelled$color, hjust = 0, vjust = 0, size = 2.4,
      family = "sans", check_overlap = TRUE
    ))
}

# a features track's spans as the span mark draws them, the lanes counted down
# from the track's top
feature_track <- function(track, layout) {
  ggplot2::ggplot() +
    ggplot2::scale_y_continuous(limits = c(0, track$height), expand = c(0, 0)) +
    span_layers(track, 3, layout,
                x = function(x) x / layout$col_width + 0.5,
                y = function(y) track$height - y)
}

# an arc's apex grows with the square root of its span, as the viewer's does;
# an arc a scale over the values colors maps its key through that scale
arc_track <- function(track, layout) {
  a <- track$arcs
  p <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  if (is.null(a)) return(p)
  a$color[is.na(a$color)] <- "#000000de"
  span <- max(a$end - a$start, 1)
  t <- seq(0, 1, length.out = 33)
  curves <- do.call(rbind, lapply(seq_len(nrow(a)), function(i) {
    apex <- 0.95 * sqrt((a$end[i] - a$start[i]) / span)
    data.frame(
      x = a$start[i] + (a$end[i] - a$start[i]) * t,
      y = 4 * apex * t * (1 - t),
      group = i, key = a$key[i], color = a$color[i]
    )
  }))
  m <- keyed_marks(curves, track$title, layout, "colour")
  p + list(
    if (nrow(m$own) > 0) {
      ggplot2::geom_path(
        data = m$own,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
        colour = m$own$color, linewidth = 0.3
      )
    },
    if (nrow(m$mapped) > 0) {
      ggplot2::geom_path(
        data = m$mapped,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$group,
                     colour = .data$key),
        linewidth = 0.3, show.legend = TRUE
      )
    },
    m$scale
  )
}

# each character over its cell, the cell filled where the track's scale or the
# residue scheme colors it
text_track <- function(track, layout, size) {
  cells <- track$cells
  p <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  if (is.null(cells)) return(p)
  filled <- cells[!is.na(cells$fill), ]
  filled$color <- filled$fill
  m <- keyed_marks(filled, track$title, layout)
  rect <- ggplot2::aes(xmin = .data$col - 0.5, xmax = .data$col + 0.5,
                       ymin = 0, ymax = 1)
  keyed_rect <- rect
  keyed_rect$fill <- quote(.data$key)
  p + list(
    if (nrow(m$own) > 0) {
      ggplot2::geom_rect(data = m$own, rect, fill = m$own$fill)
    },
    if (nrow(m$mapped) > 0) {
      ggplot2::geom_rect(data = m$mapped, keyed_rect, show.legend = TRUE)
    },
    m$scale,
    ggplot2::geom_text(
      data = cells, ggplot2::aes(x = .data$col, y = 0.5, label = .data$letter),
      colour = cells$text, size = size, family = "sans"
    )
  )
}

logo_track <- function(track) {
  l <- track$letters
  p <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(limits = c(0, track$max_bits), expand = c(0, 0))
  if (is.null(l)) return(p)
  l <- l[l$ymax - l$ymin > track$max_bits / 200, ]
  shapes <- logo_polygons(l)
  if (is.null(shapes)) {
    return(p + ggplot2::geom_rect(
      data = l,
      ggplot2::aes(xmin = .data$col - 0.45, xmax = .data$col + 0.45,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = ifelse(is.na(l$color), "grey40", l$color)
    ))
  }
  p + ggplot2::geom_polygon(
    data = shapes,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$group,
                 subgroup = .data$contour),
    fill = shapes$fill
  )
}

# each letter's outline in the system sans font, scaled to its slice of the
# column the way the viewer stretches a glyph; NULL without systemfonts
logo_polygons <- function(l) {
  outlines <- lapply(unique(l$letter), glyph_unit)
  names(outlines) <- unique(l$letter)
  if (nrow(l) == 0 || any(vapply(outlines, is.null, logical(1)))) return(NULL)
  pieces <- lapply(seq_len(nrow(l)), function(i) {
    o <- outlines[[l$letter[i]]]
    data.frame(
      x = l$col[i] - 0.45 + 0.9 * o$x,
      y = l$ymin[i] + (l$ymax[i] - l$ymin[i]) * o$y,
      contour = o$contour, group = i,
      fill = if (is.na(l$color[i])) "grey40" else l$color[i]
    )
  })
  do.call(rbind, pieces)
}

glyph_unit <- function(letter) {
  if (!is.null(the$glyphs[[letter]])) return(the$glyphs[[letter]])
  if (!requireNamespace("systemfonts", quietly = TRUE) ||
      !exists("glyph_outline", asNamespace("systemfonts"))) {
    return(NULL)
  }
  font <- systemfonts::match_fonts("sans", weight = "bold")
  info <- systemfonts::glyph_info(letter, path = font$path, index = font$index)
  o <- tryCatch(
    systemfonts::glyph_outline(info$index[1], font$path, font$index, size = 100),
    error = function(e) NULL
  )
  if (is.null(o) || nrow(o) == 0) return(NULL)
  unit <- data.frame(
    x = (o$x - min(o$x)) / diff(range(o$x)),
    y = (o$y - min(o$y)) / diff(range(o$y)),
    contour = o$contour
  )
  if (is.null(the$glyphs)) the$glyphs <- list()
  the$glyphs[[letter]] <- unit
  unit
}
