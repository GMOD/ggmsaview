#' Add a layer to an MSA viewer
#'
#' Composes a viewer with a layer, in the way ggplot2 and ggtree compose a plot
#' with a geom. \code{msaview()} takes every layer as an argument too, and these
#' build the same viewer one piece at a time, each call reading a data frame.
#'
#' @param e1 An \code{msaview} htmlwidget.
#' @param e2 A layer from \code{\link{geom_msa_highlight}},
#'   \code{\link{geom_msa_clade}}, \code{\link{geom_msa_track}},
#'   \code{\link{geom_msa_domains}},
#'   \code{\link{geom_msa_rowdata}}, \code{\link{geom_msa_strip}},
#'   \code{\link{geom_msa_features}}, \code{\link{geom_msa_bar}},
#'   \code{\link{geom_msa_point}}, \code{\link{geom_msa_text}},
#'   \code{\link{geom_msa_nodepoint}}, \code{\link{geom_msa_tiplab}},
#'   \code{\link{scale_residue_color}},
#'   \code{\link{scale_row_color}}, \code{\link{theme_msa}},
#'   \code{\link{coord_msa}}, \code{\link{geom_msa_selection}},
#'   \code{\link{filter_msa_rows}}, \code{\link{arrange_msa_rows}} or
#'   \code{\link{facet_msa_rows}}.
#' @return The viewer, carrying the layer.
#'
#' @examples
#' \dontrun{
#' library(ggmsaview)
#'
#' variants <- data.frame(
#'   row = c("Human", "Human"),
#'   start = c(248, 273),
#'   end = c(248, 273),
#'   label = c("R248Q", "R273H")
#' )
#'
#' msaview(msa = "p53.aln", tree = "p53.nh") +
#'   geom_msa_domains("p53-domains.gff") +
#'   geom_msa_highlight(variants) +
#'   scale_residue_color("clustal", channel = "color") +
#'   theme_msa("dark", col_width = 14)
#' }
#' @export
"+.msaview" <- function(e1, e2) {
  if (!inherits(e2, "msaview_layer")) {
    stop("the right side of + must be an ggmsaview layer, such as ",
         "geom_msa_track() or theme_msa(), got a '", class(e2)[1], "'")
  }
  for (prop in names(e2$set)) {
    e1$x$props[[prop]] <- e2$set[[prop]]
  }
  # an unnamed list of one serializes as a JSON array, so appending needs no I()
  for (prop in names(e2$append)) {
    e1$x$props[[prop]] <- c(e1$x$props[[prop]], e2$append[[prop]])
  }
  for (prop in names(e2$into)) {
    record <- e1$x$props[[prop]] %||% list()
    for (member in names(e2$into[[prop]]$set)) {
      record[[member]] <- e2$into[[prop]]$set[[member]]
    }
    for (member in names(e2$into[[prop]]$append)) {
      record[[member]] <- c(record[[member]], e2$into[[prop]]$append[[member]])
    }
    e1$x$props[[prop]] <- record
  }
  e1
}

# `set` replaces a prop, `append` adds to a list-valued one, and `into` does
# either to the members of an object-valued prop, as
# `list(rowTransform = list(set = list(facet = ...)))`. A NULL in `set` would
# drop the prop, so the constructors leave it out instead.
msa_layer <- function(set = list(), append = list(), into = list()) {
  structure(list(set = drop_null(set), append = append, into = into),
            class = "msaview_layer")
}

drop_null <- function(x) x[!vapply(x, is.null, logical(1))]

#' @export
print.msaview_layer <- function(x, ...) {
  props <- c(names(x$set), names(x$append), names(x$into))
  cat("<ggmsaview layer>", paste(props, collapse = ", "), "\n")
  invisible(x)
}

# Reads one field off a data frame row, or NULL when the column is absent or
# the cell is NA, so a frame carrying only start and end produces highlights
# with only those fields, and an NA color takes the default.
field_of <- function(df, i, name) {
  value <- if (name %in% names(df)) df[[name]][[i]]
  if (length(value) == 1 && is.na(value)) NULL else value
}

#' Highlight columns, residues or rows
#'
#' Marks a span with a bordered band and an optional label. A data frame gives
#' one highlight per row, read from the columns \code{start}, \code{end},
#' \code{row}, \code{label}, \code{color} and \code{mark}, whichever it
#' carries. The arguments give one highlight on their own.
#'
#' Coordinates are 1-based and inclusive. \code{start} and \code{end} alone
#' count alignment columns; with \code{row} they count that row's own residues,
#' and the viewer places the band on the columns those residues land in.
#'
#' @param data A data frame of highlights, or NULL to use the arguments.
#' @param columns Alignment columns (1-based) to highlight, as a numeric
#'   vector. The viewer draws each run of adjacent columns as one band.
#' @param start,end The span, 1-based and inclusive.
#' @param row A row name, making the span residues of that row.
#' @param rows Row names to tint whole rows.
#' @param label Text drawn beside the band.
#' @param color A CSS color for the band.
#' @param mark \code{"band"} (the default) fills the span; \code{"bracket"}
#'   outlines it and draws a labeled bracket in a strip above the alignment,
#'   with nested spans stacked into lanes.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' variants <- data.frame(row = "Human", start = 248, end = 248, label = "R248Q")
#' msaview(msa = "p53.aln") + geom_msa_highlight(variants)
#'
#' msaview(msa = "p53.aln") + geom_msa_highlight(start = 100, end = 288)
#' }
#' @export
geom_msa_highlight <- function(data = NULL, start = NULL, end = NULL,
                               row = NULL, rows = NULL, label = NULL,
                               color = NULL, columns = NULL, mark = NULL) {
  if (!is.null(columns) && is.null(data) && is.null(start) && is.null(rows)) {
    return(msa_layer(
      set = list(highlightColumns = convert_highlight_columns(columns))
    ))
  }
  if (is.data.frame(data)) {
    if (!all(c("start", "end") %in% names(data))) {
      stop("a highlight data frame must have 'start' and 'end' columns, got ",
           paste(names(data), collapse = ", "))
    }
    highlights <- lapply(seq_len(nrow(data)), function(i) {
      drop_null(list(
        start = data$start[[i]],
        end = data$end[[i]],
        row = field_of(data, i, "row"),
        label = field_of(data, i, "label"),
        color = field_of(data, i, "color"),
        mark = field_of(data, i, "mark")
      ))
    })
  } else if (!is.null(data)) {
    stop("geom_msa_highlight takes a data frame, got a '", class(data)[1], "'")
  } else if (is.null(rows) && (is.null(start) || is.null(end))) {
    stop("geom_msa_highlight needs a data frame, a start and end, or rows")
  } else {
    highlights <- list(drop_null(list(
      start = start, end = end, row = row, rows = rows,
      label = label, color = color, mark = mark
    )))
  }
  msa_layer(
    append = list(highlights = convert_highlights(highlights)),
    set = list(highlightColumns = convert_highlight_columns(columns))
  )
}

#' Mark a clade of the tree
#'
#' \code{mark = "highlight"} fills the rows of a clade with a translucent
#' rectangle, from the clade's common ancestor across the tree and the
#' alignment, which is ggtree's \code{geom_hilight}.
#' \code{mark = "bracket"} draws a bar beside the rows carrying
#' \code{label}, which is \code{geom_cladelab}, or \code{geom_strip} over a
#' \code{range}. \code{mark = "collapse"} and \code{mark = "focus"} open the
#' viewer with the clade collapsed, or with the rest of the tree hidden.
#' \code{mark = "rotate"} swaps the order of the clade's children, which is
#' ggtree's \code{rotate}. All three need \code{mrca}.
#'
#' \code{mrca} names tips whose most recent common ancestor is the clade, and
#' \code{range} names the first and last tip of a run in display order, which
#' need not be monophyletic. \code{tips} is the number of tips the clade
#' covers: a clade that resolves to a different number is dropped, so a
#' re-rooted tree loses the mark rather than putting it on another clade.
#'
#' @param mrca Tip names whose common ancestor is the clade.
#' @param tips The number of tips the clade covers.
#' @param range The first and last tip of a run, in display order.
#' @param color A CSS color for the rectangle, or for a bracket's bar and
#'   label.
#' @param label Text naming the clade, drawn by the bracket mark and by a
#'   highlight carrying one.
#' @param mark What to draw over the clade: \code{"highlight"},
#'   \code{"bracket"}, \code{"collapse"}, \code{"focus"} or
#'   \code{"rotate"}.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_clade(c("Gs/TW/TNC1/2015", "Ck/TW/a174/2015"), tips = 47,
#'                  color = "#fff3c4", label = "2.3.4.4 H5Nx") +
#'   geom_msa_clade(c("Gs/TW/TNC1/2015", "Ck/TW/a174/2015"), tips = 47,
#'                  mark = "bracket", label = "2.3.4.4 H5Nx")
#' }
#' @export
geom_msa_clade <- function(mrca = NULL, tips = NULL, range = NULL,
                           color = NULL, label = NULL,
                           mark = "highlight") {
  clade <- drop_null(list(
    mrca = mrca, range = range, tips = tips, mark = mark,
    color = color, label = label
  ))
  msa_layer(append = list(clades = convert_clades(list(clade))))
}

#' A track of your own numbers, text, arcs or spans
#'
#' Draws a track above the alignment, beside the conservation tracks the viewer
#' computes. A data frame supplies the values from the column \code{value}
#' names, and a numeric vector supplies them directly. \code{kind = "bar"}
#' draws each value as a bar, and \code{kind = "strip"} as a cell the full
#' track height, a heatmap row of a per-residue score such as pLDDT.
#'
#' A bar, strip, text or arc track colors its own values through a scale: a
#' bar's or a strip's numbers, a text track's characters, or the
#' \code{value} column of an arc data frame. The scale arguments
#' \code{\link{scale_row_color}} takes describe it, and the track's
#' \code{name} titles its legend unless \code{title} does. Without any of
#' them a bar, text or arc track draws in \code{color}, and a strip colors
#' each distinct value from the default palette. An arc's own \code{color}
#' column wins over the scale.
#'
#' A bar's height reads its value through a linear scale from
#' \code{domain_min} (default 0) to \code{domain_max} (default 1), or with
#' \code{log = TRUE} a log scale whose unset ends follow the track's positive
#' values, clipped at \code{domain_quantile}. On a bar track these three
#' arguments set the height scale, and a color ramp's ends follow the values;
#' on every other kind they set the color ramp's ends.
#'
#' \code{kind = "features"} draws spans with the mark
#' \code{\link{geom_msa_features}} draws in each row: an arrow where a span
#' has a \code{strand}, a box otherwise, labeled where the text fits. Its data
#' frame has one span per row, with \code{start} and \code{end}, 1-based and
#' inclusive, and optionally \code{strand} (\code{"+"} or \code{"-"}),
#' \code{color}, which wins over any scale, and other columns, which are
#' fields \code{fill} and \code{label} read. Spans that overlap pack into
#' lanes.
#'
#' \code{row} indexes the values or spans by the residues of that row. Without
#' it they index alignment columns.
#'
#' @param data A data frame holding the values or the spans, a numeric vector
#'   of values, or NULL when \code{text} or \code{arcs} carries the track.
#' @param name The track's label, shown beside it.
#' @param id A stable id for the track. Taken from \code{name} when absent.
#' @param kind \code{"bar"}, \code{"strip"}, \code{"text"}, \code{"arc"} or
#'   \code{"features"}.
#' @param value The column of \code{data} holding the values.
#' @param text One character per column, for a text track.
#' @param range The colors or palette name of the track's color scale, as
#'   \code{\link{scale_row_color}} takes them. A named vector, such as
#'   \code{c("1" = "#dddddd")} for a text track, colors the values it names,
#'   its names read as \code{domain}. For a features track, the scale of
#'   \code{fill}.
#' @param arcs Pairs joined by a curve, as a data frame with \code{start} and
#'   \code{end} columns and an optional \code{color} and \code{value}, the
#'   number the color scale reads.
#' @param row A row name, indexing the values by that row's residues.
#' @param domain_max For a bar track, the value drawn at full height, default
#'   1. For any other kind, the high end of a color ramp.
#' @param domain_min For a bar track, the value drawn at no height, default 0.
#'   For any other kind, the low end of a color ramp. Unset, the lowest value.
#' @param domain_quantile For a bar track, the quantile an unset end of a log
#'   height scale follows. For any other kind, the quantile an unset end of a
#'   color ramp follows.
#' @param log \code{TRUE} draws a bar track's heights on a log scale.
#' @param color The one color a bar, strip, text or arc track draws in when no
#'   scale argument is given, or the fill of every span the \code{fill} scale
#'   leaves uncolored.
#' @param height The track's pixel height.
#' @param fill The span field to fill a features track's spans by.
#' @param label The span field drawn inside each span.
#' @param position How a features track's overlapping spans share its height:
#'   \code{"identity"} (the default) packs them into lanes, and
#'   \code{"strandpile"} puts forward spans above reverse ones.
#' @inheritParams scale_row_color
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' dnds <- data.frame(value = runif(200, 0, 2))
#' msaview(msa = "cds.aln") +
#'   geom_msa_track(dnds, name = "dN/dS", domain_max = 2, color = "#c0392b")
#'
#' msaview(msa = "cds.aln") +
#'   geom_msa_track(dnds, name = "dN/dS", domain_max = 2, type = "threshold",
#'                  domain = 1, range = c("#2c7fb8", "#d7301f"))
#'
#' plddt <- data.frame(value = c(96.2, 94.8, 71.3, 48.0, 33.5))
#' msaview(msa = "insulin.aln") +
#'   geom_msa_track(plddt, name = "pLDDT", kind = "strip", row = "Human",
#'                  type = "threshold", domain = c(50, 70, 90),
#'                  range = c("#ff7d45", "#ffdb13", "#65cbf3", "#0053d6"),
#'                  labels = c("very low", "low", "confident", "very high"))
#'
#' contacts <- data.frame(start = c(10, 40), end = c(120, 160),
#'                        value = c(0.9, 0.4))
#' msaview(msa = "kinase.aln") +
#'   geom_msa_track(name = "contacts", kind = "arc", arcs = contacts,
#'                  type = "linear", range = "viridis")
#'
#' ss <- data.frame(start = c(33, 48, 91), end = c(43, 50, 97),
#'                  type = c("helix", "strand", "helix"))
#' msaview(msa = "insulin.aln") +
#'   geom_msa_track(ss, name = "Secondary structure", kind = "features",
#'                  row = "Human", fill = "type", label = "type",
#'                  range = c(helix = "#e15759", strand = "#edc948"))
#' }
#' @export
geom_msa_track <- function(data = NULL, name = NULL, id = NULL, kind = "bar",
                           value = "value", text = NULL, range = NULL,
                           arcs = NULL, row = NULL, domain_max = NULL,
                           color = NULL, height = NULL, fill = NULL,
                           label = NULL, position = NULL, domain = NULL,
                           title = NULL, type = NULL, domain_min = NULL,
                           domain_mid = NULL, domain_quantile = NULL,
                           reverse = NULL, labels = NULL,
                           missing_label = NULL, log = NULL) {
  if (is.null(name)) {
    stop("geom_msa_track needs a name")
  }
  if (!is.null(log) && kind != "bar") {
    stop("log sets a bar track's height scale; a color ramp takes ",
         "type = \"log\"")
  }
  values <- NULL
  features <- NULL
  if (is.data.frame(data)) {
    if (kind == "features") {
      features <- data
    } else if (kind == "arc" && is.null(arcs)) {
      arcs <- data
    } else if (!value %in% names(data)) {
      stop("the track data frame has no '", value, "' column, got ",
           paste(names(data), collapse = ", "))
    } else {
      values <- data[[value]]
    }
  } else if (is.numeric(data)) {
    values <- data
  } else if (is.character(data) && length(data) == 1) {
    text <- data
  } else if (!is.null(data)) {
    stop("geom_msa_track takes a data frame, a numeric vector or a string, ",
         "got a '", class(data)[1], "'")
  }

  encoding <- if (kind == "features") {
    drop_null(list(
      fill = if (!is.null(fill)) {
        c(list(field = fill), color_members(
          range, domain, title, type, domain_min, domain_max, domain_mid,
          domain_quantile, reverse, labels, missing_label
        ))
      },
      text = if (!is.null(label)) list(field = label)
    ))
  } else {
    bar <- kind == "bar"
    color_scale <- color_members(
      range, domain, title, type,
      if (!bar) domain_min, if (!bar) domain_max, domain_mid,
      if (!bar) domain_quantile, reverse, labels, missing_label
    )
    drop_null(list(
      color = if (length(color_scale) > 0) color_scale,
      y = if (bar) value_scale(domain_min, domain_max, domain_quantile, log)
    ))
  }
  track <- drop_null(list(
    id = id %||% slug(name), name = name, kind = kind, values = values,
    data = text, arcs = arcs, features = features, row = row, color = color,
    height = height, position = position,
    encoding = if (length(encoding) > 0) encoding
  ))
  msa_layer(append = list(columnTracks = convert_column_tracks(list(track))))
}

slug <- function(name) {
  gsub("^-+|-+$", "", gsub("[^a-z0-9]+", "-", tolower(name)))
}

#' Domain annotations over the alignment
#'
#' Draws each annotation as a box across the columns it covers, in a color per
#' feature type, with a legend. Takes a GFF3 file, a URL the viewer fetches,
#' GFF3 text, or a data frame of features with \code{row} (or
#' \code{seqname}), \code{start} and \code{end} columns, whose other
#' columns are fields an encoding can read. See the \code{gff} argument of
#' \code{\link{msaview}}.
#'
#' \code{react-msaview-cli interpro} builds the file from InterPro's
#' precomputed matches.
#'
#' @param gff The annotations.
#' @param gff_url A URL the viewer fetches the GFF from, the same as passing
#'   the URL as \code{gff}.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "p53.aln") + geom_msa_domains("p53-domains.gff")
#' }
#' @export
geom_msa_domains <- function(gff = NULL, gff_url = NULL) {
  if (!is.null(gff_url)) {
    return(msa_layer(set = list(gffUrl = check_url(gff_url))))
  }
  if (is.null(gff)) {
    stop("geom_msa_domains needs a gff or a gff_url")
  }
  if (is_url(gff)) {
    return(msa_layer(set = list(gffFilehandle = uri_location(gff))))
  }
  if (is.data.frame(gff)) {
    return(msa_layer(set = list(features = convert_features(gff))))
  }
  msa_layer(set = list(gff = convert_gff(gff)))
}

#' The residue color scale, the channel it paints and the field it reads
#'
#' Writes the viewer's \code{residueColor}, the color object of the residue
#' mark. Each argument left unset keeps what \code{msaview(residue_color =)}
#' or an earlier \code{scale_residue_color()} gave it, so a layer can change
#' the channel alone and keep the scheme.
#'
#' @param scheme A color scheme name, such as \code{"clustal"} or
#'   \code{"nucleotide"}, or a named vector of colors keyed by residue letter,
#'   such as \code{c(K = "#1f77b4", R = "#1f77b4")}, which leaves every letter
#'   it does not name uncolored. \code{msaview()}'s \code{residue_color} lists
#'   the schemes.
#' @param channel What the scale paints: \code{"fill"} (default) colors each
#'   cell's background, \code{"color"} colors the letters, which lets a domain
#'   overlay mark each span with a bar under the row.
#' @param field What a cell reads: \code{"feature"} (default) takes the color
#'   a domain overlay gives the span covering the cell, while the overlay draws
#'   one, and \code{"residue"} keeps the residue's color and marks each domain
#'   with a bar under its row.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "p53.aln") +
#'   geom_msa_domains("p53-domains.gff") +
#'   scale_residue_color("clustal", channel = "color")
#'
#' msaview(msa = "p53.aln") +
#'   geom_msa_domains("p53-domains.gff") +
#'   scale_residue_color("clustal", field = "residue")
#'
#' msaview(msa = "p53.aln") +
#'   scale_residue_color(c(K = "#1f77b4", R = "#1f77b4", D = "#d62728"))
#' }
#' @export
scale_residue_color <- function(scheme = NULL, channel = NULL, field = NULL) {
  color <- convert_residue_color(list(
    channel = channel, field = field, scale = scheme
  ))
  msa_layer(into = if (!is.null(color)) {
    list(residueColor = list(set = color))
  } else {
    list()
  })
}

#' The viewer's chrome and cell size
#'
#' @param theme \code{"light"}, \code{"dark"}, or a list of MUI theme options
#'   merged over the JBrowse theme.
#' @param col_width,row_height Pixels per column and per row.
#' @param draw_tree Draw the phylogeny. FALSE leaves a gutter holding the row
#'   labels alone.
#' @param show_branch_len Draw branch lengths. FALSE draws a cladogram.
#' @param tree_area_width Width of the tree and label gutter in pixels. Unset,
#'   the gutter fits the labels and the tree.
#' @param allowed_gappyness Hide columns that are at least this percent gaps.
#' @param hide_header Leave out the viewer's toolbar.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "p53.aln") + theme_msa("dark", col_width = 14, hide_header = TRUE)
#' }
#' @export
theme_msa <- function(theme = NULL, col_width = NULL, row_height = NULL,
                      draw_tree = NULL, show_branch_len = NULL,
                      tree_area_width = NULL, allowed_gappyness = NULL,
                      hide_header = NULL) {
  msa_layer(set = list(
    theme = theme, colWidth = col_width, rowHeight = row_height,
    drawTree = draw_tree, showBranchLen = show_branch_len,
    treeAreaWidth = tree_area_width,
    allowedGappyness = allowed_gappyness, hideHeader = hide_header
  ))
}

#' The order and root of the tree
#'
#' \code{order = "ladderize"} is ggtree's \code{ladderize}, and
#' \code{root} does what treeio's \code{root} and phytools'
#' \code{midpoint.root} do to the tree before ggtree draws it. The rows of
#' the alignment follow the tips.
#'
#' @param order \code{"branchLength"}, \code{"input"}, \code{"ladderize"}
#'   or \code{"ladderizeReverse"}. See \code{tree_order} in
#'   \code{\link{msaview}}.
#' @param root \code{"midpoint"}, or the outgroup's tip names.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   coord_tree(order = "ladderize", root = "midpoint")
#' }
#' @export
coord_tree <- function(order = NULL, root = NULL) {
  msa_layer(set = list(
    treeOrder = check_tree_order(order),
    treeRoot = convert_tree_root(root)
  ))
}

#' The span the viewer opens on
#'
#' Zooms and scrolls to a span once the alignment loads. \code{start} and
#' \code{end} alone count alignment columns; with \code{row} they count that
#' row's own residues.
#'
#' @param start,end The span, 1-based and inclusive.
#' @param row A row name, making the span residues of that row.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "p53.aln") + coord_msa(100, 288, row = "Human")
#' }
#' @export
coord_msa <- function(start, end, row = NULL) {
  msa_layer(set = list(
    region = convert_region(drop_null(list(start = start, end = end, row = row)))
  ))
}

#' A block of the alignment, selected
#'
#' Selects alignment columns \code{start} to \code{end} across the rows
#' \code{rows} names, or across every row without it, and draws the block as
#' a blue band. A shift-drag in the viewer replaces it, and in Shiny
#' \code{input$<output id>_selection} holds the block the reader selected.
#'
#' \code{trim = TRUE} hides every column outside the selection, the way
#' \code{allowed_gappyness} hides gappy ones. The tracks, the highlights and
#' \code{\link{msaview_plot}} draw the remaining columns, and the column
#' numbers stay the file's. The trim leaves every row in place;
#' \code{\link{msaview_save}(selection_only = TRUE)} draws the selected block
#' alone, from the first selected row to the last.
#'
#' @param start,end The columns, 1-based and inclusive.
#' @param rows A character vector of row names, or NULL for every row.
#' @param trim Hide every column outside the selection.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "p53.aln") +
#'   geom_msa_selection(100, 120, rows = c("Human", "Mouse"))
#'
#' msaview(msa = "p53.aln") + geom_msa_selection(94, 292, trim = TRUE)
#' }
#' @export
geom_msa_selection <- function(start, end, rows = NULL, trim = FALSE) {
  msa_layer(set = list(
    selection = convert_selection(
      drop_null(list(start = start, end = end, rows = rows))
    ),
    trimToSelection = if (isTRUE(trim)) TRUE
  ))
}

#' Structures the rows map onto
#'
#' Says which residue of which structure each residue of a row is, for a page
#' showing a structure beside the alignment. Takes the list of mappings
#' \code{msaview()} takes.
#'
#' @param data A list of mappings, each a list with \code{row},
#'   \code{structure} and \code{segments}.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' mapping <- list(
#'   row = "Human",
#'   structure = list(id = "1TUP", kind = "experimental", asymId = "A"),
#'   segments = list(list(rowStart = 94, rowEnd = 289,
#'                        structStart = 94, structEnd = 289))
#' )
#' msaview(msa = "p53.aln") + geom_msa_structure(list(mapping))
#' }
#' @export
geom_msa_structure <- function(data) {
  msa_layer(append = list(residueMappings = convert_residue_mappings(data)))
}

#' Extra fields per row
#'
#' Carries a table of fields per alignment row, which the channels in
#' \code{\link{scale_row_color}} read. A data frame gives one row per alignment
#' row, named by its \code{label} or \code{row} column, which is the shape a
#' ggtree \code{tibble(label = , trait = )} already has.
#'
#' @param data A data frame of fields per row, or a named list of them.
#' @param key The column naming each row. Taken from \code{label}, \code{row} or
#'   \code{name} when absent.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' lineages <- data.frame(
#'   label = c("A/duck/Anhui/1/2013", "A/chicken/Taiwan/a174/2015"),
#'   clade = c("2.3.4.4b", "2.3.2.1c")
#' )
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_rowdata(lineages) +
#'   scale_row_color("clade")
#' }
#' @export
geom_msa_rowdata <- function(data, key = NULL) {
  msa_layer(set = list(rowData = convert_row_data(data, key)))
}

#' Color a mark by a field of the row or feature table
#'
#' Colors one of the marks the viewer always draws by a field of a table.
#' \code{"tipLabel"} colors each tip label in the tree, \code{"rowTint"}
#' washes the row across the tree gutter and the alignment, and
#' \code{"branch"} colors a tree edge whose tips all share one value, which is
#' ggtree's \code{groupClade} with the group read from the table, and
#' \code{"rowLabel"} names the field each tip label prints, which is ggtree's
#' \code{geom_tiplab(aes(label = ))}; all four read
#' \code{\link{geom_msa_rowdata}}'s table, and a row with no value in the
#' \code{rowLabel} field prints its row name. \code{"nodePoint"} fills the
#' bubble on each internal node of the tree from its field \code{support},
#' the label a Newick file writes after a node's closing paren, which is
#' ggtree's \code{geom_nodepoint(aes(color = support))}; a node with no label
#' keeps the plain bubble. \code{"featureFill"} colors
#' each span of the annotation overlay and \code{"featureLabel"} names the
#' field drawn inside a span, both from the features the GFF carries, where a
#' field is an annotation property (\code{accession}, \code{name},
#' \code{featureType}) or a GFF attribute such as \code{Name} or
#' \code{gene}. A feature carrying a GFF \code{color=} attribute keeps that
#' color.
#'
#' \code{range} is a palette name, written to the scale's \code{scheme}, or a
#' vector of colors, written to its \code{range}. \code{domain} lists values
#' in legend order, each taking the \code{range} color at its index. A vector
#' of colors colors as many values as it holds, so a value a \code{domain} and
#' \code{range} of one length leave out keeps the plain mark. Records
#' declaring the same scale and \code{title} share one scale and one legend.
#' \code{"featureLabel"} and \code{"rowLabel"} read only \code{field}.
#'
#' \code{type = "linear"} or \code{"log"} reads each value as a number and
#' colors it from a ramp between \code{domain_min} and \code{domain_max}, and
#' its legend is a gradient bar, which \code{\link{msaview_plot}} draws as a
#' \code{scale_fill_gradientn()}. An end left unset follows the field's
#' values, or with \code{domain_quantile} below 1 the quantile at each end.
#' \code{type = "threshold"} cuts the number line at the values \code{domain}
#' lists and colors each interval from \code{range}, which is one longer.
#' \code{type = "identity"} reads each value as the CSS color it names, for a
#' table whose producer chose the colors; its legend lists the colors
#' \code{domain} names, under the names \code{labels} gives them.
#' The viewer reports a contradiction among these, such as a
#' \code{domain_max} below \code{domain_min}, in the legend, and
#' \code{msaview_plot()} prints it as a message.
#'
#' @param field The field to read.
#' @param channel \code{"tipLabel"} (default), \code{"rowTint"},
#'   \code{"branch"}, \code{"rowLabel"}, \code{"nodePoint"},
#'   \code{"featureFill"} or \code{"featureLabel"}.
#' @param range A palette name, or a vector of colors. A named vector, such as
#'   \code{c(H5 = "#e41a1c")}, colors the values it names, its names read as
#'   \code{domain}. The palettes of an ordinal scale are \code{"ggplot"},
#'   \code{"set1"}, \code{"dark2"}, \code{"okabeito"} and \code{"tableau"};
#'   the ramps of the others are \code{"viridis"} (the default),
#'   \code{"magma"}, \code{"inferno"}, \code{"cividis"}, \code{"reds"},
#'   \code{"blues"}, \code{"redblue"} and \code{"purpleorange"}. Two or more
#'   colors make a ramp of their own.
#' @param domain Field values in legend order, a threshold scale's cuts in
#'   ascending order, or the colors an identity scale's legend lists.
#' @param title The legend's title. Default: the field.
#' @param type \code{"ordinal"} (the default), \code{"linear"}, \code{"log"},
#'   \code{"threshold"} or \code{"identity"}.
#' @param domain_min,domain_max The low and the high end of a ramp. Unset, the
#'   field's lowest and highest values.
#' @param domain_mid The value a ramp's middle color marks, for a diverging
#'   ramp such as \code{"redblue"}.
#' @param domain_quantile The quantile an unset ramp end follows, from 0.5 to
#'   1. Default 1, the extremes.
#' @param reverse \code{TRUE} turns the ramp round.
#' @param labels The legend's name for each \code{domain} value of an ordinal
#'   or identity scale, or for each interval of a threshold scale, in order.
#' @param missing_label The legend's name for the grey a record with no value,
#'   or one that is no number, draws in.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_rowdata(lineages) +
#'   scale_row_color("clade", range = "set1") +
#'   scale_row_color("clade", channel = "rowTint",
#'                   range = c("2.3.4.4b" = "#e41a1c"))
#'
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_rowdata(lineages) +
#'   scale_row_color("year", type = "linear", range = "viridis",
#'                   domain_min = 2000, domain_max = 2024)
#'
#' msaview(msa = "genes.aln", gff = "genes.gff") +
#'   scale_row_color("gene", channel = "featureFill", range = "set1") +
#'   scale_row_color("Name", channel = "featureLabel")
#'
#' msaview(msa = "model.aln", gff = "plddt.gff") +
#'   scale_row_color("plddt", channel = "featureFill", type = "threshold",
#'                   domain = c(50, 70, 90),
#'                   range = c("#ff7d45", "#ffdb13", "#65cbf3", "#0053d6"),
#'                   labels = c("Very low", "Low", "Confident", "Very high"))
#'
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_rowdata(lineages) +
#'   scale_row_color("strain", channel = "rowLabel") +
#'   scale_row_color("support", channel = "nodePoint", type = "threshold",
#'                   domain = c(70, 95),
#'                   range = c("#d7191c", "#fdae61", "#1a9641"))
#'
#' msaview(msa = "brca1.aln", gff = "clinvar.gff") +
#'   scale_row_color("color", channel = "featureFill", type = "identity",
#'                   domain = c("#d7191c", "#2c7bb6"),
#'                   labels = c("pathogenic", "benign"), title = "ClinVar")
#' }
#' @export
scale_row_color <- function(field, channel = "tipLabel", range = NULL,
                            domain = NULL, title = NULL, type = NULL,
                            domain_min = NULL, domain_max = NULL,
                            domain_mid = NULL, domain_quantile = NULL,
                            reverse = NULL, labels = NULL,
                            missing_label = NULL) {
  encoding <- c(
    list(channel = channel, field = field),
    if (!channel %in% text_channels) {
      color_members(range, domain, title, type, domain_min, domain_max,
                    domain_mid, domain_quantile, reverse, labels,
                    missing_label)
    }
  )
  msa_layer(append = list(encodings = convert_encodings(list(encoding))))
}

#' Color each node point by its support
#'
#' Fills the bubble on each internal node of the tree from the node's support
#' value, the label a Newick file writes after its closing paren: \code{95}
#' and \code{80} in \code{((A,B)95,(C,D)80);}. This is ggtree's
#' \code{geom_nodepoint(aes(color = support))}, and
#' \code{\link{scale_row_color}(channel = "nodePoint")} with the field
#' \code{support}. The viewer reads each label as a string, so a threshold
#' scale bands bootstrap percentages and an ordinal scale colors named clades.
#' A node with no label keeps the plain bubble. \code{\link{msaview_plot}}
#' draws a point on each node the scale colors.
#'
#' @param field The node field to read. The viewer holds one, \code{support}.
#' @inheritParams scale_row_color
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln", tree = "h5_bootstrap.nh") +
#'   geom_msa_nodepoint(type = "threshold", domain = c(70, 95),
#'                      range = c("#d7191c", "#fdae61", "#1a9641"),
#'                      labels = c("< 70", "70-95", ">= 95"),
#'                      title = "bootstrap")
#' }
#' @export
geom_msa_nodepoint <- function(field = "support", range = NULL, domain = NULL,
                               title = NULL, type = NULL, domain_min = NULL,
                               domain_max = NULL, domain_mid = NULL,
                               domain_quantile = NULL, reverse = NULL,
                               labels = NULL, missing_label = NULL) {
  scale_row_color(field, "nodePoint", range, domain, title, type, domain_min,
                  domain_max, domain_mid, domain_quantile, reverse, labels,
                  missing_label)
}

#' Print a field of the row table as each tip label
#'
#' Prints the \code{\link{geom_msa_rowdata}} field \code{field} as each tip's
#' label in the tree, so a tree keyed by accession shows strain names. A row
#' with no value in the field prints its row name. This is ggtree's
#' \code{geom_tiplab(aes(label = field))}, and
#' \code{\link{scale_row_color}(field, channel = "rowLabel")}. The tree area
#' fits the labels it prints, and \code{\link{msaview_plot}} names the rows
#' with them.
#'
#' @param field The field to print.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_rowdata(data.frame(label = accessions, strain = strains)) +
#'   geom_msa_tiplab("strain")
#' }
#' @export
geom_msa_tiplab <- function(field) {
  scale_row_color(field, "rowLabel")
}

# A layer function's scale arguments as the color object: `scale` and
# `legend`, spelt as Vega-Lite and jbrowse-components spell them. A named
# `range`, ggplot2's `values = c(a = "red")`, gives the names as `domain`. One
# unnamed string that is not a color is a palette or ramp name, written to
# `scale$scheme`; colors go to `scale$range`, an array at any length.
color_members <- function(range = NULL, domain = NULL, title = NULL,
                          type = NULL, domain_min = NULL, domain_max = NULL,
                          domain_mid = NULL, domain_quantile = NULL,
                          reverse = NULL, labels = NULL,
                          missing_label = NULL) {
  check_scale_type(type)
  if (!is.null(names(range))) {
    domain <- domain %||% names(range)
    range <- unname(as.character(unlist(range)))
  }
  scheme <- is.character(range) && length(range) == 1 &&
    !grepl("^#|\\(", range)
  scale <- drop_null(list(
    type = type,
    domain = if (!is.null(domain)) I(scale_domain(domain)),
    range = if (!is.null(range) && !scheme) I(as.character(range)),
    scheme = if (scheme) range,
    domainMin = domain_min, domainMax = domain_max, domainMid = domain_mid,
    domainQuantile = domain_quantile, reverse = reverse
  ))
  legend <- drop_null(list(
    title = title,
    labels = if (!is.null(labels)) I(as.character(labels)),
    missingLabel = missing_label
  ))
  drop_null(list(
    scale = if (length(scale) > 0) scale,
    legend = if (length(legend) > 0) legend
  ))
}

# A bar's or a point's value scale as `{scale: {type, domainMin, domainMax,
# domainQuantile}}`, NULL when no argument sets a member
value_scale <- function(domain_min = NULL, domain_max = NULL,
                        domain_quantile = NULL, log = NULL) {
  scale <- drop_null(list(
    type = if (isTRUE(log)) "log",
    domainMin = domain_min, domainMax = domain_max,
    domainQuantile = domain_quantile
  ))
  if (length(scale) > 0) list(scale = scale)
}

#' A column of colored cells beside the tree
#'
#' Draws one cell per alignment row between the tree and the alignment, colored
#' by a field of \code{\link{geom_msa_rowdata}}'s table, which is ggtree's
#' \code{gheatmap}. Eight of these make the tip-aligned matrix a surveillance
#' figure puts beside its phylogeny.
#'
#' The scale takes the arguments \code{\link{scale_row_color}} takes, so a
#' numeric field such as a collection year draws on a ramp with
#' \code{type = "linear"}. Strips declaring the same scale and \code{title},
#' such as eight segment columns over one set of lineage colors, share one
#' scale over every value their fields take, and one legend.
#'
#' @param field The field to color by.
#' @param range A palette name, or a vector of colors. A named vector, such as
#'   \code{c(H5 = "#e41a1c")}, colors the values it names, its names read as
#'   \code{domain}.
#' @param width The column's width in pixels. Default: the row height.
#' @param header The name drawn above the column. Default: the field.
#' @inheritParams scale_row_color
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_rowdata(segments) +
#'   geom_msa_strip("HA", range = "set1", width = 12) +
#'   geom_msa_strip("NA", header = "NA segment") +
#'   geom_msa_strip("year", type = "linear", range = "viridis")
#' }
#' @export
geom_msa_strip <- function(field, range = NULL, width = NULL,
                           header = NULL, title = NULL, domain = NULL,
                           type = NULL, domain_min = NULL, domain_max = NULL,
                           domain_mid = NULL, domain_quantile = NULL,
                           reverse = NULL, labels = NULL,
                           missing_label = NULL) {
  panel <- c(
    list(
      kind = "strip",
      encoding = list(
        fill = c(list(field = field), color_members(
          range, domain, title, type, domain_min, domain_max, domain_mid,
          domain_quantile, reverse, labels, missing_label
        ))
      )
    ),
    drop_null(list(width = width, header = header))
  )
  msa_layer(append = list(rowPanels = convert_row_panels(list(panel))))
}

#' A column of gene arrows beside the tree
#'
#' Draws the spans a GFF carries in a panel of their own, one row per
#' alignment row, as arrows where a feature has a strand. \code{x = "column"}
#' draws them in the alignment's columns; \code{x = "position"} draws each row
#' in its own residue positions on one linear scale, so a genome with no
#' alignment has an x. With \code{align} the rows shift so that the named gene
#' starts at one x, which is gggenes' \code{make_alignment_dummies}, and a row
#' carrying no such gene keeps its own origin.
#'
#' \code{fill} and \code{label} name fields of the features: an
#' \code{Annotation} field (\code{name}, \code{accession},
#' \code{featureType}) or a GFF attribute of column 9, such as \code{Name}.
#' \code{fill} colors each arrow's body, as ggplot2 and gggenes'
#' \code{fill} aesthetic does, and \code{label} is drawn inside it. Without
#' \code{fill} the spans take the colors the alignment's overlay gives them.
#' The \code{fill} scale takes the arguments \code{\link{scale_row_color}}
#' takes.
#'
#' @param x \code{"position"} (default) or \code{"column"}.
#' @param fill The feature field to fill the spans by.
#' @param label The feature field drawn inside each span.
#' @param range A palette name, or a vector of colors, for \code{fill}. A
#'   named vector colors the values it names, its names read as
#'   \code{domain}.
#' @param align A gene name to line the rows up on.
#' @param width The column's width in pixels. Default: 200.
#' @param header The name drawn above the column.
#' @param domain Values of \code{fill} in legend order, or a threshold
#'   scale's cuts.
#' @param title The legend's title. Default: the \code{fill} field.
#' @param position How overlapping features share a row's height:
#'   \code{"identity"} (the default) packs them into lanes, and
#'   \code{"strandpile"} puts forward features above reverse ones.
#' @inheritParams scale_row_color
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(tree = "marker.nh", gff = "neighborhoods.gff") +
#'   geom_msa_features(fill = "Name", label = "Name", align = "genE",
#'                     width = 320, header = "neighborhood")
#' }
#' @export
geom_msa_features <- function(x = "position", fill = NULL, label = NULL,
                              range = NULL, align = NULL, width = NULL,
                              header = NULL, domain = NULL, title = NULL,
                              position = NULL, type = NULL,
                              domain_min = NULL, domain_max = NULL,
                              domain_mid = NULL, domain_quantile = NULL,
                              reverse = NULL, labels = NULL,
                              missing_label = NULL) {
  encoding <- drop_null(list(
    fill = if (!is.null(fill)) {
      c(list(field = fill), color_members(
        range, domain, title, type, domain_min, domain_max, domain_mid,
        domain_quantile, reverse, labels, missing_label
      ))
    },
    text = if (!is.null(label)) list(field = label)
  ))
  panel <- drop_null(list(
    kind = "features", x = x,
    encoding = if (length(encoding) > 0) encoding,
    transform = if (!is.null(align)) list(list(type = "align", on = align)),
    width = width, header = header, position = position
  ))
  msa_layer(append = list(rowPanels = convert_row_panels(list(panel))))
}

#' A bar or a dot per row beside the tree
#'
#' Draws a numeric \code{\link{geom_msa_rowdata}} field in a panel between the
#' tree and the alignment, a bar per row from 0 to its value, or a dot at its
#' value, which is ggtree's \code{geom_facet} with \code{geom_col} or
#' \code{geom_point}. The value runs along the panel's width, and an axis
#' above the panel, titled by \code{header}, labels it. A row with no number in
#' the field draws no mark.
#'
#' Each end of the value scale is \code{domain_min} or \code{domain_max}
#' where given, and otherwise follows the field's values, clipped at
#' \code{domain_quantile}. A bar panel with an end unset takes 0 into its
#' scale, so its bars start at 0. \code{log = TRUE} spaces the values by
#' orders of magnitude, reads only the positive ones, and draws each bar from
#' the scale's low end.
#'
#' \code{fill} on a bar panel, or \code{color} on a point panel, names a field
#' each mark is colored by: the same field, or another, such as a host. Its
#' scale takes the arguments \code{\link{scale_row_color}} takes, apart from
#' the ends, which set the value scale here, so a color ramp follows the
#' field's values. Panels and encodings declaring the same scale and title list
#' one legend. Without a color field the marks draw grey.
#'
#' @param x The numeric field to draw.
#' @param fill,color The field to color each bar or dot by.
#' @param width The panel's width in pixels. Default: 80.
#' @param header The axis title drawn above the panel. Default: \code{x}.
#' @param domain_min,domain_max The ends of the value scale. Unset, they
#'   follow the field's values.
#' @param domain_quantile The quantile an unset end of the value scale
#'   follows, from 0.5 to 1. Default 1, the extremes.
#' @param log \code{TRUE} draws the values on a log scale.
#' @param range The colors or palette name of the color scale, as
#'   \code{\link{scale_row_color}} takes them.
#' @param domain Values of the color field in legend order, or a threshold
#'   scale's cuts.
#' @param title The color legend's title. Default: the color field.
#' @param size A dot's diameter in pixels, at most the row height. Default: 6.
#' @inheritParams scale_row_color
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' isolates <- data.frame(label = c("A", "B", "C"),
#'                        genome = c(29903, 30119, 29751),
#'                        year = c(2019, 2020, 2021),
#'                        host = c("bat", "human", "pangolin"))
#' msaview(msa = "cov.aln", tree = "cov.nh") +
#'   geom_msa_rowdata(isolates) +
#'   geom_msa_bar("genome", fill = "host", range = "dark2",
#'                header = "genome (bp)") +
#'   geom_msa_point("year", size = 5)
#' }
#' @export
geom_msa_bar <- function(x, fill = NULL, width = NULL, header = NULL,
                         domain_min = NULL, domain_max = NULL,
                         domain_quantile = NULL, log = NULL, range = NULL,
                         domain = NULL, title = NULL, type = NULL,
                         domain_mid = NULL, reverse = NULL, labels = NULL,
                         missing_label = NULL) {
  value_panel(
    "bar", x, fill,
    value_scale(domain_min, domain_max, domain_quantile, log),
    color_members(range, domain, title, type, domain_mid = domain_mid,
                  reverse = reverse, labels = labels,
                  missing_label = missing_label),
    width = width, header = header
  )
}

#' @rdname geom_msa_bar
#' @export
geom_msa_point <- function(x, color = NULL, size = NULL, width = NULL,
                           header = NULL, domain_min = NULL,
                           domain_max = NULL, domain_quantile = NULL,
                           log = NULL, range = NULL, domain = NULL,
                           title = NULL, type = NULL, domain_mid = NULL,
                           reverse = NULL, labels = NULL,
                           missing_label = NULL) {
  value_panel(
    "point", x, color,
    value_scale(domain_min, domain_max, domain_quantile, log),
    color_members(range, domain, title, type, domain_mid = domain_mid,
                  reverse = reverse, labels = labels,
                  missing_label = missing_label),
    width = width, header = header, size = size
  )
}

#' A field printed beside each row
#'
#' Prints a \code{\link{geom_msa_rowdata}} field in a panel between the tree
#' and the alignment, cut in the middle where it runs longer than the panel,
#' which is ggtree's \code{geom_facet} with \code{geom_text}. \code{color}
#' names a field each text is colored by, through a scale that takes the
#' arguments \code{\link{scale_row_color}} takes; under
#' \code{type = "identity"} each text draws in the color its field names.
#' Without \code{color} the text draws in the theme's text color.
#'
#' @param field The field to print.
#' @param color The field to color each text by.
#' @param width The panel's width in pixels. Default: 80.
#' @param header The name drawn above the panel. Default: \code{field}.
#' @param range The colors or palette name of the color scale, as
#'   \code{\link{scale_row_color}} takes them.
#' @param domain Values of \code{color} in legend order, or a threshold
#'   scale's cuts.
#' @param title The color legend's title. Default: the \code{color} field.
#' @inheritParams scale_row_color
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln", tree = "h5.nh") +
#'   geom_msa_rowdata(lineages) +
#'   geom_msa_text("clade", color = "clade", range = "set1", width = 70)
#' }
#' @export
geom_msa_text <- function(field, color = NULL, width = NULL, header = NULL,
                          range = NULL, domain = NULL, title = NULL,
                          type = NULL, domain_min = NULL, domain_max = NULL,
                          domain_mid = NULL, domain_quantile = NULL,
                          reverse = NULL, labels = NULL,
                          missing_label = NULL) {
  value_panel(
    "text", field, color, NULL,
    color_members(range, domain, title, type, domain_min, domain_max,
                  domain_mid, domain_quantile, reverse, labels,
                  missing_label),
    width = width, header = header
  )
}

# A bar, point or text panel: the field it reads, as `encoding.x` with its
# value scale or as `encoding.text`, and `encoding.color` over `color_field`
# where one is named
value_panel <- function(kind, field, color_field, scale, color_scale,
                        width = NULL, header = NULL, size = NULL) {
  if (is.null(field)) {
    stop("a ", kind, " panel needs a field to draw")
  }
  if (is.null(color_field) && length(color_scale) > 0) {
    stop("a ", kind, " panel's color scale needs a field to color by, as ",
         if (kind == "bar") "fill" else "color", " = ")
  }
  read <- c(list(field = field), scale)
  encoding <- drop_null(list(
    x = if (kind != "text") read,
    text = if (kind == "text") read,
    color = if (!is.null(color_field)) c(list(field = color_field), color_scale)
  ))
  panel <- drop_null(list(kind = kind, encoding = encoding, width = width,
                          header = header, size = size))
  msa_layer(append = list(rowPanels = convert_row_panels(list(panel))))
}

#' Draw each row as its differences from one row
#'
#' Every other row draws as its differences from \code{row}, with matching
#' residues as \code{.}, the reading aid comparative papers use.
#'
#' @param row The row to compare against.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "p53.aln") + stat_msa_diff("Human")
#' }
#' @export
stat_msa_diff <- function(row) {
  msa_layer(set = list(relativeTo = sanitize_names_or_null(row)))
}

#' Keep the rows whose field passes a test
#'
#' Keeps the rows whose \code{\link{geom_msa_rowdata}} field takes one of
#' \code{values}, or lies in \code{range}, which is dplyr's \code{filter}
#' over the row table. Each call adds a test, and a row stays when it passes
#' every one. A row with no value for the field fails. With a tree, the
#' viewer prunes the tree to the tips it keeps.
#'
#' A filter keeping no row is not applied, and the viewer's header warns.
#'
#' @param field The field to test.
#' @param values The values to keep, compared as strings, so \code{2015}
#'   matches \code{"2015"}.
#' @param range The lowest and highest value to keep, inclusive. \code{NA}
#'   leaves an end open, as \code{c(2010, NA)}.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln") +
#'   geom_msa_rowdata(isolates) +
#'   filter_msa_rows("host", values = c("duck", "goose")) +
#'   filter_msa_rows("year", range = c(2010, NA))
#' }
#' @export
filter_msa_rows <- function(field, values = NULL, range = NULL) {
  test <- drop_null(list(field = field, oneOf = values, range = range))
  msa_layer(into = list(rowTransform = list(
    append = list(filter = list(convert_row_filter(test)))
  )))
}

#' Order the rows by their fields
#'
#' Orders the rows by \code{\link{geom_msa_rowdata}} fields, the first
#' field first and each further one ordering the rows the ones before it leave
#' tied, which is dplyr's \code{arrange}. Wrap a field in \code{desc()} to
#' sort it descending, or give \code{order}. Values sort in natural order:
#' numbers by value, and digit runs by magnitude, so \code{2.3.2.1a} sorts
#' before \code{2.3.4.4b}. Rows with no value sort last.
#'
#' A tree sets the row order, so with a tree the viewer leaves the rows in
#' the tree's order and its header warns. A later \code{arrange_msa_rows()}
#' replaces an earlier one.
#'
#' @param ... Field names, each a string, a character vector of them, or
#'   \code{desc(field)}. \code{desc()} is read from the call, so dplyr's need
#'   not be attached.
#' @param order \code{"ascending"} or \code{"descending"}, one for every
#'   field or one per field. It replaces what \code{desc()} says.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "h5.aln") +
#'   geom_msa_rowdata(isolates) +
#'   arrange_msa_rows("host", desc("year"))
#'
#' msaview(msa = "h5.aln") +
#'   geom_msa_rowdata(isolates) +
#'   arrange_msa_rows("year", order = "descending")
#' }
#' @export
arrange_msa_rows <- function(..., order = NULL) {
  env <- parent.frame()
  keys <- unlist(lapply(as.list(substitute(list(...)))[-1], function(e) {
    descending <- is.call(e) &&
      utils::tail(as.character(e[[1]]), 1) == "desc"
    fields <- eval(if (descending) e[[2]] else e, env)
    if (!is.character(fields)) {
      stop("arrange_msa_rows takes field names as strings, got ",
           deparse(e), call. = FALSE)
    }
    lapply(fields, function(field) {
      list(field = field, order = if (descending) "descending")
    })
  }), recursive = FALSE)
  if (length(keys) == 0) {
    stop("arrange_msa_rows needs a field to sort by")
  }
  if (!is.null(order)) {
    if (!length(order) %in% c(1, length(keys))) {
      stop("arrange_msa_rows takes one order, or one per field (",
           length(keys), "), got ", length(order))
    }
    order <- rep_len(order, length(keys))
    keys <- Map(function(key, o) c(key["field"], order = o), keys, order)
  }
  msa_layer(into = list(rowTransform = list(
    set = list(sort = unname(lapply(keys, function(key) {
      convert_row_sort(drop_null(key))
    })))
  )))
}

#' Band the rows by a field
#'
#' Stacks the rows in bands by a \code{\link{geom_msa_rowdata}} field, each
#' labeled with its value in the gutter the \code{"bracket"} clade mark
#' draws in, which is the row axis of ggplot2's \code{facet_grid}. The values
#' \code{sort} lists come first, then the others in natural order, and the
#' rows with no value last, under no label. \code{\link{msaview_plot}}
#' draws each band's bracket and label.
#'
#' With a tree the rows stay in the tree's order, and a band is labeled only
#' where its rows are one clade of the tree.
#'
#' @param field The field to band by.
#' @param sort Values whose bands come first, in that order.
#' @return A layer to add to a viewer with \code{+}.
#'
#' @examples
#' \dontrun{
#' msaview(msa = "kinases.aln") +
#'   geom_msa_rowdata(groups) +
#'   facet_msa_rows("group", sort = c("TK", "TKL"))
#' }
#' @export
facet_msa_rows <- function(field, sort = NULL) {
  msa_layer(into = list(rowTransform = list(
    set = list(facet = convert_row_facet(drop_null(list(field = field,
                                                         sort = sort))))
  )))
}
