#!/usr/bin/env Rscript
# The README's gallery: each example drawn twice from the package's own data,
# by JBrowseMSA's SVG export (msaview_save(), through react-msaview-cli) and by
# msaview_plot(), side by side. Writes the images, the README's table and
# docs/gallery.md, which holds each recipe. Run from the repo root.
#   Rscript scripts/gallery.R [example ...]  redraw all, or the examples named
#   Rscript scripts/gallery.R --docs         rewrite the table and docs only
#   Rscript scripts/gallery.R --check        run each recipe, then as --docs
# The export runs react-msaview-cli from a JBrowseMSA checkout where
# `pnpm build` has run, at JBROWSEMSA (default ../react-msaview), or the command
# MSAVIEW_CLI names. Needs rsvg-convert and ImageMagick for the webp.
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
jbm <- Sys.getenv("JBROWSEMSA", file.path("..", "react-msaview"))
cli <- strsplit(Sys.getenv("MSAVIEW_CLI",
  paste("node", file.path(jbm, "packages/cli/dist/index.js"))), " ")[[1]]
width <- 1200
dir <- "gallery"

# Each recipe's last value is the widget; `plot` is the msaview_plot() call
# drawing it, `autoplot(w)` where it is absent
examples <- list(
  insulin_tree = list(
    title = "Preproinsulin across nine vertebrates, beside its tree",
    caption = "The B and A chains of mature insulin stay conserved under the conservation bars, while the C-peptide between them varies.",
    r = r"-(
w <- msaview(msa = e("insulin.aln"), tree = e("insulin.nh"))
)-"),
  insulin_diff = list(
    title = "The same alignment drawn as its differences from human",
    caption = "Every residue matching the human row prints as a dot, so the substitutions stand out in the C-peptide and the signal peptide.",
    r = r"-(
w <- msaview(msa = e("insulin.aln"), tree = e("insulin.nh"),
  relative_to = "Human")
)-"),
  insulin_clades = list(
    title = "Clade brackets over the insulin tree",
    caption = "A bracket names the five mammals and another the three non-mammals, each in the gutter beside the tree.",
    r = r"-(
w <- msaview(msa = e("insulin.aln"), tree = e("insulin.nh")) +
  geom_msa_clade(range = c("Pig", "Mouse"), tips = 5, mark = "bracket",
    label = "Mammals", color = "#1565c0") +
  geom_msa_clade(range = c("Xenopus", "Chicken"), tips = 3, mark = "bracket",
    label = "Non-mammals")
)-"),
  kinase_domains = list(
    title = "Ten Src-family kinases with their InterPro domains",
    caption = "The SH3, SH2 and kinase domains line up in every member, the unique N-terminal region before them varying in length.",
    r = r"-(
w <- msaview(msa = e("kinase.aln"), tree = e("kinase.nh"),
  gff = e("kinase-domains.gff"), col_width = 1.5)
)-"),
  trna_structure = list(
    title = "tRNA with its consensus secondary structure",
    caption = "The four stems of the cloverleaf draw as base-pair arcs over the dot-bracket row, above conservation.",
    r = r"-(
w <- msaview(msa = e("trna.stock"))
)-"),
  gene_cluster = list(
    title = "A gene cluster across six genomes, with each genome's clade",
    caption = "Each genome's genes as arrows in its own coordinates, genC reversed in every one and genB missing from Genome_5, beside a strip of its clade.",
    r = r"-(
clades <- data.frame(name = paste0("Genome_", 1:6),
  clade = c("A", "A", "B", "B", "C", "C"))
w <- msaview(msa = e("gene-cluster.stock"), gff = e("gene-cluster.gff"),
  col_width = 1) +
  geom_msa_rowdata(clades) +
  geom_msa_strip("clade", header = "Clade") +
  geom_msa_features(x = "position", color = "Name", label = "Name",
    header = "Genes")
)-")
)
examples <- lapply(examples, function(x) {
  x$r <- trimws(x$r)
  x
})

png_of <- function(name, side) file.path(dir, sprintf("%s-%s.png", name, side))
webp_of <- function(name, side) sub("png$", "webp", png_of(name, side))

recipe_widget <- function(example) {
  env <- new.env()
  env$e <- function(f) system.file("extdata", f, package = "ggmsaview")
  eval(parse(text = example$r), env)
  env$w
}

# a fifth of an inch a row and a little over half an inch a track, which
# ggplot2's text sizes read well at
figure_height <- function(w) {
  lay <- msaview_layout(w)
  keys <- unique(unlist(c(lay$features$shapes$key,
    lapply(lay$row_panels, function(p) c(p$cells$key, p$shapes$key)))))
  max(0.9 + 0.2 * nrow(lay$rows) + 0.6 * sum(vapply(lay$tracks,
    function(t) t$on, logical(1))), 0.6 + 0.3 * length(keys))
}

draw_viewer <- function(name, w) {
  svg <- tempfile(fileext = ".svg")
  msaview_save(w, svg, width = width, cli = cli)
  status <- system2("rsvg-convert", c("-w", width, "-o",
    shQuote(png_of(name, "viewer")), shQuote(svg)))
  if (status != 0) stop(name, ": rsvg-convert failed")
}

draw_r <- function(name, w) {
  ggsave(png_of(name, "r"), autoplot(w), width = width / 120,
    height = figure_height(w), dpi = 120, bg = "white")
}

to_webp <- function(name) {
  for (side in c("viewer", "r")) {
    system2("magick", c(png_of(name, side), "-quality", 85,
      webp_of(name, side)))
    file.remove(png_of(name, side))
  }
}

image_pair <- function(name, e, prefix) sprintf(
  "<tr><td width=\"50%%\"><img src=\"%s%s\" alt=\"%s in JBrowseMSA. %s\" width=\"100%%\"></td><td width=\"50%%\"><img src=\"%s%s\" alt=\"%s drawn by ggmsaview. %s\" width=\"100%%\"></td></tr>",
  prefix, webp_of(name, "viewer"), e$title, e$caption, prefix,
  webp_of(name, "r"), e$title, e$caption)

gallery_table <- function(examples) {
  row <- function(name, e) c(
    sprintf("<tr><td colspan=\"2\"><b><a href=\"docs/gallery.md#%s\">%s</a></b></td></tr>",
      name, e$title),
    image_pair(name, e, ""))
  c("<!-- gallery -->", "<table>",
    "<tr><th width=\"50%\">JBrowseMSA</th><th width=\"50%\">ggmsaview</th></tr>",
    unlist(Map(row, names(examples), examples)), "</table>", "<!-- /gallery -->")
}

update_readme <- function(examples) {
  lines <- readLines("README.md")
  from <- grep("^<!-- gallery -->$", lines)
  to <- grep("^<!-- /gallery -->$", lines)
  stopifnot(length(from) == 1, length(to) == 1)
  writeLines(c(lines[seq_len(from - 1)], gallery_table(examples),
    lines[-seq_len(to)]), "README.md")
}

recipe <- function(name, e) c(
  sprintf("<a id=\"%s\"></a>", name), "",
  paste("##", e$title), "",
  "<table>", image_pair(name, e, "../"), "</table>", "",
  e$caption, "",
  "```r", e$r, "autoplot(w)", "msaview_save(w, \"figure.svg\")", "```", "")

write_recipes <- function(examples) writeLines(c(
  "# Gallery recipes", "",
  "Each figure in the README's gallery, drawn from the package's own data by",
  "JBrowseMSA's SVG export and by ggmsaview.", "",
  "- Each recipe runs as written after",
  "",
  "  ```r",
  "  library(ggmsaview)",
  "  e <- function(f) system.file(\"extdata\", f, package = \"ggmsaview\")",
  "  ```",
  "",
  "- `autoplot(w)` draws the right-hand figure with ggplot2, and",
  "  `msaview_save(w, \"figure.svg\")` the left-hand one, through",
  "  react-msaview-cli, which needs Node.",
  "- `scripts/gallery.R` draws both and writes this page, so edit an example",
  "  there.",
  "",
  unlist(Map(recipe, names(examples), examples))),
  file.path("docs", "gallery.md"))

wanted <- commandArgs(trailingOnly = TRUE)
if (identical(wanted, "--check")) {
  for (name in names(examples)) {
    message(name)
    grDevices::pdf(NULL)
    print(autoplot(recipe_widget(examples[[name]])))
    grDevices::dev.off()
  }
} else if (!identical(wanted, "--docs")) {
  dir.create(dir, showWarnings = FALSE)
  todo <- if (length(wanted)) examples[wanted] else examples
  for (name in names(todo)) {
    message(name)
    w <- recipe_widget(todo[[name]])
    draw_viewer(name, w)
    draw_r(name, w)
    to_webp(name)
  }
}
dir.create("docs", showWarnings = FALSE)
update_readme(examples)
write_recipes(examples)
