import './polyfills.ts'

import { colord } from 'colord'

import {
  arcColorValue,
  letterValue,
  numberText,
} from '@lib/columnTrackColor.ts'
import {
  axisFontSize,
  axisTicks,
  valuePosition,
} from '@lib/components/axis.ts'
import { featureSpanMarks } from '@lib/components/msa/drawFeatureSpans.ts'
import { domainOverlaySource } from '@lib/components/msa/renderBoxFeatureCanvasBlock.ts'
import { tileColorFn } from '@lib/components/msa/tileColor.ts'
import { featurePanelSource } from '@lib/components/rowpanels/renderFeaturePanel.ts'
import { barTrackValues } from '@lib/components/tracks/drawTracks.ts'
import { opaqueColor } from '@lib/components/tracks/highlightBands.ts'
import {
  bracketBarWidth,
  bracketGap,
  cladeLabelBoxes,
} from '@lib/components/tree/cladeBrackets.ts'
import { ellipsize } from '@lib/components/tree/ellipsize.ts'
import { featureField } from '@lib/featureFields.ts'
import { layoutTree } from '@lib/hierarchy.ts'
import { CHAR_WIDTH_RATIO } from '@lib/measureTextCanvas.ts'
import MSAModelF from '@lib/model.ts'
import { fieldChannel } from '@lib/scales.ts'
import { columnLogoStack } from '@lib/sequenceLogo.ts'
import { nodeField } from '@lib/treeDisplay.ts'
import { adjustColorForContrast } from '@lib/util.ts'
import { viewerSnapshot } from '@lib/viewerSnapshot.ts'

import type { MSAViewerProps } from '@lib/components/MSAViewer.tsx'
import type {
  SpanBand,
  SpanSource,
} from '@lib/components/msa/drawFeatureSpans.ts'
import type { LaidOutNode } from '@lib/hierarchy.ts'
import type { LegendItem, LegendSection } from '@lib/legendSpec.ts'
import type {
  Annotation,
  BasicTrack,
  ColumnSpan,
  NodeWithIds,
  PanelAxis,
  ResolvedStripPanel,
  ResolvedTextPanel,
  ResolvedValuePanel,
} from '@lib/types.ts'

type Table = Record<string, unknown[]>

function table(...names: string[]): Table {
  return Object.fromEntries(names.map(name => [name, []]))
}

function push(t: Table, row: Record<string, unknown>) {
  for (const key in t) {
    t[key]!.push(row[key] ?? null)
  }
}

// R reads hex, with alpha, but not rgb() or hsl()
function hex(color: string | null | undefined) {
  return color == null ? undefined : colord(color).toHex()
}

function luminance(color: string) {
  const { r, g, b } = colord(color).toRgb()
  const channel = (v: number) => {
    const c = v / 255
    return c <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4
  }
  return 0.2126 * channel(r) + 0.7152 * channel(g) + 0.0722 * channel(b)
}

// MUI's getContrastText at its default threshold of 3, as the viewer's
// letters take it
function contrastText(fill: string) {
  return 1.05 / (luminance(fill) + 0.05) >= 3 ? '#ffffff' : '#000000de'
}

type Model = ReturnType<ReturnType<typeof MSAModelF>['create']>

/**
 * A column track's marks as R draws them. Each bar, strip cell, arc and text
 * cell carries the color the viewer paints it in, and where a scale over the
 * track's own values gave that color, the value it read as `key`, which the
 * legend titled `title` lists. A mark whose value the legend hides is left
 * out, as the viewer leaves it undrawn.
 */
function trackTable(model: Model, track: BasicTrack['model']) {
  const { kind } = track
  const colors = model.columnTrackColors.get(track.id)
  const scaled = colors?.scale !== undefined
  const title = colors?.scale?.title
  if (kind === 'bar') {
    const raw = model.columnTrackContent.get(track.id)?.raw
    const bars = table('col', 'value', 'key', 'color')
    barTrackValues(model, track.id).forEach((value, i) => {
      if (!colors?.hiddenColumns?.[i]) {
        const own = colors?.columns?.[i]
        push(bars, {
          col: i + 1,
          value,
          key: scaled && own !== undefined ? numberText(raw?.[i]) : undefined,
          color: hex(own ?? track.barColor ?? 'gray'),
        })
      }
    })
    return { title, bars }
  }
  if (kind === 'strip') {
    const raw = model.columnTrackContent.get(track.id)?.raw
    const cells = table('col', 'key', 'color')
    colors?.columns?.forEach((color, i) => {
      if (color !== undefined) {
        push(cells, {
          col: i + 1,
          key: scaled ? numberText(raw?.[i]) : undefined,
          color: hex(color),
        })
      }
    })
    return { title, cells }
  }
  if (kind === 'arc') {
    const arcs = table('start', 'end', 'key', 'color')
    for (const arc of track.arcs ?? []) {
      const value = arcColorValue(arc, !!colors?.identity)
      if (!colors?.isHidden(value)) {
        const fromScale =
          arc.color === undefined && scaled ? colors.colorOf(value) : undefined
        push(arcs, {
          start: arc.start + 1,
          end: arc.end + 1,
          key: fromScale === undefined ? undefined : value,
          color: hex(arc.color ?? fromScale ?? track.arcColor),
        })
      }
    }
    return { title, arcs }
  }
  if (kind === 'logo') {
    const { colStats, alphabetMaxBits: maxBits, colorScheme } = model
    const letters = table('col', 'letter', 'ymin', 'ymax', 'color')
    for (let col = 0; col < colStats.numColumns; col++) {
      let bottom = 0
      for (const { letter, bits } of columnLogoStack(colStats, col, maxBits)) {
        push(letters, {
          col: col + 1,
          letter,
          ymin: bottom,
          ymax: bottom + bits,
          color: hex(colorScheme[letter]),
        })
        bottom += bits
      }
    }
    return { maxBits, letters }
  }
  if (kind === 'highlights') {
    const bands = table('start', 'end', 'lane', 'label', 'color')
    for (const band of track.bands ?? []) {
      push(bands, {
        start: band.startCol + 1,
        end: band.endCol + 1,
        lane: band.lane,
        label: band.label,
        color: band.color ? hex(opaqueColor(band.color)) : undefined,
      })
    }
    return { bands }
  }
  if (kind === 'features') {
    return featureTrackTables(model, track)
  }
  if (kind === 'text') {
    const scheme = colors?.letters ?? model.colorScheme
    const cells = table('col', 'letter', 'key', 'fill', 'text')
    const data = track.data ?? ''
    for (let col = 0; col < data.length; col++) {
      const letter = data[col]!
      const value = letterValue(letter)
      if (value !== undefined && !colors?.isHidden(value)) {
        const fill = model.bgColor ? hex(scheme[value]) : undefined
        push(cells, {
          col: col + 1,
          letter,
          key: scaled && fill ? value : undefined,
          fill,
          text: fill ? contrastText(fill) : '#000000de',
        })
      }
    }
    return { title, cells }
  }
  return undefined
}

// QuickJS has no fonts, so a label fits by the export's 0.6em fallback
const measureText = (text: string, fontSize: number) =>
  text.length * fontSize * CHAR_WIDTH_RATIO

// the source puts a row's top edge at y - rowHeight, its row center rowHeight
// below the top of the figure
const rowsDown = (model: Model) => (y: number) =>
  (y + model.rowHeight / 2) / model.rowHeight

/**
 * The span mark as polygons and labels, in the source's pixels. `yOf` maps a
 * pixel y to the caller's: rows down from the first row's top edge for the
 * overlay and a panel, pixels down for a track. `x` stays in pixels for the
 * caller to scale. `key` is the value the fill scale reads, which a legend
 * item names as its `value`: the field's value, or with no field the
 * accession the domain palette keys on.
 */
function spanTables<T extends SpanBand>(
  model: Model,
  source: SpanSource<T>,
  field: string | undefined,
  title: string,
  yOf = rowsDown(model),
) {
  const { rowHeight } = model
  const keyOf = (a: Annotation) =>
    field === undefined ? a.accession : featureField(a, field)
  const shapes = table('group', 'x', 'y', 'key', 'fill')
  const labels = table('x', 'y', 'label', 'color', 'size')
  let group = 0
  const rows = featureSpanMarks({
    ...source,
    measureText,
    xMin: -Infinity,
    xMax: Infinity,
  })
  for (const row of rows) {
    for (const { band, fill, box, outline } of row.shapes) {
      group++
      const points: [number, number][] = box
        ? [
            [box.x, box.y],
            [box.x + box.w, box.y],
            [box.x + box.w, box.y + box.h],
            [box.x, box.y + box.h],
          ]
        : outline!
      for (const [x, y] of points) {
        push(shapes, {
          group,
          x,
          y: yOf(y),
          key: keyOf(band.annotation),
          fill: hex(fill),
        })
      }
    }
    for (const { text, x, y, fill } of row.labels) {
      push(labels, {
        x,
        y: yOf(y),
        label: text,
        color: contrastText(hex(fill)!),
        size: row.fontSize / rowHeight,
      })
    }
  }
  return { title, shapes, labels }
}

// drawFeatureTrack's layout in drawTracks.ts: a pixel clear above and below
// each lane, and a head as long as the lane is tall up to 16 pixels
const laneGap = 2
const maxHeadLength = 16

/** a features track's spans as the span mark draws them, y in track pixels */
function featureTrackTables(model: Model, track: BasicTrack['model']) {
  const scale = model.featureTrackScales.get(track.id)
  const spec = model.columnTracks.find(t => t.id === track.id)
  const { colWidth } = model
  const height = track.height
  const spans = track.spans ?? []
  return scale
    ? spanTables<ColumnSpan>(
        model,
        {
          rows: spans.length > 0 ? [{ y: height, bands: spans }] : [],
          layout: {
            height: laneCount => Math.max(1, height / laneCount - laneGap),
            top: (_y, lane, h) => lane * (h + laneGap) + laneGap / 2,
            headLength: h => Math.min(h, maxHeadLength),
          },
          xOf: span => [span.startCol * colWidth, span.endCol * colWidth],
          colors: scale.colors,
          labelOf: span => scale.labels?.get(span.annotation),
        },
        fieldChannel(spec?.encoding?.fill)?.field,
        scale.scale?.title ?? track.name,
        y => y,
      )
    : undefined
}

// drawAxis's spacing in components/axis.ts: the ticks tried first, and the
// room kept between two labels once they are measured
const firstSpacing = 20
const labelPadding = 8

/**
 * A panel axis as drawAxis draws it, in the panel's pixels: each tick's
 * position and its label, spaced so the widest label fits between two
 */
function axisGuide({ scale, extent, format }: PanelAxis) {
  const length = Math.abs(extent[1] - extent[0])
  const at = valuePosition(scale, extent)
  const labelsOf = ({ values, step }: { values: number[]; step: number }) =>
    values.map(value => ({
      value,
      text: format(value, scale.type === 'log' ? value : step),
    }))
  const first = labelsOf(axisTicks(scale, length, firstSpacing))
  const room =
    Math.max(0, ...first.map(l => measureText(l.text, axisFontSize))) +
    labelPadding
  const ticks =
    room > firstSpacing ? labelsOf(axisTicks(scale, length, room)) : first
  return {
    at: ticks.map(t => at(t.value)),
    label: ticks.map(t => t.text),
  }
}

/**
 * A bar or point panel's marks in its pixels, as renderValuePanel draws them:
 * a bar from the axis position of 0, or the scale's low end where 0 lies
 * outside it, to its value, and a dot at its value. `key` is the value the
 * color field holds, which the panel's legend lists.
 */
function valuePanelMarks(model: Model, panel: ResolvedValuePanel) {
  const { scale, extent } = panel.axis
  const at = valuePosition(scale, extent)
  const base = scale.type === 'log' ? extent[0] : at(0)
  const marks = table('row', 'x', 'base', 'key', 'color')
  model.leaves.forEach((leaf, i) => {
    const { name } = leaf.data
    const value = panel.values.get(name)
    if (value !== undefined && !panel.hiddenRows.has(name)) {
      push(marks, {
        row: i + 1,
        x: at(value),
        base,
        key: panel.colorField && model.rowData[name]?.[panel.colorField],
        color: hex(panel.colors.get(name) ?? 'gray'),
      })
    }
  })
  return {
    title: panel.legendTitle ?? panel.colorField,
    size: Math.min(panel.size, model.rowHeight),
    marks,
    axis: axisGuide(panel.axis),
  }
}

// renderTextPanel's padding either side of the text
const textPadding = 2

/** a text panel's text per row, cut in the middle to the panel's width */
function textPanelMarks(model: Model, panel: ResolvedTextPanel) {
  const { fontSize } = model
  const texts = table('row', 'text', 'key', 'color')
  model.leaves.forEach((leaf, i) => {
    const { name } = leaf.data
    const value = panel.texts.get(name)
    if (value !== undefined && !panel.hiddenRows.has(name)) {
      push(texts, {
        row: i + 1,
        text: ellipsize(
          t => measureText(t, fontSize),
          value,
          panel.width - 2 * textPadding,
        ).text,
        key: panel.colorField && model.rowData[name]?.[panel.colorField],
        color: hex(panel.colors.get(name)),
      })
    }
  })
  return {
    title: panel.legendTitle ?? panel.colorField,
    padding: textPadding,
    texts,
  }
}

function stripCells(model: Model, panel: ResolvedStripPanel) {
  const cells = table('row', 'key', 'color')
  model.leaves.forEach((leaf, i) => {
    const { name } = leaf.data
    const color = panel.colors.get(name)
    if (color) {
      push(cells, {
        row: i + 1,
        key: model.rowData[name]?.[panel.field],
        color: hex(color),
      })
    }
  })
  return { title: panel.legendTitle ?? panel.field, cells }
}

function rowPanelTables(model: Model) {
  return model.resolvedRowPanels.map(panel => ({
    kind: panel.kind,
    header: panel.header,
    width: panel.width,
    ...(panel.kind === 'strip'
      ? stripCells(model, panel)
      : panel.kind === 'text'
        ? textPanelMarks(model, panel)
        : panel.kind === 'features'
          ? {
              axis: panel.axis ? axisGuide(panel.axis) : undefined,
              ...spanTables(
                model,
                featurePanelSource(model, panel, model.leaves),
                panel.field,
                panel.legendTitle ?? panel.field ?? 'Domains',
              ),
            }
          : valuePanelMarks(model, panel)),
  }))
}

// the ruler is left to the x axis, which numbers the same columns
function trackTables(model: Model) {
  const on = new Set(model.turnedOnTracks.map(t => t.model.id))
  return model.tracks
    .filter(t => t.model.kind !== 'ruler')
    .map(({ model: track }) => ({
      id: track.id,
      name: track.name,
      kind: track.kind,
      on: on.has(track.id),
      height: track.height,
      ...trackTable(model, track),
    }))
}

// the gutter at the right of the tree area, in its own px: each bracket's bar
// as renderCladeBrackets draws it and each label in the box cladeLabelBoxes
// gives it, its top in rows
function cladeGutter(model: Model) {
  const { cladeGutterWidth: width, effectiveTreeAreaWidth, rowHeight } = model
  if (width === 0) {
    return undefined
  }
  const bars = table('first', 'last', 'color')
  for (const { rows, mark, markColor } of model.resolvedClades) {
    if (mark === 'bracket') {
      push(bars, {
        first: rows[0] + 1,
        last: rows[1] + 1,
        color: hex(markColor),
      })
    }
  }
  const labels = table('label', 'left', 'top', 'height', 'horizontal', 'color')
  for (const box of cladeLabelBoxes(model)) {
    push(labels, {
      label: box.label,
      left: box.left - (effectiveTreeAreaWidth - width),
      top: box.top / rowHeight,
      height: box.height / rowHeight,
      horizontal: box.horizontal,
      color: hex(box.color),
    })
  }
  return { width, barX: bracketGap + bracketBarWidth / 2, bars, labels }
}

const rampEnd = (label: string) => Number(label.replace(/^[≤≥]/, ''))

/**
 * Each legend section as R reads it: its swatches as parallel columns, less
 * the ones whose value the reader hid, which the export leaves out too, and a
 * ramp's gradient with its ends as numbers and whether the scale reads them on
 * a log axis. A section names no scale, so the scale whose items it lists
 * gives its kind.
 */
function legendTables(model: Model) {
  const scales = [...model.resolvedScales.values()]
  const same = (a: LegendItem[], b: LegendItem[]) =>
    JSON.stringify(a) === JSON.stringify(b)
  return model.legends.map((section: LegendSection) => {
    const swatches = section.items.filter(
      item => item.color !== undefined && !item.hidden,
    )
    const ramp = section.items.find(item => item.gradient)?.gradient
    const kind = scales.find(s => same(s.items, section.items))?.kind
    return {
      id: section.id,
      title: section.title,
      note: section.note ?? null,
      value: swatches.map(item => item.value ?? item.label),
      label: swatches.map(item => item.label),
      color: swatches.map(item => hex(item.color)),
      gradient: ramp
        ? {
            offset: ramp.stops.map(stop => stop.offset),
            color: ramp.stops.map(stop => hex(stop.color)),
            minLabel: ramp.minLabel,
            maxLabel: ramp.maxLabel,
            min: rampEnd(ramp.minLabel),
            max: rampEnd(ramp.maxLabel),
            log: kind === 'log',
          }
        : null,
    }
  })
}

function treeTables(model: Model) {
  // an edge takes the `branch` channel's color of the node at its far end, a
  // node's own connector the color of the node
  const segments = table('x', 'y', 'xend', 'yend', 'color')
  const branch = model.branchColors
  const colorOf = (node: LaidOutNode<NodeWithIds>) =>
    hex(branch?.get(node.data.id))
  const collapsed = table('x', 'y', 'group')
  // each internal node the nodePoint channel colors, keyed by its support
  const nodes = table('x', 'y', 'key', 'color')
  const nodeColors = model.nodePointColors
  if (model.tree.noTree) {
    return undefined
  }
  const { root, rootToTipLength } = layoutTree(model.root, {
    rowHeight: 1,
    treeWidth: 1,
  })
  const byLength = model.showBranchLen && rootToTipLength > 0
  const xOf = (node: LaidOutNode<NodeWithIds>) =>
    byLength ? node.branchLength : node.depthFraction
  // layoutTree puts the first tip at 0.5; the tables number rows from 1
  const yOf = (node: LaidOutNode<NodeWithIds>) => node.row + 0.5
  const stack = [root]
  let group = 0
  while (stack.length > 0) {
    const node = stack.pop()!
    const x = xOf(node)
    if (node.children?.length) {
      const nodeColor = nodeColors?.get(node.data.id)
      if (nodeColor) {
        push(nodes, {
          x,
          y: yOf(node),
          key: nodeField(node.data, 'support'),
          color: hex(nodeColor),
        })
      }
      push(segments, {
        x,
        y: yOf(node.children[0]!),
        xend: x,
        yend: yOf(node.children.at(-1)!),
        color: colorOf(node),
      })
      for (const child of node.children) {
        const y = yOf(child)
        push(segments, { x, y, xend: xOf(child), yend: y, color: colorOf(child) })
        stack.push(child)
      }
    } else if (byLength && node.collapsedBranchLength !== undefined) {
      group++
      const y = yOf(node)
      for (const [px, py] of [
        [x, y],
        [node.collapsedBranchLength, y - 0.35],
        [node.collapsedBranchLength, y + 0.35],
      ]) {
        push(collapsed, { x: px, y: py, group })
      }
    }
  }
  return {
    segments,
    collapsed,
    nodes,
    nodeTitle: model.encodingFor('nodePoint')?.title ?? null,
    byLength,
    scaleBar: byLength ? (model.treeScaleBar ?? null) : null,
  }
}

function figure(propsJson: string) {
  const props = JSON.parse(propsJson) as MSAViewerProps
  const model = MSAModelF().create(viewerSnapshot(props))
  model.setWidth(1000)
  const { leaves, columns, numColumns, relativeTo } = model
  const names = leaves.map(leaf => leaf.data.name)
  const rowOf = new Map(names.map((name, i) => [name, i + 1]))
  const tints = model.rowTints
  const tintOf = new Map(
    tints ? model.rowNames.map((name, i) => [name, tints[i]] as const) : [],
  )

  const colorAt = tileColorFn(model)
  const colors = table('col', 'letter', 'fill', 'text', 'letterColor')
  const seen = new Set<string>()
  const dynamic = model.dynamicColorSchemeName !== undefined
  for (const seq of columns.values()) {
    for (let col = 0; col < seq.length; col++) {
      const letter = seq[col]!
      const key = dynamic ? `${col}:${letter}` : letter
      if (!seen.has(key)) {
        seen.add(key)
        const fill = hex(colorAt(col, letter))
        if (fill) {
          push(colors, {
            col: dynamic ? col + 1 : undefined,
            letter,
            fill,
            text: contrastText(fill),
            letterColor: adjustColorForContrast(fill, '#ffffff'),
          })
        }
      }
    }
  }

  const domains = table(
    'row',
    'start',
    'end',
    'label',
    'accession',
    'fill',
    'stroke',
    'text',
  )
  for (const [name, bands] of model.domainBands) {
    const row = rowOf.get(name)
    if (row !== undefined) {
      for (const { annotation, startCol, endCol } of bands) {
        const color = model.featureColors.get(annotation)
        push(domains, {
          row,
          start: startCol + 1,
          end: endCol,
          label: annotation.name,
          accession: annotation.accession,
          fill: hex(color?.fill),
          stroke: hex(color?.stroke),
          text: color ? contrastText(hex(color.fill)!) : undefined,
        })
      }
    }
  }

  const highlights = table('start', 'end', 'row', 'label', 'color', 'mark')
  for (const h of model.resolvedHighlights) {
    const rows = h.rowIndices.length
      ? h.rowIndices.map(r => r + 1)
      : [undefined]
    for (const row of rows) {
      push(highlights, {
        start: h.startCol === undefined ? undefined : h.startCol + 1,
        end: h.endCol === undefined ? undefined : h.endCol + 1,
        row,
        label: h.label,
        color: hex(h.color),
        mark: h.mark,
      })
    }
  }

  const clades = table('first', 'last', 'mark', 'color', 'label')
  for (const c of model.resolvedClades) {
    push(clades, {
      first: c.rows[0] + 1,
      last: c.rows[1] + 1,
      mark: c.mark,
      // a highlight fills with the wash the viewer draws, a bracket with its
      // mark's own color
      color: hex(c.mark === 'highlight' ? c.color : (c.markColor ?? c.color)),
      label: c.label,
    })
  }

  return JSON.stringify({
    numColumns,
    rows: {
      name: names,
      label: names.map(name => model.rowLabelOf(name)),
      seq: names.map(name => columns.get(name) ?? ''),
      labelColor: names.map(name => hex(model.tipLabelColors?.get(name)) ?? null),
      tint: names.map(name => hex(tintOf.get(name)) ?? null),
    },
    colors,
    reference: relativeTo ?? null,
    fillCells: model.bgColor && !model.actuallyShowDomains,
    letterMode: !model.bgColor,
    lettersOverDomains:
      model.actuallyShowDomains &&
      !model.subFeatureRows &&
      !model.domainUnderline,
    treeAreaWidth: model.effectiveTreeAreaWidth,
    cladeGutter: cladeGutter(model) ?? null,
    tree: treeTables(model) ?? null,
    domains,
    highlights,
    clades,
    legends: legendTables(model),
    scaleProblems: model.scaleProblems,
    fileColumns: Array.from(
      { length: numColumns },
      (_, col) => model.visibleColToGlobalCol(col) + 1,
    ),
    rowHeight: model.rowHeight,
    colWidth: model.colWidth,
    fontSize: model.fontSize,
    features: model.actuallyShowDomains
      ? spanTables(
          model,
          domainOverlaySource(model, model.leaves),
          model.featureFillEncoding?.field,
          model.featureFillEncoding?.title ?? 'Domains',
        )
      : null,
    rowPanels: rowPanelTables(model),
    tracks: trackTables(model),
  })
}

;(globalThis as Record<string, unknown>).msaviewFigure = figure
