#!/usr/bin/env node
/**
 * Builds the JavaScript ggmsaview ships, from a JBrowseMSA checkout
 * (../react-msaview, or $JBROWSEMSA) where `pnpm install && pnpm build` has run:
 *
 * - inst/htmlwidgets/lib/react-msaview.umd.js, the viewer the widget mounts
 * - inst/js/msaview-core.js, js/bridge.ts and the viewer's model, which R runs
 *   in QuickJS to lay out a ggplot2 figure
 * - inst/msaviewer-props.json, the MSAViewer props test-props.R holds
 *   msaview() to
 *
 * Each file names the react-msaview version and the commit it was built from,
 * and the bundles are committed, so install_github() needs no build step.
 */
import { execFileSync } from 'node:child_process'
import fs from 'node:fs'
import { createRequire } from 'node:module'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..')
const jbm = path.resolve(
  process.env.JBROWSEMSA ?? path.join(root, '..', 'react-msaview'),
)
const lib = path.join(jbm, 'packages/lib')
const umd = path.join(lib, 'bundle/index.js')
if (!fs.existsSync(umd)) {
  console.error(
    `${umd} not found: run \`pnpm install && pnpm build\` in ${jbm}, or set JBROWSEMSA`,
  )
  process.exit(1)
}

const { version } = JSON.parse(
  fs.readFileSync(path.join(lib, 'package.json'), 'utf8'),
)
const commit = execFileSync('git', ['-C', jbm, 'rev-parse', 'HEAD'], {
  encoding: 'utf8',
}).trim()
const banner = `/* react-msaview ${version}, GMOD/JBrowseMSA ${commit}, built by scripts/build-js.mjs: do not edit */\n`

const widget = path.join(root, 'inst/htmlwidgets/lib/react-msaview.umd.js')
fs.mkdirSync(path.dirname(widget), { recursive: true })
fs.writeFileSync(widget, banner + fs.readFileSync(umd, 'utf8'))

// htmlwidgets keys its dependency cache on this version, so it moves with the
// bundle or a page serves a stale copy alongside a new widget
const yaml = path.join(root, 'inst/htmlwidgets/msaview.yaml')
fs.writeFileSync(
  yaml,
  fs
    .readFileSync(yaml, 'utf8')
    .replace(/^(\s*version:\s*)\S+/m, `$1${version}`),
)

// the bridge imports the viewer's source as `@lib/...`; the model's few
// imports from the @jbrowse/core/util barrel resolve to a shim, since the
// barrel pulls in React and MUI, and renderToSvg, which the model imports
// lazily and an iife would inline with react-dom, to a stub
const js = path.join(root, 'js')
const shims = {
  '@jbrowse/core/util': path.join(js, 'coreUtilShim.ts'),
  '@jbrowse/core/util/io': path.join(js, 'coreUtilShim.ts'),
}
const esbuild = createRequire(path.join(lib, 'package.json'))('esbuild')
await esbuild.build({
  entryPoints: [path.join(js, 'bridge.ts')],
  outfile: path.join(root, 'inst/js/msaview-core.js'),
  bundle: true,
  nodePaths: [path.join(lib, 'node_modules')],
  format: 'iife',
  target: 'es2020',
  minify: true,
  legalComments: 'none',
  banner: { js: banner.trimEnd() },
  define: { 'process.env.NODE_ENV': '"production"' },
  logLevel: 'warning',
  plugins: [
    {
      name: 'ggmsaview-core',
      setup(build) {
        build.onResolve({ filter: /^@lib\// }, args => ({
          path: path.join(lib, 'src', args.path.slice('@lib/'.length)),
        }))
        build.onResolve({ filter: /^@jbrowse\/core\/util(\/io)?$/ }, args => ({
          path: shims[args.path],
        }))
        build.onResolve({ filter: /renderToSvg\.tsx$/ }, () => ({
          path: path.join(js, 'renderToSvgStub.ts'),
        }))
      },
    },
  ],
})

const viewer = fs.readFileSync(
  path.join(lib, 'src/components/MSAViewer.tsx'),
  'utf8',
)
const body = /export interface MSAViewerProps \{([\s\S]*?)\n\}/.exec(viewer)
if (!body) {
  throw new Error('no MSAViewerProps interface in MSAViewer.tsx')
}
const props = [...body[1].matchAll(/^ {2}(\w+)\??:/gm)].map(m => m[1])
fs.writeFileSync(
  path.join(root, 'inst/msaviewer-props.json'),
  `${JSON.stringify({ version, commit, props }, null, 2)}\n`,
)

console.log(`built react-msaview ${version} from ${commit.slice(0, 10)}`)
