# Example data

The gallery and the examples read these files. Each is a copy of a file in
JBrowseMSA's `packages/examples/data`, where the scripts that made it live.

| File | What it is | Source |
| --- | --- | --- |
| `insulin.aln`, `insulin.nh` | Preproinsulin across nine vertebrates, aligned with ClustalW, and its neighbor-joining tree | UniProt, by JBrowseMSA's `scripts/examples-gen` |
| `kinase.aln`, `kinase.nh`, `kinase-domains.gff` | Ten Src-family kinases, their tree and their InterPro domains | UniProt and InterPro, by `scripts/examples-gen` |
| `trna.stock` | A 24-sequence subset of the tRNA seed alignment, with its consensus secondary structure and tree | Rfam RF00005 |
| `gene-cluster.stock`, `gene-cluster.gff` | A synthetic colinear gene cluster across six genomes, for illustration | JBrowseMSA's `scripts/gene-cluster/generate.mjs` |
