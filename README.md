R/multiparental
==

An `R` package for exploring, visualising and analysing haplotype data from "multiparental populations", ie. populations arising from breeding schemes such as the [Collaborative Cross](http://csbio.unc.edu/CCstatus/index.py), [Diversity Outbred](http://jaxmice.jax.org/strain/009376.html) or [Heterogeneous Stock](http://mus.well.ox.ac.uk/mouse/HS/) in which the founder haplotypes are known.  Currently it is mouse-centric but not mouse-specific; it could be used with eg. the [MAGIC population](http://mus.well.ox.ac.uk/19genomes/magic.html) in *Arabidopsis*.

The core data structure is built on the [`GenomicRanges::GRanges` class](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html), which supports fast overlap queries and other operations on genomic intervals (more generally, integer ranges.)

Due to unresolved issues around function exports from `GenomicRanges`, `multiparental` only works if `GenomicRanges` has also been loaded into the `R` environment.

Functions
--
* fast queries across large sets (>1000) of inferred haplotypes
* `ggplot2`-based plotting of genome mosaics
* haplotype-based QTL-mapping under several phenotype models

Interop
--
The package interacts with Karl Broman's `R/qtl` (http://www.rqtl.org) package for QTL-mapping, and my `megamuga_hmm` package (coming soon) for performing haplotype inference from dense genotypes.

TODO
--
* add case-control QTL-mapping (via 2xk contingency tables)
* break off visualisation functions into separate package (`genomeviz`)
* figure out dependency on `GenomicRanges`
