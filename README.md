# miscgene
An R package that contains some helper function for accessing biomaRt and Reactome. The main functions currently include:

- symbol_to_entrez(); uses biomaRt to obtain entrezgene IDs for a set of Gene Symbols.
- entrez_to_symbol(); obtains Gene Symbols for a set of entrezgene IDs.
- entrez_to_attribute(); obtains user-specified attributes for a set of entrezgene IDs.
- get_reactome_pathways(); uses reactome.db to obtain a list of Reactome pathways. 

# Installation
miscgene can be installed from GitHub using the `devtools` package:

``` r
# install.packages('devtools')
devtools::install_github('tgrimes/miscgene')
```