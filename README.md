## A generic hidden Markov model for multi-parent populations

[![DOI](https://zenodo.org/badge/380058622.svg)](https://zenodo.org/badge/latestdoi/380058622)

This repository contains the source for the paper

> Broman KW (2022) A generic hidden Markov model
> for multi-parent populations. [G3
> (Bethesda)](https://academic.oup.com/g3journal) 12:jkab396
> [![PubMed](https://kbroman.org/icons16/pubmed-icon.png)](https://pubmed.ncbi.nlm.nih.gov/34791211/)
> [![pdf](https://kbroman.org/icons16/pdf-icon.png)](https://academic.oup.com/g3journal/article-pdf/12/2/jkab396/42382435/jkab396.pdf)
> [![GitHub](https://kbroman.org/icons16/github-icon.png)](https://github.com/kbroman/Paper_GenericHMM)
> [![doi](https://kbroman.org/icons16/doi-icon.png)](https://doi.org/10.1093/g3journal/jkab396)

- [`gen_hmm.Rnw`](gen_hmm.Rnw) - source document (LaTeX + knitr)
- [`gen_hmm_bib`](gen_hmm.bib) - BibTex bibliography
- [`DOapp/`](DOapp) - application to Diversity Outbred mice - [`do_analysis.html`](https://kbroman.org/Paper_GenericHMM/DOapp/do_analysis.html)
- [`CCapp/`](CCapp) - application to Collaborative Cross mice - [`cc_analysis.html`](https://kbroman.org/Paper_GenericHMM/CCapp/cc_analysis.html)
- [`Makefile`](Makefile) - GNU Make file

### Required R packages

- [R/qtl2](https://kbroman.org/qtl2)
- [qtl2convert](https://github.com/kbroman/qtl2convert)
- [qtl2fst](https:/github.com/kbroman/qtl2fst)
- [broman](https://github.com/kbroman/broman)
- [data.table](https://rdatatable.gitlab.io/data.table/)
- [readxl](https://readxl.tidyverse.org)
- [here](https://here.r-lib.org)

### License

The content in this repository is licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

[![CC BY 4.0](https://licensebuttons.net/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/)
