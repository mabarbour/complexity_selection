[![DOI](https://zenodo.org/badge/108833263.svg)](https://zenodo.org/badge/latestdoi/108833263)

# complexity_selection
This repository contains data and code to reproduce analyses reported in the manuscript "Loss of consumers constrains phenotypic evolution in the resulting food web" (in press at *Evolution Letters*). The original version of this manuscript was posted on [*bioRxiv*](https://www.biorxiv.org/content/10.1101/527416v1).

For a pretty HTML document displaying all code to reproduce these results, download *reproduce_analyses.html* from the **analyses** folder. For an ugly version that is viewable on GitHub, click on the "reproduce_analyses.md" file.

# Repository structure

**manuscript**: contains *manuscript.Rmd* that knits together *manuscript.pdf*. It also contains a style file for *The American Naturalist* (*amnat.bst*), a bibtex file of references (*references.bib*), and an intermediate *.tex* file. *Am-Nat-preamble-latex.tex* contains useful code for generating a pdf output that resembles the Latex template for *The American Naturalist* in Rmarkdown.

- **prior_versions**: contains prior versions of this manuscript.

**analyses**: contains code to reproduce all analyses. It also contains the figures that are embedded within *manuscript.pdf*. The main workhorse is *reproduce_analyses.Rmd*, which is called on when knitting the manuscript together with *manuscript.Rmd*.

- **reproduce_analyses_cache** and **reproduce_analyses_files**: contain all of the components for generating *reproduce_analyses.html*. This information was stored because it takes greater than 1 h to re-run all of the analyses from scratch.

**data**: includes the dataset analyzed for this manuscript. 

- **raw_data_management**: contains raw data and code used to generate the final dataset. If content of this folder is not clear after reading *finalize_manuscript.dataset.R*, *gall_contents.info.md*, and *revising_dataset.R*, then please open an issue in this repository to notify me. 

The names of all other files should be self-explanatory. If they are not, please open an issue in this repository.
