
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scrub <img src="tools/logo.png" align="right" style="padding-left:10px;background-color:white;" width="150px" />

<!-- badges: start -->

[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)
[![R build
status](https://github.com/IDEELResearch/scrub/workflows/R-CMD-check/badge.svg)](https://github.com/IDEELResearch/scrub/actions)
<!-- badges: end -->

## {scrub} - Spatial Cleaned Resistance and Unified Base

This is a working R compendium (think R package but for reproducible
analysis). A good overview on research compendiums, see the [R for
Reproducible Research](https://annakrystalli.me/rrresearch/index.html)
course.

### Installation

    git clone https://github.com/IDEELResearch/scrub.git
    cd scrub
    open scrub.Rproj

`devtools::install_dev_deps()` will install all required packages, as
specified in the Imports in DESCRIPTION. (At a later date when analysis
is finalised, `renv` can be used to create a reproducible R environment
that anyone can use by calling `renv::restore` to set up package
dependencies.)

### Overview

The structure within analysis is as follows:

    R/                            # Packaged R functions 

    analysis/
        |
        ├── 01_xxxxx /            # analysis scripts used for generating figures
        |
        ├── plots/                # location of figures produced by the analysis scripts
        |
        ├── tables/               # location of any tables produced by the analysis scripts
        |
        ├── data_raw/             # data obtained from elsewhere and treated read-only    
        |
        ├── data_derived/         # intermediate data generated during the analysis

- Analysis scripts are to be run in the numbered order they are
  included. If there are shared numbers, then any order of these scripts
  works.

- Data that is **read only**, e.g. data shared from elsewhere and not
  generated using code in this repository, is stored in `data_raw`

- Data that is generated using code in this repository is stored in
  `data_derived`.

- Outputs from the analysis scripts, such as plots and tables are stored
  in `plots` and `tables` respectively.

### Milestones (This can probably be deleted once all in Github issues)

- V1 of scrub csv table that is ready to go into stave
- Finish systematic lit review abstract screening for inclusion
- Finish full text extraction of lit review screening  
- Checking and cross validating the extracted seekdeep server data
- Pf7k, WWARN, WHO pulling:
- Meta study formatting correctly
- Cleaning and deduplication
- Write a validation suite (lat long) - jeff/cecile sits in python
- Final merged two
- Deduplication
- Writing up what we did to get this data

### Compendium DOI:

<https://zenodo.org/record/XXX>

The files at the URL above will generate the results as found in the
publication.

### The R package

This repository is organized as an R package. There are only a few R
functions exported in this package - the majority of the R code is in
the analysis directory. The R package structure is here to help manage
dependencies, to take advantage of continuous integration, and so we can
keep file and data management simple. For any R packages that are used
frequently in this repository, they are documented in `R/` and are used
in the analysis folder using `devtools::load_all()`.

To download the package source as you see it on GitHub, for offline
browsing, use this line at the shell prompt (assuming you have Git
installed on your computer):

``` r
git clone https://github.com/IDEELResearch/scrub.git
```

Once the download is complete, open the `scrub.Rproj` in RStudio to
begin working with the package and compendium files. We will endeavour
to keep all package dependencies required listed in the DESCRIPTION.

<!-- To add this once all the analysis is done -->
<!-- In addition, once analysis is completed, we will use `renv` to track package dependencies for reproducibility. Please use `renv::restore` to restore the state of the project and see https://rstudio.github.io/renv/articles/renv.html for more information. -->

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2024, copyright
holder: OJ Watson

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse
