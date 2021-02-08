
<!-- README.md is generated from README.Rmd. Please edit that file -->

# herpesdrg

<!-- badges: start -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
![GitHub repo
size](https://img.shields.io/github/repo-size/ojcharles/herpesdrg.svg)
![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/ojcharles/herpesdrg.svg)
<!-- badges: end -->

## Overview

herpesdrg is a R package to enable antiviral drug resistance genotyping,
with Herpes Simplex Virus 1, HSV 2, Human Cytomegalovirus and Varicella
zoster virus sequencing data. Accepted inputs are FASTA (whole genomes &
fragments) which will be mapped to RefSeq NC\_001806.2. NGS variant data
assembled to NC\_001806.2 is accepted in VCF \>= ver4.0 & Varscan2 tab
formats.

#### Database

All data extracted from literature sources and contains the
relationships between:

  - Virus
  - Gene
  - Mutation
  - Drug susceptibility fold change from sensitive strain or “Resistant”
    / “Polymorphism” as in publication
  - The method for Resistance Phenotyping and strain generation
  - Publication reference

#### Web service

A user-friendly Shiny Applications has been bundled with this package.
The same application is available over the internet here
<http://cmv-resistance.ucl.ac.uk/herpesdrg/> where the terms of use are
contained.

## Installation

You can install the current version from
[GitHub](https://github.com/ojcharles/herpesdrg) with:

``` r
# install.packages("devtools")
devtools::install_github("ojcharles/herpesdrg")
```

Dependencies for FASTA file handling are MAFFT and SNP-Sites available
preferably via conda. snp-sites \>= 2.3 has been tested.

``` bash
conda config --add channels bioconda
conda install snp-sites
conda install mafft
```

## Usage

### vcf data

``` r
library("herpesdrg")

## call resistant variants
my_sample = system.file("testdata", "HCMV_A10.vcf", package = "herpesdrg")

data = call_resistance(infile = my_sample, all_mutations = F, virus = "HCMV") # options are HSV1 HSV2 HCMV and VZV

data[ , c("change", "freq", "Ganciclovir", "Cidofovir", "Foscarnet", "Letermovir", "Maribavir", "ref_link")]
#>       change   freq Ganciclovir Cidofovir Foscarnet Letermovir Maribavir
#> 1 UL54_D588N   5.2%           2       1.3       2.8                   NA
#> 2 UL54_D588N   5.2%         3.8       2.7         9                   NA
#> 3 UL97_C592G 10.77%           3                                       NA
#> 4 UL97_C592G 10.77%           3                                       NA
#> 5 UL97_C592G 10.77%           3                                       NA
#> 6 UL97_H411Y  1.69%         0.5                                       12
#> 7 UL97_H411Y  1.69%                                                   18
#> 8 UL97_T409M 37.13%                                                   90
#> 9 UL97_T409M 37.13%         0.9                                       81
#>                                                                         ref_link
#> 1                          https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059355/
#> 2 https://www.sciencedirect.com/science/article/pii/S1386653201001603?via%3Dihub
#> 3                          https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2876423/
#> 4                                           https://aac.asm.org/content/55/1/382
#> 5 https://www.sciencedirect.com/science/article/pii/S0166354219304115?via%3Dihub
#> 6                                     https://jvi.asm.org/content/82/1/246.short
#> 7 https://www.sciencedirect.com/science/article/pii/S0166354219304115?via%3Dihub
#> 8 https://www.sciencedirect.com/science/article/pii/S0166354219304115?via%3Dihub
#> 9                           https://academic.oup.com/jid/article/196/1/91/844651


## call all variants
mutations_all = call_resistance(infile = my_sample, all_mutations = T, virus = "HCMV")

#to view all mutations in resistance genes we can filter
mutations_res = mutations_all[mutations_all$GENEID %in% c("UL54", "UL97", "UL27", "UL51", "UL56", "UL89"),]

head(mutations_res[,c(1,8,21,32:40)],)
#>           change   freq    CONSEQUENCE Ganciclovir Aciclovir Cidofovir
#> 996    UL51_T47M 99.88%  nonsynonymous        <NA>      <NA>      <NA>
#> 997   UL51_Y112Y   100%     synonymous        <NA>      <NA>      <NA>
#> 999     UL54_709   1.3% not translated        <NA>      <NA>      <NA>
#> 1000    UL54_883 21.15% not translated        <NA>      <NA>      <NA>
#> 1001    UL54_885 91.02% not translated        <NA>      <NA>      <NA>
#> 1002 UL54_A1108T 99.05%  nonsynonymous        <NA>      <NA>      <NA>
#>      Foscarnet Brincidofovir Letermovir Brivudine Pencyclovir Tomeglovir
#> 996       <NA>            NA       <NA>      <NA>        <NA>         NA
#> 997       <NA>            NA       <NA>      <NA>        <NA>         NA
#> 999       <NA>            NA       <NA>      <NA>        <NA>         NA
#> 1000      <NA>            NA       <NA>      <NA>        <NA>         NA
#> 1001      <NA>            NA       <NA>      <NA>        <NA>         NA
#> 1002      <NA>            NA       <NA>      <NA>        <NA>         NA


## call nonsynonymous variants
# are there any non-synonymous (DNA variants that result in a change of amino acid) variants in resistance genes
mutations_res_nonsyn = mutations_all[mutations_all$GENEID == "UL54" & mutations_all$CONSEQUENCE == "nonsynonymous",]


# here the 6 mutations are nonsynonymous, with no identified resistance effect.
head(mutations_res_nonsyn[,c("change", "freq", "CONSEQUENCE", "Ganciclovir", "Foscarnet", "ref_title")],)
#>           change   freq   CONSEQUENCE Ganciclovir Foscarnet
#> 1002 UL54_A1108T 99.05% nonsynonymous        <NA>      <NA>
#> 1004  UL54_A692V   100% nonsynonymous        <NA>      <NA>
#> 1007  UL54_D588N   5.2% nonsynonymous           2       2.8
#> 1008  UL54_D588N   5.2% nonsynonymous         3.8         9
#> 1018  UL54_L897S 99.82% nonsynonymous        <NA>      <NA>
#> 1019 UL54_N1116K  1.16% nonsynonymous        <NA>      <NA>
#>                                                                                                                               ref_title
#> 1002                                                                                                                               <NA>
#> 1004                                                                                                                               <NA>
#> 1007                              Phenotypic Diversity of Cytomegalovirus DNA Polymerase Gene Variants Observed after Antiviral Therapy
#> 1008 Variations in the cytomegalovirus DNA polymerase and phosphotransferase genes in relation to foscarnet and ganciclovir sensitivity
#> 1018                                                                                                                               <NA>
#> 1019                                                                                                                               <NA>
```

### fasta sequences

``` r
# herpesdrg accepts sequence fragments or whole-genomes, one fasta sequence at a time.

# load example data
my_sequence = system.file("testdata", "HSV1_F716L.fasta", package = "herpesdrg")

dat = call_resistance(infile = my_sequence, all_mutations = F,virus = "HSV1")
#> [1] "fasta found"

head(dat[,c("change", "freq", "REFCODON", "VARCODON", "virus","Aciclovir", "Cidofovir", "Pencyclovir", "ref_doi")])
#>       change freq REFCODON VARCODON virus Aciclovir Cidofovir Pencyclovir
#> 1 UL30_F716L 100%      TTC      CTC  HSV1         6       0.1         2.2
#>                              ref_doi
#> 1 doi: 10.1128/AAC.49.2.606-611.2005
```

### other features

``` r
## view the full database
db = herpesdrg_data()
head(db$aa_change)
#> [1] "D301N" "C304S" "N408D" "N408K" "N408S" "N410K"


## run the shiny application
# ShinyHerpes()
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on the GitHub Issues page. For questions and other
discussions feel free to contact. [Oscar Charles -
maintainer](mailto:oscar.charles.18@ucl.ac.uk)
