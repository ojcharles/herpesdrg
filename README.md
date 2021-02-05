
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
my_sample = system.file("testdata", "HSV1_F716L.vcf", package = "herpesdrg")

data = call_resistance(infile = my_sample, all_mutations = T, virus = "HSV1") # options are HSV1 HSV2 HCMV and VZV

data[ , c("change", "freq", "Aciclovir", "Pencyclovir", "Foscarnet")]
#>        change freq Aciclovir Pencyclovir Foscarnet
#> 1 UL30_A1235A 100%      <NA>        <NA>      <NA>
#> 2    UL30_F2F 100%      <NA>        <NA>      <NA>
#> 3  UL30_F716L 100%         6         2.2        10
#> 4 UL30_G1006G 100%      <NA>        <NA>      <NA>
#> 5  UL31_M290T 100%      <NA>        <NA>      <NA>


## call all variants
mutations_all = call_resistance(infile = my_sample, all_mutations = T)

#to view all mutations in resistance genes we can filter
mutations_res = mutations_all[mutations_all$GENEID %in% c("UL30"),]

head(mutations_res[,c(1,8,21,32:40)])
#>        change freq   CONSEQUENCE Ganciclovir Aciclovir Cidofovir Foscarnet
#> 1 UL30_A1235A 100%    synonymous        <NA>      <NA>      <NA>      <NA>
#> 2    UL30_F2F 100%    synonymous        <NA>      <NA>      <NA>      <NA>
#> 3  UL30_F716L 100% nonsynonymous                     6       0.1        10
#> 4 UL30_G1006G 100%    synonymous        <NA>      <NA>      <NA>      <NA>
#>   Brincidofovir Letermovir Brivudine Pencyclovir Tomeglovir
#> 1            NA       <NA>      <NA>        <NA>         NA
#> 2            NA       <NA>      <NA>        <NA>         NA
#> 3            NA                              2.2         NA
#> 4            NA       <NA>      <NA>        <NA>         NA


## call nonsynonymous variants
# are there any non-synonymous (DNA variants that result in a change of amino acid) variants in resistance genes
mutations_res_nonsyn = mutations_res[mutations_res$CONSEQUENCE == "nonsynonymous",]


# here the top 3 mutations are nonsynonymous, with no identified resistance effect.
head(mutations_res_nonsyn[,c(1,8,21,32:40)])
#>       change freq   CONSEQUENCE Ganciclovir Aciclovir Cidofovir Foscarnet
#> 3 UL30_F716L 100% nonsynonymous                     6       0.1        10
#>   Brincidofovir Letermovir Brivudine Pencyclovir Tomeglovir
#> 3            NA                              2.2         NA
```

### fasta sequences

``` r
# herpesdrg accepts sequence fragments or whole-genomes, one fasta sequence at a time.

# load example data
my_sequence = system.file("testdata", "HSV1_F716L.vcf", package = "herpesdrg")

dat = call_resistance(infile = my_sequence, all_mutations = T)

head(dat)
#>        change    seqnames start   end width strand         id freq RefCount
#> 1 UL30_A1235A NC_001806.2 66511 66511     1      + single run 100%        0
#> 2    UL30_F2F NC_001806.2 62812 62812     1      + single run 100%        0
#> 3  UL30_F716L NC_001806.2 64952 64952     1      + single run 100%        0
#> 4 UL30_G1006G NC_001806.2 65824 65824     1      + single run 100%        0
#> 5  UL31_M290T NC_001806.2 66511 66511     1      - single run 100%        0
#>   VarCount VarAllele varAllele CDSLOC.start CDSLOC.end CDSLOC.width PROTEINLOC
#> 1        1         G         G         3705       3705            1       1235
#> 2        1         C         C            6          6            1          2
#> 3        1         C         C         2146       2146            1        716
#> 4        1         C         C         3018       3018            1       1006
#> 5        1         G         C          869        869            1        290
#>   QUERYID TXID  CDSID GENEID   CONSEQUENCE REFCODON VARCODON REFAA VARAA
#> 1       4   15 18, 62   UL30    synonymous      GCA      GCG     A     A
#> 2       1   15     18   UL30    synonymous      TTT      TTC     F     F
#> 3       2   15     18   UL30 nonsynonymous      TTC      CTC     F     L
#> 4       3   15     18   UL30    synonymous      GGA      GGC     G     G
#> 5       4   62 18, 62   UL31 nonsynonymous      ATG      ACG     M     T
#>   aachange mutation_id virus genotype gene aa_change Ganciclovir Aciclovir
#> 1   A1235A          NA  <NA>     <NA> <NA>      <NA>        <NA>      <NA>
#> 2      F2F          NA  <NA>     <NA> <NA>      <NA>        <NA>      <NA>
#> 3    F716L        1508  HSV1      TAS UL30     F716L                     6
#> 4   G1006G          NA  <NA>     <NA> <NA>      <NA>        <NA>      <NA>
#> 5    M290T          NA  <NA>     <NA> <NA>      <NA>        <NA>      <NA>
#>   Cidofovir Foscarnet Brincidofovir Letermovir Brivudine Pencyclovir Tomeglovir
#> 1      <NA>      <NA>            NA       <NA>      <NA>        <NA>         NA
#> 2      <NA>      <NA>            NA       <NA>      <NA>        <NA>         NA
#> 3       0.1        10            NA                              2.2         NA
#> 4      <NA>      <NA>            NA       <NA>      <NA>        <NA>         NA
#> 5      <NA>      <NA>            NA       <NA>      <NA>        <NA>         NA
#>   Maribavir Cyclopropavir
#> 1        NA            NA
#> 2        NA            NA
#> 3        NA            NA
#> 4        NA            NA
#> 5        NA            NA
#>                                                                                                                                                                                            ref_title
#> 1                                                                                                                                                                                               <NA>
#> 2                                                                                                                                                                                               <NA>
#> 3 Genotypic Characterization of the DNA Polymerase and Sensitivity to Antiviral Compounds of Foscarnet-Resistant Herpes Simplex Virus Type 1 (HSV-1) Derived from a Foscarnet-Sensitive HSV-1 Strain
#> 4                                                                                                                                                                                               <NA>
#> 5                                                                                                                                                                                               <NA>
#>                                               ref_link
#> 1                                                 <NA>
#> 2                                                 <NA>
#> 3 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC547286/
#> 4                                                 <NA>
#> 5                                                 <NA>
#>                              ref_doi   strain_generation test_method co_gene
#> 1                               <NA>                <NA>        <NA>    <NA>
#> 2                               <NA>                <NA>        <NA>    <NA>
#> 3 doi: 10.1128/AAC.49.2.606-611.2005 Lab induced strain          PRA        
#> 4                               <NA>                <NA>        <NA>    <NA>
#> 5                               <NA>                <NA>        <NA>    <NA>
#>   co_aa created_date created_by note status    X
#> 1  <NA>         <NA>       <NA> <NA>   <NA> <NA>
#> 2  <NA>         <NA>       <NA> <NA>   <NA> <NA>
#> 3         04/02/2021  OJCharles           A     
#> 4  <NA>         <NA>       <NA> <NA>   <NA> <NA>
#> 5  <NA>         <NA>       <NA> <NA>   <NA> <NA>
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
