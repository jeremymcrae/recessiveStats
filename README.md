### recessiveStats: analysis of recessive variation

#### Estimating cumulative allele frequency rare lof and rare functional variants in each gene.
- ExAC
- DDD unaffected parents
- INTERVAL (?)

#### Estimating probability of HPO matching
- Generate MaxIC scores.
- Estimate p values of phenotype matching observed between N patients drawn from
  M families by simulation.

#### Estimate probability of ‘suspected syndrome’ matching
- Inspect ‘suspected syndrome’ field and prepare it for manually curation to
  generate comma separated lists of correctly spelled clinical syndromes.
- Manually curate dataset, any uncertainties are flagged clinical review.
- Generate pairwise matrix of similarity based on exact matches between
  clinical syndromes and frequency with which each syndrome is suspected.

#### Installation

The package can be installed or updated using R 3.1.0 or greater with:
```R
library(devtools) # if necessary install with install.packages("devtools")
devtools::install_github("jeremymcrae/recessiveStats")

# Alternatively, clone the repository, run R 3.1 from within the top level of
# the repository and use the devtools::build() to build the package for other R
# versions.
```

#### Usage

```R
# load the package
library(recessiveStats)

hgnc ="DNAH14"
chrom = "1"
lof_lof = 5
lof_func = 6
analyse_inherited_enrichment(hgnc, chrom, lof_lof, lof_func)
```
