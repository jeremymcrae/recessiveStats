### recessiveStats: analysis of recessive variation

This package analyses genes suspected of having a recessive mode of inheritance.
We test for enrichment of compound heterozygous variants in genes and similarity
of suspected syndromes between probands. A combined test is calculated by
incorporating results from testing for similarity of Human Phenotype Ontology
(HPO) terms among probands (see https://github.com/jeremymcrae/hpo_similarity).

##### Enrichment of recessive variants in genes
We identified probands with rare compound heterozygous variants in genes
(where at least one of the compound het pair has a loss-of-function consequence),
or rare homozygous loss-of-function variants. These variants were identified
with https://github.com/jeremymcrae/clinical-filter. Probands who had diagnoses
were excluded from further analysis, as any additional recessive variants cannot
contribute to their disorder. The numbers of probands in each of the functional
categories were tallied for each gene. Only one proband per family per gene was
permitted, so that the probands contributing to the tallies were independent
evidence.

The numbers of probands were tallied in genes for two functional categories:
- **LoF/Lof:** homozygous loss-of-function (LoF) and compound hets where both
    variants were LoF.
- **LoF/Func:** compound hets where only one variant is LoF.

The baseline rate of rare LoF/LoF and rare LoF/Func variation was determined
from two sources:
- [ExAC](http://exac.broadinstitute.org/) version 0.3 datasets.
- DDD unaffected parents (for each gene, the parents of the probands with
    variants in that gene were excluded).

The frequency of rare (allele frequency < 0.01) loss-of-function and rare
functional variants in each gene was calculated (gene coordinates were obtained
from gencode v19). The frequencies of the variants in each gene for these
classes were summed, to get an overall rate of rare functional and rare LoF
variation for each gene.

The baseline rate of LoF/LoF events in a gene is the rate of rare LoF variation
squared, whereas the baseline rate of LoF/Func is the rate of Lof variation
multiplied by the sum of the rate of LoF variation and functional variation
(LoF/LoF * (LoF/LoF + LoF/Func)).

The probability of getting more than or equal to the number of the observed
inherited events in each class was estimated using a binomial model where
n=number of families with undiagnosed probands and k=baseline rate of variation.

##### Suspected syndrome matches
- Inspect ‘suspected syndrome’ field and prepare it for manually curation to
  generate comma separated lists of correctly spelled clinical syndromes.
- Manually curate dataset, any uncertainties are flagged clinical review.
- Generate pairwise matrix of similarity based on exact matches between
  clinical syndromes and frequency with which each syndrome is suspected.

##### Installation
The package can be installed or updated using R 3.1.0 or greater with:
```R
library(devtools) # if necessary install with install.packages("devtools")
devtools::install_github("jeremymcrae/recessiveStats")

# Alternatively, clone the repository, run R 3.1 from within the top level of
# the repository and use the devtools::build() to build the package for other R
# versions.
```

##### Usage
```R
# load the package
library(recessiveStats)

hgnc ="DNAH14"
chrom = "1"
lof_lof = 5
lof_func = 6
probands = c("DDDP00000X", "DDDP00000Y", "DDDP00000Z")
analyse_inherited_enrichment(hgnc, chrom, lof_lof, lof_func, probands)
```
