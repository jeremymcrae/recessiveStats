[![Build Status](https://travis-ci.org/jeremymcrae/recessiveStats.svg?branch=master)]
(https://travis-ci.org/jeremymcrae/recessiveStats)

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

The baseline rate of biallelic LoF events in a gene is the rate of rare LoF
variation squared, whereas the baseline rate of LoF/Func is the rate of Lof
variation multiplied by the sum of the rate of LoF variation and functional
variation (LoF * (LoF + Func)).

The probability of getting more than or equal to the number of the observed
inherited events in each class was estimated using a binomial model where
n=number of families with undiagnosed probands and k=baseline rate of variation.

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
cohort_n = 1000

hgnc ="DNAH14"
chrom = "1"
biallelic_lof = 5
biallelic_func = 10
lof_func = 6
probands = c("DDDP00000X", "DDDP00000Y", "DDDP00000Z")
analyse_inherited_enrichment(hgnc, chrom, biallelic_lof, biallelic_func, lof_func, probands=probands, cohort_n=cohort_n)
```

Rather than naming a gene, you can give a chromosome range (but define the gene
symbol as NULL, otherwise all of the variants that don't match the gene symbol
are removed):
```R
start=225083964
end=225586996
analyse_inherited_enrichment(hgnc=NULL, chrom, biallelic_lof, biallelic_func, lof_func, start=start, end=end, probands=probands, cohort_n=cohort_n)
```

You can also take the autozygosity into account. Calculate the proportion of
probands who have an autozygous region overlapping the gene `bcftools roh` is
recommended). Then you can include the rate as:
```R
RATE=0.005
analyse_inherited_enrichment(hgnc, chrom, biallelic_lof, biallelic_func, lof_func, probands=probands, autozygous_rate=RATE)
```

Also, if your probands are of multiple ethnicities, you can account for
differences in allele frequencies between ethnicities by specifying the number
of probands that would be classified as belonging to each ExAC population. For
example:
```R
cohort_n = list("AFR"=100, "EAS"=50, "NFE"=800, "SAS"=50)
analyse_inherited_enrichment(hgnc, chrom, biallelic_lof, biallelic_func, lof_func, probands=probands, cohort_n=cohort_n)
```

The populations available in ExAC are:

 code | description
----- | --------------------
 AFR  | African/African American
 AMR  | American
 EAS  | East Asian
 FIN  | Finnish
 NFE  | Non-Finnish European
 OTH  | Other
 SAS  | South Asian
