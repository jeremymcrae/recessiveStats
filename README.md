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

hgnc = "DNAH14"
chrom = "1"
counts = list()
counts$biallelic_lof = 5
counts$biallelic_func = 10
counts$lof_func = 6

# define the variants in the control population, including allele counts
variants = read.table(header = TRUE, text = "
    chrom  pos   AC  AN    CQ
    1      1000  1   1000  missense_variant
    1      1000  5   1000  stop_gained
    1      1000  8   1000  stop_lost
    1      1000  20  1000  synonymous_variant")

analyse_inherited_enrichment(counts, variants, cohort_n)

# For ease of use, you can also load variant counts from ExAC with:
variants = get_exac_variants_for_gene(hgnc, chrom)

# The ExAC variants object should be a list of tables, one for each of the ExAC
# populations, e.g. 'AFR', 'EAS', 'NFE'. You'll need to pick the population that
# matches your cohort. For this example we'll use the 'NFE'.
analyse_inherited_enrichment(counts, variants[['NFE']], cohort_n)

# alternatively, provide your path to the ExAC VCF e.g.
variants = get_exac_variants_for_gene(hgnc, chrom, fileName='PATH_TO_VCF')
```

Rather than naming a gene, you can give a chromosome range (but define the gene
symbol as NULL, otherwise all of the variants that don't match the gene symbol
are removed):
```R
start=225083964
end=225586996
variants = get_exac_variants_for_gene(hgnc=NULL, chrom='1', start=start, end=end)
analyse_inherited_enrichment(counts, variants[['NFE']], cohort_n)

# define your own MAF threshold for variant inclusion
analyse_inherited_enrichment(counts, variants[['NFE']], threshold=0.005, cohort_n)
```

You can also take the autozygosity into account. Calculate the proportion of
probands who have an autozygous region overlapping the gene `bcftools roh` is
recommended). Then you can include the rate as:
```R
RATE=0.005
analyse_inherited_enrichment(counts, variants[['NFE']], cohort_n, autozygosity=RATE)
```

Also, if your probands are of multiple ethnicities, you can account for
differences in allele frequencies between ethnicities by specifying the number
of probands that would be classified as belonging to each ExAC population. For
example:
```R
cohort_n = list("AFR"=100, "EAS"=50, "NFE"=800, "SAS"=50)
analyse_inherited_enrichment(counts, variants, cohort_n)
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
