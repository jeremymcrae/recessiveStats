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
