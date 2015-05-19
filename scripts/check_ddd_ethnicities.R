# script to plot PCA data from 1000 genomes and DDD probands, and count the
# number of DDD probands in each 1000 genomes continental group

library(Cairo)
library(mclust)
library(ggplot2)
library(grid)

# define the data paths
KGENOMES_PCA_PATH = "data-raw/pca.1000g.20150511.txt"
DDD_PCA_PATH = "data-raw/pca.ddd_projection.20150511.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"
SANGER_IDS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/person_sanger_decipher.txt"
KGENOMES_INFO_PATH = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt"
DIAGNOSED_PROBANDS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/ddd_likely_diagnosed.txt"

# load the DDD cohort files, so we can identify which samples are the probands,
# and map from their DDD ID to Decipher and sanger IDs
families = read.table(FAMILIES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
families = families[families$dad_id != 0, ]
sanger_ids = read.table(SANGER_IDS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
families = merge(families, sanger_ids, by.x="individual_id", by.y="person_stable_id", all.x=TRUE)
families = families[, c("individual_id", "sanger_id")]
diagnosed = read.table(DIAGNOSED_PROBANDS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
families = families[!families$individual_id %in% diagnosed$person_id, ]

# load the PCA results for the DDD and thousand genomes
kgenomes = read.table(KGENOMES_PCA_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
ddd = read.table(DDD_PCA_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
ddd = ddd[ddd$sample.id %in% families$sanger_id, ]

# load the 1000 genomes sample information, so we can identify which sample
# belongs to which population
kgenomes_info = read.table(KGENOMES_INFO_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
kgenomes_info = kgenomes_info[, c("Sample", "Population")]

# recode all the specific 1000 Genomes populations into continental populations
kgenomes_info$Population[kgenomes_info$Population %in% c("MXL", "PUR", "CLM", "PEL")] = "AMR"
kgenomes_info$Population[kgenomes_info$Population %in% c("GIH", "PJL", "BEB", "STU", "ITU")] = "SAS"
kgenomes_info$Population[kgenomes_info$Population %in% c("CHB", "JPT", "CHS", "CDX", "KHV")] = "EAS"
kgenomes_info$Population[kgenomes_info$Population %in% c("CEU", "TSI", "FIN", "GBR", "IBS")] = "EUR"
kgenomes_info$Population[kgenomes_info$Population %in% c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")] = "AFR"
names(kgenomes_info) = c("sample.id", "population")

# make sure that we have the appropriate contintental lable for each 1000
# Genomes sample. Also, drop the AMR population, since they are sparsely
# distributed around the SAS and EUR populations.
kgenomes = merge(kgenomes, kgenomes_info, by="sample.id", all.x=TRUE)
kgenomes = kgenomes[kgenomes$population != "AMR", ]

# run a clustering algorithm on the 1000 Genomes PCA data
clustering = MclustDA(kgenomes[, c("pca1", "pca2")], kgenomes$population, modelType="MclustDA")

# predict the group that each DDD proband belongs to
ddd_predictions = predict(clustering, ddd[, c("pca1", "pca2")])
ddd$classification = ddd_predictions$classification

# plot the PCA clustering, to assess how well the DDD probands have been predicted
Cairo(file="1000_genomes_pca.pdf", type="pdf", height=15, width=15, units="cm")
p = ggplot(kgenomes, aes(x=pca1, y=pca2, color=population))
p = p + geom_point()
p = p + theme_classic()
p = p + theme(panel.border=element_rect(colour="black", fill=NA, size=1))
p = p + geom_point(data=ddd, aes(x=pca1, y=pca2, color=factor(ddd$classification)), alpha=0.40)
p = p + theme(legend.position = c(0, 1), legend.justification = c(0, 1))
p = p + theme(axis.ticks.length=unit(0.4, "cm"))
print(p)
dev.off()

print(table(ddd$classification))
