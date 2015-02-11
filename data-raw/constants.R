

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
PHENOTYPES_PATH = file.path(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")
EXAC_PATH = "/lustre/scratch113/resources/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz"

devtools::use_data(DATAFREEZE_DIR, PHENOTYPES_PATH, EXAC_PATH, internal=TRUE, overwrite=TRUE)
