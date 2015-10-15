

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
DDD_VCFS_DIR = "/lustre/scratch114/projects/ddd/release/20140912/final"
PHENOTYPES_PATH = file.path(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")
BCFTOOLS = "/software/hgi/pkglocal/bcftools-1.0/bin/bcftools"

devtools::use_data(DATAFREEZE_DIR, PHENOTYPES_PATH, DDD_VCFS_DIR, BCFTOOLS,
    internal=TRUE, overwrite=TRUE)
