

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
PHENOTYPES_PATH = file.path(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")

devtools::use_data(DATAFREEZE_DIR, PHENOTYPES_PATH, internal=TRUE, overwrite=TRUE)
