
export BCFTOOLS_PLUGINS=/nfs/users/nfs_j/jm33/apps/bcftools/plugins
TEMP_DIR=/lustre/scratch113/projects/ddd/users/jm33/bcfs
VCF_LIST=${TEMP_DIR}/vcf_list.txt
vcfs=`ls ${TEMP_DIR}/*.vcf.gz`

# alternatively, use vcftools to merge the VCFs
vcf-merge \
    --ref-for-missing 0/0 \
    ${vcfs} \
| bcftools view \
    --output-type b \
    --output-file ${TEMP_DIR}/ddd_4k.vcftools.bcf
