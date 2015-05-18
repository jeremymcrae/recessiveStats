
export BCFTOOLS_PLUGINS=/nfs/users/nfs_j/jm33/apps/bcftools/plugins
TEMP_DIR=/lustre/scratch113/projects/ddd/users/jm33/bcfs
VCF_LIST=${TEMP_DIR}/vcf_list.txt
vcfs=`ls ${TEMP_DIR}/*.vcf.gz`
echo ${vcfs} > ${VCF_LIST}
sed -i -- 's/ /\n/g' ${VCF_LIST}

bcftools merge \
    --file-list ${VCF_LIST} \
| bcftools +missing2ref \
    --output-type b \
    --output ${TEMP_DIR}/ddd_4k.bcftools.bcf
