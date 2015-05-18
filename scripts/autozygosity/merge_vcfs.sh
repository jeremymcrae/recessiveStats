
TEMP_DIR=/lustre/scratch113/projects/ddd/users/jm33/bcfs
VCF_LIST=${TEMP_DIR}/vcf_list.txt
vcfs=`ls ${TEMP_DIR}/*.vcf.gz`
echo ${vcfs} > ${VCF_LIST}
sed -i -- 's/ /\n/g' ${VCF_LIST}

bsub -o merging.bjob_output.txt \
    -R "select[mem>13000] rusage[mem=13000]" -M 13000 \
    -q basement \
    bcftools merge \
        --file-list ${VCF_LIST} \
        --output-type b \
        --output ${TEMP_DIR}/ddd_4k.bcf


bsub -o merging_vcftools.bjob_output.txt \
    -R "select[mem>15000] rusage[mem=15000]" -M 15000 \
    vcf-merge \
        --ref-for-missing 0/0 \
        ${vcfs}
