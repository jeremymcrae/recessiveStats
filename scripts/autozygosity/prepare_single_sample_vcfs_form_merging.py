""" script to generate a multi-sample VCF (or BCF) from the thousands of single
sample VCFs for the unaffected probands from the DDD cohort. We use the
annotated uber-VCFs (containing variant calls from multiple callers (but
primarily called by GATK)).

Unfortunately, the existing VCFs cannot be directly merged using either vcftools
(vcf-merge) or bcftools (bcftool merge):
 - vcftools: vcftools merge eventually fails if given sufficent samples to merge.
             This is possibly due to multiallelic variants, but the error given
             is "Could not determine the ploidy (nals=3, nvals=5)".
 - bcftools: bcftools merge cannot natively handle the DDD single-sample VCFs.
             I think this is due to differences between the field formats
             specified in the VCF header, versus the actual INFO and FORMAT
             fields. Even trying to merge two VCFs fails with "Segmentation fault".

My solution is to strip out all of the non genotype data. The means removing all
INFO fields, and all but for the GT field from the FORMAT. I use vcf-annotate to
drop all the redunant fields (since "bcftools annotate" segfaults on the single
sample VCFs.) I exclude all non-PASS variants with "bcftools view". The stripped
down VCFs can then be merged with "bcftools merge".
"""

import os
import subprocess
import time
import random

BCFTOOLS = "/software/hgi/pkglocal/bcftools-1.2/bin/bcftools"
VCFANNOTATE = "/software/hgi/pkglocal/vcftools-0.1.11/bin/vcf-annotate"
TABIX = "/software/hgi/pkglocal/tabix-git-1ae158a/bin/tabix"
TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bcfs/"
EXCLUDE_FIELDS = "FORMAT/DP,FORMAT/AD,INFO/AFR_AF,INFO/Allele,INFO/Amino_acids,INFO/AMR_AF,INFO/AN,INFO/BaseQRankSum,INFO/BIOTYPE,INFO/CALLSOURCE,INFO/CANONICAL,INFO/CCC,INFO/CCDS,INFO/cDNA_position,INFO/CDS_position,INFO/CIEND,INFO/CIPOS,INFO/CLIN_SIG,INFO/ClippingRankSum,INFO/CNSOLIDATE_PASS,INFO/Codons,INFO/COMMON1KG,INFO/COMMONBACKWARDS,INFO/COMMONFORWARDS,INFO/Consequence,INFO/Conserved,INFO/CONVEX,INFO/CONVEX_PASS,INFO/CONVEXSCORE,INFO/CQ,INFO/CQ_TRANSCRIPT_COUNT_CQ,INFO/CQ_TRANSCRIPT_COUNT_GENE,INFO/culprit,INFO/DB,INFO/DDD_AF,INFO/DENOVO-INDEL,INFO/DENOVO-SNP,INFO/DISTANCE,FORMAT/DNM_CONFIG_child_mom_dad,INFO/DP,FORMAT/DP_CHILD,FORMAT/DP_FATHER,FORMAT/DP_MOTHER,FORMAT/DP4_CHILD,INFO/DP4_FATHER,INFO/DP4_MOTHER,INFO/DS,INFO/EAS_AF,INFO/END,INFO/Enhancer,INFO/ENSG,INFO/ENSP,INFO/ENSR,INFO/ENST,INFO/ESP_AF,INFO/EUR_AF,INFO/Existing_variation,INFO/EXON,INFO/FS,FORMAT/gatk_DP,FORMAT/GQ,INFO/GQ_MEAN,INFO/GQ_STDDEV,INFO/HaplotypeScore,INFO/Heart,INFO/HETSNPS,INFO/HGNC,INFO/HGNC_ALL,INFO/HGNC_ID,INFO/HGVSc,INFO/HGVSp,INFO/HIGH_INF_POS,INFO/HWP,INFO/IMPRECISE,INFO/InbreedingCoeff,INFO/inFatherVCF,FORMAT/INHERITANCE,FORMAT/INHERITANCEP,INFO/inMotherVCF,INFO/INTERNALFREQ,INFO/INTRON,INFO/MADL2R,INFO/MAX_AF,INFO/MEANLR2,FORMAT/MIN_DP,INFO/MLEAC,INFO/MLEAF,INFO/MOTIF_NAME,INFO/MOTIF_POS,INFO/MOTIF_SCORE_CHANGE,INFO/MQ,INFO/MQ0,INFO/MQRankSum,INFO/NCC,INFO/NEGATIVE_TRAIN_SITE,INFO/NOVEL,INFO/NUMBERPROBESCONVEX,FORMAT/PL,INFO/PolyPhen,INFO/POSITIVE_TRAIN_SITE,FORMAT/PP_DNM,INFO/Protein_position,INFO/PUBMED,INFO/QD,INFO/RARE,INFO/RAREBACKWARDS,INFO/RAREFORWARDS,INFO/RC50INTERNALFREQ,INFO/ReadPosRankSum,FORMAT/samtools_DP,INFO/SAS_AF,FORMAT/SB,INFO/segmentaldup,INFO/SIFT,INFO/SOMATIC,INFO/STRAND,INFO/SVLEN,INFO/SVTYPE,INFO/SYMBOL,INFO/SYMBOL_SOURCE,INFO/TRF,INFO/UK10K_cohort_AF,INFO/UK10KFREQ,INFO/VQSLOD,INFO/AC"
DIAGNOSED_PATH = "/lustre/scratch113/projects/ddd/users/jm33/ddd_likely_diagnosed.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"

def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that form a unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start. Note that
            the list can be empty, in which case there are no dependencies.
        memory: minimum memory requirements (in megabytes)
        requeue_code: exit codes under which it is acceptable to requeue the
            job, or None
        logfile: path for bjob output
    
    Returns:
        nothing
    """
    
    if job_id is None:
        job_id = get_random_string()
    
    job = "-J \"{0}\"".format(job_id)
    
    mem = ""
    if memory is not None:
        mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
    
    requeue = ""
    if requeue_code is not None:
        requeue = "-Q 'EXCLUDE({0})'".format(requeue_code)
    
    dependent = ""
    if dependent_id is not None:
        if type(dependent_id) == list:
            dependent_id = " && ".join(dependent_id)
        dependent = "-w '{0}'".format(dependent_id)
    
    log = "bjob_output.txt"
    if logfile is not None:
        log = logfile
    
    preamble = ["bsub", job, dependent, requeue, "-q", "normal", "-o", log, mem]
    command = ["bash", "-c", "\""] + command + ["\""]
    
    command = " ".join(preamble + command)
    subprocess.call(command, shell=True)

def is_number(string):
    """ check whether a string can be converted to a number
    
    Args:
        string: value as a string, could be a number
        
    Returns:
        True/False for whether the value can be converted to a number
    """
    
    try:
        number = float(string)
    except ValueError:
        return False
    
    return True

def get_random_string():
    """ make a random string, which we can use for bsub job IDs, so that
    different jobs do not have the same job IDs.
    
    Returns:
        random 8 character string
    """
    
    # set up a random string to associate with the run
    hash_string = "{0:8x}".format(random.getrandbits(32))
    hash_string = hash_string.strip()
    
    # don't allow the random strings to be equivalent to a number, since
    # the LSF cluster interprets those differently from letter-containing
    # strings
    while is_number(hash_string):
        hash_string = "{0:8x}".format(random.getrandbits(32))
        hash_string = hash_string.strip()
    
    return hash_string

def get_vcfs(diagnosed_path, families_path):
    """ get a list of the VCF paths for the unaffected DDD probands.
    """
    
    # identify the set of probands who have diagnoses
    diagnosed = set([])
    with open(diagnosed_path) as handle:
        for line in handle:
            diagnosed.add(line.strip())
    
    vcfs = set([])
    with open(families_path) as handle:
        for line in handle:
            line = line.strip().split("\t")
            # remove parents, diagnosed probands, and the header
            if line[2] == "0" or line[1] in diagnosed or line[2] == "dad_id":
                continue
            
            vcfs.add(line[6])
    
    return sorted(list(vcfs))

def generate_mergeable_vcfs(vcfs, temp_dir):
    """ converts sample VCFs to VCFs that can be merged with bcftools
    
    Args:
        vcfs: list of paths to single-sample vcfs.
        temp_dir: path to temporary directory to to place the new VCFs in.
    
    Returns:
        list of converted VCF paths, and bsub job IDs for the compute jobs that
        are converting them.
    """
    
    stripped_vcfs = []
    job_ids = []
    for vcf in vcfs:
        # identify the sample ID from the VCF (so that we can give the stripped
        # down VCF a sensible filename)
        command = [BCFTOOLS, "query", "--list-samples", vcf]
        sample_id = subprocess.check_output(command).strip()
        
        # get a path to put the VCF in
        stripped_vcf = os.path.join(temp_dir, "{}.vcf.gz".format(sample_id))
        
        # set up the vcftools and bcftools command that will strip out all
        # unecessary INFO and FORMAT fields, then remove all failing variant
        # lines. Finally, tabix index the resulting VCF.
        command = [VCFANNOTATE, \
            "--remove", EXCLUDE_FIELDS, \
            vcf, "|", \
            BCFTOOLS, "view", \
            "--apply-filters", "PASS", \
            "--output-type"," z", \
            "--output-file", stripped_vcf,
            ";", TABIX, "-p", "vcf", "-f", stripped_vcf]
        
        job_id = get_random_string()
        submit_bsub_job(command, job_id)
        time.sleep(0.5)
        stripped_vcfs.append(stripped_vcf)
        job_ids.append(job_id)
    
    return stripped_vcfs, job_ids

def main():
    vcfs = get_vcfs(DIAGNOSED_PATH, FAMILIES_PATH)
    stripped = generate_mergeable_vcfs(vcfs, TEMP_DIR)

if __name___ == "__main__":
    main()
