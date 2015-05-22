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
drop all the redundant fields (since "bcftools annotate" segfaults on the single
sample VCFs.) I exclude all non-PASS variants with "bcftools view". The stripped
down VCFs can then be merged with "bcftools merge".
"""

import os
import sys
import glob
import subprocess
import tempfile
import shutil
import itertools
import random
import time

IS_PYTHON3 = sys.version_info[0] == 3

if not IS_PYTHON3:
    itertools.zip_longest = itertools.izip_longest

BCFTOOLS = "/software/hgi/pkglocal/bcftools-1.2/bin/bcftools"
VCFANNOTATE = "/software/hgi/pkglocal/vcftools-0.1.11/bin/vcf-annotate"
TABIX = "/software/hgi/pkglocal/tabix-git-1ae158a/bin/tabix"
TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bcfs"

# define all the fields to strip from the VCFs (this prevents problems when
# merging)
format_fields = ["AD", "DNM_CONFIG_child_mom_dad", "DP", "DP4_CHILD",
    "DP_CHILD", "DP_FATHER", "DP_MOTHER", "gatk_DP", "GQ", "INHERITANCE",
    "INHERITANCEP", "MIN_DP", "PL", "PP_DNM", "samtools_DP", "SB"]
format_fields = ["FORMAT/{}".format(x) for x in format_fields]
info_fields = ["AC", "AFR_AF", "Allele", "Amino_acids", "AMR_AF", "AN",
    "BaseQRankSum", "BIOTYPE", "CALLSOURCE", "CANONICAL", "CCC", "CCDS",
    "cDNA_position", "CDS_position", "CIEND", "CIPOS", "CLIN_SIG",
    "ClippingRankSum", "CNSOLIDATE_PASS", "Codons", "COMMON1KG",
    "COMMONBACKWARDS", "COMMONFORWARDS", "Consequence", "Conserved",
    "CONVEX", "CONVEX_PASS", "CONVEXSCORE", "CQ", "CQ_TRANSCRIPT_COUNT_CQ",
    "CQ_TRANSCRIPT_COUNT_GENE", "culprit", "DB", "DDD_AF", "DENOVO-INDEL",
    "DENOVO-SNP", "DISTANCE", "DP", "DP4_FATHER", "DP4_MOTHER", "DS", "EAS_AF",
    "END", "Enhancer", "ENSG", "ENSP", "ENSR", "ENST", "ESP_AF", "EUR_AF",
    "Existing_variation", "EXON", "FS", "GQ_MEAN", "GQ_STDDEV",
    "HaplotypeScore", "Heart", "HETSNPS", "HGNC", "HGNC_ALL", "HGNC_ID",
    "HGVSc", "HGVSp", "HIGH_INF_POS", "HWP", "IMPRECISE", "InbreedingCoeff",
    "inFatherVCF", "inMotherVCF", "INTERNALFREQ", "INTRON", "MADL2R", "MAX_AF",
    "MEANLR2", "MLEAC", "MLEAF", "MOTIF_NAME", "MOTIF_POS",
    "MOTIF_SCORE_CHANGE", "MQ", "MQ0", "MQRankSum", "NCC",
    "NEGATIVE_TRAIN_SITE", "NOVEL", "NUMBERPROBESCONVEX", "PolyPhen",
    "POSITIVE_TRAIN_SITE", "Protein_position", "PUBMED", "QD", "RARE",
    "RAREBACKWARDS", "RAREFORWARDS", "RC50INTERNALFREQ", "ReadPosRankSum",
    "SAS_AF", "segmentaldup", "SIFT", "SOMATIC", "STRAND", "SVLEN", "SVTYPE",
    "SYMBOL", "SYMBOL_SOURCE", "TRF", "UK10K_cohort_AF", "UK10KFREQ", "VQSLOD"]
info_fields = ["INFO/{}".format(x) for x in info_fields]
EXCLUDE_FIELDS = ",".join(format_fields + info_fields)
DIAGNOSED_PATH = "/lustre/scratch113/projects/ddd/users/jm33/ddd_likely_diagnosed.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"

def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None, rerunnable=False):
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
    
    rerun = ""
    if rerunnable:
        rerun = "-r"
    
    log = "bjob_output.txt"
    if logfile is not None:
        log = logfile
    
    preamble = ["bsub", job, dependent, requeue, rerun, "-q", "normal", "-o", log, mem]
    command = ["bash", "-c", "\""] + command + ["\""]
    
    command = " ".join(preamble + command)
    subprocess.check_call(command, shell=True)

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

def has_current_jobs(job_ids):
    """ find whether any of a list of job IDs are current jobs on the cluster
    
    Args:
        job_ids: list of bsub job names
    
    Returns:
        True/False for whether any of the listed job IDs are current jobs on the
        cluster.
    """
    
    command = ["bjobs", "-o", "\"JOBID", "USER", "STAT", "QUEUE", "JOB_NAME", \
        "delimiter=';'\""]
    command = " ".join(command)
    jobs = subprocess.check_output(command, shell=True, stderr=open(os.devnull))
    
    # if there aren't any currently running or pending jobs, then the output
    if jobs == "":
        return False
    
    jobs = jobs.decode()
    jobs = jobs.strip().split("\n")
    
    current_jobs = set([])
    for line in jobs:
        # ignore the header line
        if line.startswith("JOBID"):
            continue
        
        line = line.split(";")
        status = line[2]
        queue = line[3]
        job_name = line[4]
        current_jobs.add(job_name)
    
    return len(set(job_ids) & current_jobs) > 0

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

def merge_vcf_pairs(paths, temp_dir):
    """ merge pairs of VCFs using bcftools
    
    Args:
        paths: list of paths for VCFs
    
    Returns:
        list of paths to merged VCFs, and list of bsub job IDs, so that we can
        track job completetion before running any additional rounds of VCF
        merging.
    """
    
    new_paths = []
    job_ids = []
    
    # get a list of the VCF paths in groups of two
    args = [iter(paths)] * 2
    pairs = itertools.zip_longest(*args)
    pairs = [[path for path in pair if path != None] for pair in pairs]
    
    for pair in pairs:
        path = tempfile.mkstemp(dir=temp_dir, prefix="tmp_merge.",
            suffix=".vcf.gz")[1]
        new_paths.append(path)
        
        job_id = get_random_string()
        
        # if we have come to the end of the list, and there is only one VCF in
        # the pair, then we only copy the VCF to the temp file, before indexing
        if len(pair) == 1:
            shutil.copyfile(pair[0], path)
            command = [BCFTOOLS, "index", "--tbi", path]
        else:
            command = [BCFTOOLS, "merge", "--output-type", "z", "--output", \
                path, pair[0],  pair[1], ";", BCFTOOLS, "index", "--tbi", path]
        
        submit_bsub_job(command, job_id)
        time.sleep(0.5)
        
        job_ids.append(job_id)
        
    return new_paths, job_ids

def merge_vcfs(vcfs, temp_dir):
    """ merge a large number of single-sample VCFs into a single multisample VCF
    
    Args:
        vcfs: list of paths to single-sample VCFs
    """
    
    # merge the initial list of VCFs, which gives a smaller list of VCFs,
    # containing >1 sample per VCF. These can be merged again, and repeated,
    # until we end up with a single VCF, containing all the initial samples.
    while len(vcfs) > 1:
        vcfs, job_ids = merge_vcf_pairs(vcfs, temp_dir)
        
        # add a delay, so that the jobs can be submitted to the cluster
        time.sleep(10)
        while has_current_jobs(job_ids):
            time.sleep(30)
    
    # finally, swap the missing genotypes to reference alleles, and convert to BCF
    command = [BCFTOOLS, "+missing2ref", \
        "--output-type", "b",
        "--output", \
            os.path.join(os.path.dirname(temp_dir)[0], "ddd_4k.bcftools.bcf"), \
        vcfs[0]]
    submit_bsub_job(command)

def main():
    vcfs = get_vcfs(DIAGNOSED_PATH, FAMILIES_PATH)
    stripped, job_ids = generate_mergeable_vcfs(vcfs, TEMP_DIR)
    
    while has_current_jobs(job_ids):
        time.sleep(30)
    
    merge_vcfs(stripped, TEMP_DIR)

if __name__ == "__main__":
    main()
