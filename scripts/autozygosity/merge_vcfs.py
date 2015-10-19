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

I found that it's best to merge the files in stages. Attempting to merge several
thousand VCFs in one step is extremely slow. My strategy iteratively merges pairs
of VCFs. So at the first step, only two samples are merged together, but in the
next round, the VCFs containg two samples, so we get a merged VCF with four
samples in it. We build up to a single VCF containing all the samples. This is
much, much quicker than merging all samples.
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
import argparse

IS_PYTHON3 = sys.version_info[0] == 3

if not IS_PYTHON3:
    itertools.zip_longest = itertools.izip_longest

BCFTOOLS = "/software/hgi/pkglocal/bcftools-1.2/bin/bcftools"
VCFANNOTATE = "/software/hgi/pkglocal/vcftools-0.1.11/bin/vcf-annotate"
TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/bcfs"

# define all the fields to strip from the VCFs (this prevents problems when
# merging)
format_fields = ["AD", "DNM_CONFIG_child_mom_dad", "DP", "DP4_CHILD",
    "DP_CHILD", "DP_FATHER", "DP_MOTHER", "gatk_DP", "GQ", "INHERITANCE",
    "INHERITANCEP", "MIN_DP", "PP_DNM", "samtools_DP", "SB"]
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
DIAGNOSED_PATH = "/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.diagnosed.2015-10-12.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt"


def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to merge all single"
        "sample VCFs for probands in a cohort into a BCF with all samples"
        "contained therein.")
    parser.add_argument("--vcf-annotate", default=VCFANNOTATE, help="path to vcf-annotate binary")
    parser.add_argument("--bcftools", default=BCFTOOLS, help="path to bcftools binary")
    parser.add_argument("--families", default=FAMILIES_PATH, \
        help="Path to listing family relationships, with VCF paths.")
    parser.add_argument("--diagnosed", help="Path to listing diagnosed probands,"
        "if you want to exclude the diagnosed probands.")
    parser.add_argument("--temp", default=TEMP_DIR, \
        help="temporary folder to store intermediate merged VCFs")
    parser.add_argument("--output", help="Path to write output bcf to.")
    
    args = parser.parse_args()
    
    return args

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

def get_vcfs(families_path, diagnosed_path=None):
    """ get a list of the VCF paths for the unaffected DDD probands.
    """
    
    # identify the set of probands who have diagnoses
    diagnosed = set([])
    if diagnosed_path is not None:
        with open(diagnosed_path) as handle:
            for line in handle:
                diagnosed.add(line.strip())
    
    vcfs = set([])
    with open(families_path) as handle:
        for line in handle:
            line = line.strip().split("\t")
            # remove parents, diagnosed probands, and the header
            if line[2] == "0" or line[3] == "0" or line[1] in diagnosed or line[2] == "dad_id":
                continue
            
            vcfs.add(line[6])
    
    return sorted(list(vcfs))

def generate_mergeable_vcfs(vcfs, temp_dir, vcfannotate, bcftools):
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
        command = [bcftools, "query", "--list-samples", vcf]
        try:
            sample_id = subprocess.check_output(command).strip()
        except:
            continue
        
        # get a path to put the VCF in
        stripped_vcf = os.path.join(temp_dir, "{}.bcf".format(sample_id))
        
        # set up the vcftools and bcftools command that will strip out all
        # unecessary INFO and FORMAT fields, then remove all failing variant
        # lines.
        command = [vcfannotate, \
            "--remove", EXCLUDE_FIELDS, \
            vcf, "|", \
            bcftools, "view", \
            "--apply-filters", "PASS", \
            "--max-alleles", "2", \
            "--output-type", "b", \
            "--output-file", stripped_vcf, \
            ";", bcftools, "index", stripped_vcf]
        
        job_id = get_random_string()
        submit_bsub_job(command, job_id, memory=100)
        time.sleep(0.25)
        stripped_vcfs.append(stripped_vcf)
        job_ids.append(job_id)
    
    return stripped_vcfs, job_ids

def merge_vcf_pairs(paths, temp_dir, bcftools):
    """ merge pairs of VCFs using bcftools
    
    Args:
        paths: list of paths for VCFs
        temp_dir: path to directory to store intermediate merged VCFs
    
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
            suffix=".bcf")[1]
        new_paths.append(path)
        
        job_id = get_random_string()
        
        # if we have come to the end of the list, and there is only one VCF in
        # the pair, then we only copy the VCF to the temp file, before indexing
        if len(pair) == 1:
            shutil.copyfile(pair[0], path)
            command = [bcftools, "index", path]
        else:
            command = [bcftools, "merge", "--merge", "none", "--output-type", "b", "--output", \
                path, pair[0],  pair[1], ";", bcftools, "index", path]
        
        submit_bsub_job(command, job_id, memory=100)
        time.sleep(0.5)
        
        job_ids.append(job_id)
        
    return new_paths, job_ids

def merge_vcfs(output_path, vcfs, temp_dir, bcftools):
    """ merge a large number of single-sample VCFs into a single multisample VCF
    
    Args:
        vcfs: list of paths to single-sample VCFs
        temp_dir: path to directory to store intermediate merged VCFs
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
    
    shutil.copyfile(vcfs[0], output_path)

def main(args):
    
    vcfs = get_vcfs(args.families, args.diagnosed)
    stripped, job_ids = generate_mergeable_vcfs(vcfs, args.temp, args.vcf_annotate, args.bcftools)
    
    while has_current_jobs(job_ids):
        time.sleep(30)
    
    merge_vcfs(args.output, stripped, args.temp, args.bcftools)

if __name__ == "__main__":
    try:
        args = get_options()
        main(args)
    finally:
        # make sure we clean up any temporary files
        temp_vcfs = glob.glob(os.path.join(args.temp, "tmp_merge*"))
        for path in temp_vcfs:
            os.remove(path)
