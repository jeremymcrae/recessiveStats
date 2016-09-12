#' filter ExAC variants by CADD score
#'
#' @param vars data frame or list of data frames with ExAC variants, output from get_exac_variants_for_gene
#' @param path path to vcf.gz file of CADD scores, default = /lustre/scratch115/projects/ddd/users/hcm/DDD/CADD_resources/ExAC_r0.3.tsv.gz
#' @param threshold CADD score cutoff for variants, default = 21
#' @param exclude boolean for whether to exclude variants without a CADD score, default=FALSE
#' @param tempdir path to temporary diectory for writing VCF files; default = /lustre/scratch115/projects/ddd/users/hcm/DDD/tempdir_for_gene_vcfs_for_annotation/
#' @param notLoFs boolean for whether to keep all LoFs regardless of CADD score, default = TRUE

#' @export
#'
#' @return data frame or list of data frame  with ExAC variants, filtered by CADD

filter_exac_variants_by_cadd<- function(vars,path="/lustre/scratch115/projects/ddd/users/hcm/DDD/CADD_resources/ExAC_r0.3.tsv.gz",threshold=21,exclude=FALSE,tempdir="/lustre/scratch115/projects/ddd/users/hcm/DDD/tempdir_for_gene_vcfs_for_annotation/",notLoFs = TRUE) {
  convert.to.list=FALSE
  #if there's only one population in vars, convert to a list
  if(class(vars)!="list"){
    vars=list(vars)
    convert.to.list=TRUE
  }
### there's the same data frame for each population, so we only need to annotate this once
  v=vars[[1]]
  vcfname=paste(tempdir,"variants_in_chr",v[1,1],"_",min(v$POS),"_",max(v$POS),".vcf",sep="")
  outvcfname=paste(tempdir,"variants_in_chr",v[1,1],"_",min(v$POS),"_",max(v$POS),".annotated.vcf",sep="")
#write out a vcf containing the variants
  write.table(v[,c("CHROM","POS","REF","ALT")],vcfname,quote=F,sep="\t",row.names=F)
#add CADD scores to this
  system(paste("python ~/autozyg_code/TabixScores.py --tabix ",path," --variants ",vcfname," --variants_out ",outvcfname," --score CADD",sep=""))
#read in the annotated VCF
  annot.v=read.delim(outvcfname,header=T,stringsAsFactors=F)
#remove the temporary vcf files
  system(paste("rm ",vcfname, " ", outvcfname,sep=""))
  v$name=paste(v[,1],v[,2],v[,3],v[,4],sep="_")
  annot.v$name=paste(annot.v[,1],annot.v[,2],annot.v[,3],annot.v[,4],sep="_")
#add annotations to the original data frame
  v=merge(v,annot.v,by.x="name",by.y="name",all.x=T)
#retain only variants over the CADD threshold
  if(!notLoFs){
    if(exclude){
      keep.v=v[(v$scaled_CADD > threshold & !is.na(v$scaled_CADD)),]
    }else {
      keep.v=v[(v$scaled_CADD > threshold & !is.na(v$scaled_CADD))|is.na(v$scaled_CADD),]
    }
  } else {
#if you want to keep all LOFs regardless of CADD
    v.lofs=v[v$CQ %in% c("transcript_ablation","splice_donor_variant","splice_acceptor_variant","stop_gained","frameshift_variant","splice_region_variant"),]
    v.not.lofs = v[!v$CQ %in% c("transcript_ablation","splice_donor_variant","splice_acceptor_variant","stop_gained","frameshift_variant","splice_region_variant"),]
    if(exclude){
      keep.v=v.not.lofs[(v.not.lofs$scaled_CADD > threshold & !is.na(v.not.lofs$scaled_CADD)),]
    }else {
      keep.v=v.not.lofs[(v.not.lofs$scaled_CADD > threshold & !is.na(v.not.lofs$scaled_CADD))|is.na(v.not.lofs$scaled_CADD),]
    }
    keep.v=data.frame(rbind(v.lofs,keep.v))
  }
  output.vars=list()
  for(i in 1:length(vars)){
    v=vars[[i]]
    v$name=paste(v[,1],v[,2],v[,3],v[,4],sep="_")
    output.vars[[i]] = v[v$name %in% keep.v$name,]
  }
  names(output.vars) = names(vars)
  #if there was just a data frame to start, convert back to data frame
  if(convert.to.list){
    output.vars=output.vars[[1]]
  }
  #return filtered variants
  return(output.vars)
}



