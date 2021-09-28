#!/usr/bin/env nextflow


/*
This workflow takes two vcf files as input (GATK- and bcftools-based),
filters down to mito genome only
filters to individual samples
computes the intersect using rtg-tools vcfeval or bcftools isec
then creates a consensus for each of the samples.
*/

gatk_vcfgz = file(params.gatk_vcfgz)
gatk_vcfgz_tbi = file("${gatk_vcfgz}.tbi")
bcftools_vcfgz = file(params.bcftools_vcfgz)
bcftools_vcfgz_tbi = file("${bcftools_vcfgz}.tbi")
mitoscaff = "NC_034374.1"
output_dir = "${workflow.launchDir}/outputs/"
ref = file("/home/humebc/projects/20210125_kopp_guinea_fowl/nummel_ref_assembly/mito.only.GCF_002078875.1_NumMel1.0_genomic.fna")
ref_dict = file("/home/humebc/projects/20210125_kopp_guinea_fowl/nummel_ref_assembly/mito.only.GCF_002078875.1_NumMel1.0_genomic.dict")
ref_fai = file("/home/humebc/projects/20210125_kopp_guinea_fowl/nummel_ref_assembly/mito.only.GCF_002078875.1_NumMel1.0_genomic.fna.fai")
ref_sdf = file("/home/humebc/projects/20210125_kopp_guinea_fowl/nummel_ref_assembly/mito.only.GCF_002078875.1_NumMel1.0_genomic.sdf")
vcf_in_ch = Channel.from(["gatk", gatk_vcfgz, gatk_vcfgz_tbi], ["bcftools", bcftools_vcfgz, bcftools_vcfgz_tbi])
sample_list_file = "/home/humebc/projects/20210125_kopp_guinea_fowl/mito_genome_production/sample_list.txt"
process_program = "rtg"

sample_list = {  
    def sample_list = []
    new File(sample_list_file).eachLine {
        line -> 
        line.split().each{
            sample_list << it;
        }
    }
    return sample_list
}()

sample_list_ch = Channel.from(sample_list)

// Filter down each of the vcf files to just the mitoscaffhold
process filter_to_mito{
    container  "biocontainers/vcftools:v0.1.16-1-deb_cv1"
    if (workflow.containerEngine == 'docker'){
        containerOptions '-u $(id -u):$(id -g)'
    }
    publishDir "${output_dir}mito_only_vcfs/${name}", mode: "copy"

    input:
    tuple val(name), path(vcfgz), path(vcfgztbi) from vcf_in_ch

    output:
    tuple val(name), path("${name}.mito.vcf.gz"), path("${name}.mito.vcf.gz.tbi") into sample_split_no_isec_ch

    script:
    """
    vcftools --chr $mitoscaff --recode --recode-INFO-all --stdout --gzvcf $vcfgz > ${name}.mito.vcf
    bgzip -c ${name}.mito.vcf > ${name}.mito.vcf.gz
    tabix -p vcf ${name}.mito.vcf.gz
    """
}

// // Now split the files by sample
// Do this with GATK as it gives us the ability to remove the alternates that are now unused and exclude the non-variants
// Then isec should hopefully work. But there is also a GATK verison of this that I am now tempted to use
//gatk SelectVariants -R /home/humebc/projects/20210125_kopp_guinea_fowl/nummel_ref_assembly/mito.only.GCF_002078875.1_NumMel1.0_genomic.fna -V /home/humebc/projects/20210125_kopp_guinea_fowl/mito_genome_production/outputs/mito_only_vcfs/gatk/gatk.mito.vcf.gz -sn RGID1_S10 -O gatk.out.vcf.gz --exclude-non-variants --remove-unused-alternates
// Acutally people suggest using rtg tool for the isec

// We need to have them down to the per sample before we do the intersect
// NB for the GATK SelectVariants to work properly there needs to be
// the FORMAT GQ header in the files, despite the fact that it is not used
// the header is : ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
// It seems that we can ask GATK to correct the header using the FixVcfHeader command
process split_to_sample{
    container "broadinstitute/gatk:4.2.2.0"
    publishDir "${output_dir}single_sample_vcfs/${name}", mode: "copy"

    input:
    each sample from sample_list_ch
    file ref
    file ref_dict
    file ref_fai
    tuple val(name), path(vcfgz), path(vcfgztbi) from sample_split_no_isec_ch

    output:
    tuple val(name), val(sample), path("${name}.${sample}.mito.vcf.gz"), path("${name}.${sample}.mito.vcf.gz.tbi") into isec_ch, consensus_no_isec_ch

    script:
    """
    gatk FixVcfHeader -I $vcfgz -O fixed.vcf.gz
    gatk SelectVariants -R $ref -sn $sample -V fixed.vcf.gz -O ${name}.${sample}.mito.vcf.gz --exclude-non-variants --remove-unused-alternates
    """

}


// Create the intersection vcf
// Here we need to end up with sample pairs where we have the gatk and the bcftools version
// of each sample so that we can then get an isec version.
if (process_program == "isec"){
    process isec_bcftools{
        container "dceoy/bcftools:latest"
        publishDir "${output_dir}single_sample_vcfs_bcftools/isec", mode: "copy", overwrite: true

        input:
        tuple val(name), val(sample), path(vcfgz), path(vcfgztbi) from isec_ch.groupTuple(by:1)

        output:
        tuple val("isec"), val(sample), path("isec.${sample}.mito.vcf.gz"), path("isec.${sample}.mito.vcf.gz.tbi") into consensus_isec_ch

        script:
        """
        mkdir isec
        bcftools isec gatk.${sample}.mito.vcf.gz bcftools.${sample}.mito.vcf.gz -p isec
        mv isec/0002.vcf isec.${sample}.mito.vcf
        bgzip -c isec.${sample}.mito.vcf > isec.${sample}.mito.vcf.gz
        tabix -p vcf isec.${sample}.mito.vcf.gz
        """
    }
}else{
    process isec_rtg{
        container "realtimegenomics/rtg-tools:3.12.1"
        publishDir "${output_dir}single_sample_vcfs_rtg/isec", mode: "copy", overwrite: true

        input:
        file ref
        file ref_sdf
        tuple val(name), val(sample), path(vcfgz), path(vcfgztbi) from isec_ch.groupTuple(by:1)

        output:
        tuple val("isec"), val(sample), path("isec.${sample}.mito.vcf.gz"), path("isec.${sample}.mito.vcf.gz.tbi") into consensus_isec_ch

        script:
        """
        rtg vcfeval -b gatk.${sample}.mito.vcf.gz -c bcftools.${sample}.mito.vcf.gz -t $ref_sdf -o rtg_out --squash-ploidy
        mv rtg_out/tp-baseline.vcf.gz isec.${sample}.mito.vcf.gz
        mv rtg_out/tp-baseline.vcf.gz.tbi isec.${sample}.mito.vcf.gz.tbi
        """
    }
}


// Create a consensus for each of the files
process consensus{
    container "dceoy/bcftools:latest"
    publishDir "${output_dir}consensus_sequences/${name}", mode: "copy", overwrite: true
    
    input:
    file ref
    tuple val(name), val(sample), path(vcfgz), path(vcfgztbi) from consensus_isec_ch.mix(consensus_no_isec_ch)

    output:
    path "${name}.${sample}.mito.consensus.fa" into consensus_out_ch

    shell:
    '''
    bcftools consensus --sample !{sample} -f !{ref} !{vcfgz} > out.fa
    awk '/>/ {print ">!{sample}_!{name}_mito";} !/>/ {print $0}' out.fa > !{name}.!{sample}.mito.consensus.fa
    '''
}