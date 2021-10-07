library(gdsfmt)
library(SNPRelate)
library(SeqArray)

path.to.kopp.bcftools.full.vcf = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp_production/20210610_bcftools_genotyping/bcftools_vcf/bcftools.call.concat.vcf.gz"
path.to.kopp.bcftools.thinned.vcf = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp_production/20210627_relatedness_output/thinned_vcf/relatedness.isec.0002.exMito.thinned.vcf"
path.to.kopp.bcftools.thinned.no.sex.vcf = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp_production/20210728_no_replicate_no_sex_chrom_vcfs/nextflow_output/thinned_vcf/relatedness.isec.0002.exMito.thinned.vcf"

snpgdsVCF2GDS(path.to.kopp.bcftools.thinned.no.sex.vcf, "kopp.bcftools.thinned.no.sex.snprelate.gds", method="biallelic.only")

kopp.bcftools.thinned.no.sex.gds = snpgdsOpen("kopp.bcftools.thinned.no.sex.snprelate.gds")

snpgdsSummary(kopp.bcftools.thinned.no.sex.gds)

set.seed(1000)

snpset <- snpgdsLDpruning(kopp.bcftools.thinned.no.sex.gds, ld.threshold=0.2, autosome.only=FALSE)

str(snpset)

snpset.id <- unlist(unname(snpset))

# This can be run using the LD snpset.id or without it if you want to keep
# all SNPs in

pca.LD <- snpgdsPCA(kopp.bcftools.thinned.no.sex.gds, snp.id=snpset.id, num.thread=100, autosome.only=FALSE)

pca.no.LD <- snpgdsPCA(kopp.bcftools.thinned.no.sex.gds, num.thread=100, autosome.only=FALSE)


pc.percent <- pca.LD$varprop*100

head(round(pc.percent, 2))

# make a data.frame
tab <- data.frame(sample.id = pca.LD$sample.id,
    EV1 = pca.LD$eigenvect[,1],    # the first eigenvector
    EV2 = pca.LD$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)

head(tab)

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
