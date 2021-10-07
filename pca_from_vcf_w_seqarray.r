library(gdsfmt)
library(SNPRelate)
library(SeqArray)

path.to.kopp.bcftools.full.vcf = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp_production/20210610_bcftools_genotyping/bcftools_vcf/bcftools.call.concat.vcf.gz"
path.to.kopp.bcftools.thinned.vcf = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp_production/20210627_relatedness_output/thinned_vcf/relatedness.isec.0002.exMito.thinned.vcf"
path.to.kopp.bcftools.thinned.no.sex.vcf = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp_production/20210728_no_replicate_no_sex_chrom_vcfs/nextflow_output/thinned_vcf/relatedness.isec.0002.exMito.thinned.vcf"

# vcf.path.to.use = path.to.kopp.bcftools.thinned.vcf

seqVCF2GDS(path.to.kopp.bcftools.thinned.no.sex.vcf, "kopp.bcftools.thinned.no.sex.seqarray.gds")

kopp.bcftools.thinned.no.sex.gds = seqOpen("kopp.bcftools.thinned.no.sex.seqarray.gds")
seqSummary(kopp.bcftools.thinned.no.sex.gds, "genotype")

# this is multiprocessed so set it to what ever number of cores you want (replace 128)
genmat <- seqParallel(128, kopp.bcftools.thinned.no.sex.gds, FUN = function(f)
    {
        s <- 0  # covariance variable with an initial value
        seqBlockApply(f, "$dosage", function(x)
            {
                p <- 0.5 * colMeans(x, na.rm=TRUE)     # allele frequencies (a vector)
                g <- (t(x) - 2*p) / sqrt(p*(1-p))      # normalized by allele frequency
                g[is.na(g)] <- 0                       # correct missing values
                s <<- s + crossprod(g)                 # update the cov matrix s in the parent environment
            }, margin="by.variant")
        s  # output
    }, .combine = "+",    # sum "s" of different processes together
    split = "by.variant")

# scaled by the number of samples over the trace
genmat <- genmat * (nrow(genmat) / sum(diag(genmat)))

# eigen-decomposition
eig <- eigen(genmat, symmetric=TRUE)

# The very basic plot
plot(eig$vectors[,1], eig$vectors[,2], xlab="PC 1", ylab="PC 2")

# The eigen vectors are accessible with eig$vectors
# The eigen values are accessible with eig$values

