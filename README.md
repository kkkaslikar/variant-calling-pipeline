# Introduction

This is a repository for a pipeline written for a variant-calling on a down-sampled data set from [this tutorial](https://melbournebioinformatics.github.io/MelBioInf_docs/tutorials/var_detect_advanced/var_detect_advanced/). The downsampled data set consisted of reads mapping only to chromosome 20.

The pipeline calls variants using two tools, the final outputs of which are later compared:

1. `bcftools`
2. `FreeBayes`

The output from both tools is saved to separate files. At this point, the pipeline bifurcates into two parallel tracks, one track for the output from each tool.

The following instructions are applicable to both tracks:

* After variant-calling, the variants (from both tools) are filtered based on variant quality and depth at the variant site. In addition, they are filtered to include only SNPs, leaving out other kinds of variants such as indels.
* Then, the filtered  variants are annotated with `snpEff`.
* The annotated variants are converted into bgzipped format and indexed using `bcftools`.

Then, the final annotated outputs from both tracks intersected with `bcftools isec` to get the shared annotated snp sites.

# Starting point

The starting point for this pipeline are BAM files which have already been aligned to the reference genome `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, obtained from Ensembl. Mapping has already been performed by bowtie2.

# Files not uploaded to repository

The following files have not been uploaded to the repository due to size restrictions from GitHub.

1. The original reads - links to these can be found here:
    2. 
    3. 
2. The BAM files.
3. The genome referece file.
4. The 
