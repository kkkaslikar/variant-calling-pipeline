import os

#%%
def list_directories(dirpath="./", givepath = False):
    '''For a given directory provided as a pathname, list the names of subdirectories directly under it, as just names or as complete path names. Provided as a list.'''
    
    if givepath == False:
        dirnames = []
        for entry in os.scandir(dirpath):
            if os.path.isdir(os.path.join(dirpath, entry.name)):
                dirnames.append(entry.name)
        return dirnames
            
    elif givepath == True:
        dirpathnames = []
        for entry in os.scandir(dirpath):
            if os.path.isdir(os.path.join(dirpath, entry.name)):
                dirpathnames.append(os.path.join(dirpath, entry.name))
        return dirpathnames

def list_files(dirpath="./", givepath = False):
    '''For a given directory provided as a pathname, list the names of file directly under it, as just names or as complete path names. Provided as a list.'''
    
    if givepath == False:
        filenames = []
        for entry in os.scandir(dirpath):
            if os.path.isfile(os.path.join(dirpath, entry.name)):
                filenames.append(entry.name)
        return filenames
            
    elif givepath == True:
        filepathnames = []
        for entry in os.scandir(dirpath):
            if os.path.isfile(os.path.join(dirpath, entry.name)):
                filepathnames.append(os.path.join(dirpath, entry.name))
        return filepathnames

#%%

genome="Homo_sapiens.GRCh38.dna.primary_assembly.fa"

samples = [f.strip(".sorted.bam") for f in list_files() if f.endswith("sorted.bam")]

tools = ["bcftools", "freebayes"]

rule all:
  input: "consensus"
  


rule get_pileup:
  input:
    gen = genome,
    aln = expand("{sample}.sorted.bam", sample = samples)
  output: expand("{sample}.bcf", sample = samples)
  threads: 4
  shell:
    "bcftools mpileup -f {input.gen} {input.aln} --threads {threads} -q 20- -O b -o {output}"


rule call_var_bcftools:
  input: expand("{sample}.bcf" , sample = samples)
  output: expand("{sample}_bcftools_called.bcf", sample = samples)
  shell:
    "bcftools call -v -m -O b -o {output} {input}"


rule call_var_freebayes:
  input:
    gen=genome,
    aln = expand("{sample}.sorted.bam" , sample = samples),
  output: expand("{sample}_freebayes_called.vcf", sample = samples)
  shell:
    "freebayes -f {input.gen} -q 20 {input.aln} > {output}"
    

rule filter_bcftools_calls:
  input: expand("{sample}_bcftools_called.bcf" , sample = samples)
  output: expand("{sample}_bcftools_called_filtered.vcf", sample = samples)
  shell:
    "bcftools filter  -i 'DP>5 & QUAL>20 & TYPE=\"snp\"' -O v -o {output} {input}"


rule filter_freebayes_calls:
  input: expand("{sample}_freebayes_called.vcf" , sample = samples)
  output: expand("{sample}_freebayes_called_filtered.vcf", sample = samples)
  shell:
    "bcftools filter  -i 'FMT/DP>5 & QUAL>20 & TYPE=\"snp\"' -O v -o {output} {input}"


rule annotate_freebayes:
  input: expand("{sample}_freebayes_called_filtered.vcf", sample = samples)
  output:
    main=expand("{sample}_freebayes_snpEff.vcf", sample = samples),
    report=expand("{sample}_freebayes_snpEff.html", sample = samples)
  shell:
    "snpEff eff -s {output.report} hg38 {input} > {output.main}"

rule annotate_bcftools:
  input: expand("{sample}_bcftools_called_filtered.vcf", sample = samples)
  output:
    main=expand("{sample}_bcftools_snpEff.vcf", sample = samples),
    report=expand("{sample}_bcftools_snpEff.html", sample = samples)
  shell:
    "snpEff eff -s {output.report} hg38 {input} > {output.main}"

rule convert_index:
  input:
    bcf=expand("{sample}_bcftools_snpEff.vcf", sample = samples),
    fb=expand("{sample}_freebayes_snpEff.vcf", sample = samples)
  output:
    bcf=expand("{sample}_bcftools_snpEff.vcf.gz", sample = samples),
    fb=expand("{sample}_freebayes_snpEff.vcf.gz", sample = samples)
  shell:
    "bcftools view {input.bcf} -O z -o {output.bcf}; bcftools index {output.bcf}; bcftools view {input.fb} -O z -o {output.fb}; bcftools index {output.fb}"

rule intersect:
  input:
    bcf=expand("{sample}_bcftools_snpEff.vcf.gz", sample = samples),
    fb=expand("{sample}_freebayes_snpEff.vcf.gz", sample = samples)
  output: directory("consensus")
  shell:
    "bcftools isec -c snps {input.fb} {input.bcf} -O b -p {output}"
