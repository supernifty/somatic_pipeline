
configfile: "cfg/config.yaml"
cluster = json.load(open("cfg/cluster.json"))

# NOTE: no mitochondria MT because they aren't in our exome
GATK_CHROMOSOMES=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')

### helper functions ###
def read_group(wildcards):
  '''
    determine read group from sample name
  '''
  # in/0656045001_BC_H5CNYDSXX_GTGTTCTA_L004_R1.fastq.gz
  # 0757079003_T_A_H5CNYDSXX_ACGTATCA_L004_R1.fastq.gz
  suffix = config["samples"][wildcards.sample][0].replace("in/{}_".format(wildcards.sample), "") # in/S1_RG_2_R1.fastq.gz
  fields = suffix.split("_") # H5CNYDSXX_GTGTTCTA_L004_R1 
  flowcell = fields[0]
  barcode = fields[1]
  lane = fields[2]
  return "@RG\tID:{sample}.{flowcell}.{barcode}.{lane}\tSM:{sample}\tPU:{flowcell}.{barcode}.{lane}\tPL:Illumina".format(flowcell=flowcell, sample=wildcards.sample, lane=lane, barcode=barcode)

def tumour_germline_dup_bams(wildcards):
  tumour_bam = 'out/{}.sorted.dups.bam'.format(wildcards.tumour)
  normal_bam = 'out/{}.sorted.dups.bam'.format(config["tumours"][wildcards.tumour])
  return [tumour_bam, normal_bam]

def tumour_germline_bams(wildcards):
  tumour_bam = 'out/{}.sorted.dups.bam'.format(wildcards.tumour)
  normal_bam = 'out/{}.sorted.dups.bam'.format(config["tumours"][wildcards.tumour])
  return [tumour_bam, normal_bam]

def germline_samples():
  samples = set(config['samples'])
  tumours = set(config['tumours'])
  return list(samples.difference(tumours))

def germline_sample(wildcards):
  return config["tumours"][wildcards.tumour]

### final outputs ###
rule all:
  input:
    expand("out/{tumour}.strelka.somatic.snvs.af.vep.vcf.gz", tumour=config['tumours']), # somatic snvs strelka
    expand("out/{tumour}.strelka.somatic.indels.vep.vcf.gz", tumour=config['tumours']), # somatic indels strelka
    expand("out/{germline}.strelka.germline.filter_gt.vep.vcf.gz", germline=germline_samples()), # germline strelka calls
    expand("out/{germline}.strelka.germline.filter_gt.vep.vcf.gz", germline=config['tumours']), # run the tumours as if they were germline
    expand("out/{sample}.oxo_metrics.txt", sample=config['samples']),
    expand("out/{sample}.artifact_metrics.txt.error_summary_metrics", sample=config['samples']),
    expand("out/{tumour}.strelka.somatic.snvs.bias.vcf.gz", tumour=config['tumours']),
    expand("out/{tumour}.mutect2.filter.bias.vcf.gz", tumour=config['tumours']), # somatic mutect2 with dkfz bias annotation
    expand("out/fastqc/{sample}/completed", sample=config['samples']), # fastqc
    expand("out/{sample}.metrics.insertsize", sample=config['samples']),
    expand("out/{sample}.metrics.alignment", sample=config['samples']),
    expand("out/{sample}.metrics.target", sample=config['samples']),
    expand("out/{tumour}.verifybamid.somatic.completed", tumour=config['tumours']),
    expand("out/{tumour}.mutect2.filter.vep.vcf.gz", tumour=config['tumours']),
    expand("out/{tumour}.intersect.filter.vcf.gz", tumour=config['tumours']), # somatic combination of strelka and mutect2, with filter

    expand("out/{tumour}.mutect2.filter.genes_of_interest.tsv", tumour=config['tumours']), # filter on genes of interest
    expand("out/{tumour}.varscan.copynumber.deletions.bed", tumour=config['tumours']), # TODO fails on mini sample
    expand("out/{tumour}.platypus.somatic.vcf.gz", tumour=config['tumours']), # platypus somatic calls
    expand("out/{tumour}.loh.bed", tumour=config['tumours']), # loh regions

    "out/mutational_signatures.combined",
    "out/mutational_signatures.filter.combined",
    "out/max_coverage.tsv",
    "out/max_trimmed_coverage.tsv",
    "out/ontarget.tsv",
    "out/germline_joint.hc.normalized.vep.vcf.gz", # gatk calls for all germline samples
    "out/tumour_joint.hc.normalized.vep.vcf.gz", # gatk calls for all tumour samples (as germline)
    "out/qc.summary.tsv",
    "out/multiqc.html", # overall general qc
    "out/ontarget.png", # combined ontarget coverage plots
    "out/ontarget_tumour.png", # somatic ontarget coverage plots
    "out/ontarget_germline.png", # germline ontarget coverage plots
    "out/mutation_rate.tsv",

### QC ###

rule qc_summary:
  input:
    expand("out/{sample}.artifact_metrics.txt.error_summary_metrics", sample=config['samples'])
  output:
    "out/qc.summary.tsv"
  log:
    stderr="log/make_summary.stderr"
  shell:
    "python src/make_summary.py --verbose --samples {input} > {output} 2>{log.stderr}"

rule fastqc:
  input:
    fastqs=lambda wildcards: config["samples"][wildcards.sample]
  output:
    "out/fastqc/{sample}/completed"
  shell:
    "{config[module_java]} && "
    "mkdir -p out/fastqc/{wildcards.sample} && "
    "tools/FastQC/fastqc --extract --outdir out/fastqc/{wildcards.sample} {input.fastqs} && "
    "touch {output}"

rule make_sequence_dict:
  input:
    reference=config["genome"]
  output:
    config["genome_dict"]
  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar CreateSequenceDictionary REFERENCE={input.reference} OUTPUT={output}"

rule make_intervals:
  input:
    bed=config["regions"],
    dict=config["genome_dict"]
  output:
    "out/regions.intervals"
  shell:  
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar BedToIntervalList INPUT={input.bed} OUTPUT={output} SEQUENCE_DICTIONARY={input.dict}"

rule qc_target:
  input:
    reference=config["genome"],
    bam="out/{sample}.sorted.dups.bam",
    intervals="out/regions.intervals"
  output:
    "out/{sample}.metrics.target"
  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar CollectHsMetrics REFERENCE_SEQUENCE={input.reference} INPUT={input.bam} OUTPUT={output} BAIT_INTERVALS={input.intervals} TARGET_INTERVALS={input.intervals}"

rule qc_alignment:
  input:
    reference=config["genome"],
    bam="out/{sample}.sorted.dups.bam"
  output:
    "out/{sample}.metrics.alignment"
  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE={input.reference} INPUT={input.bam} OUTPUT={output}"

rule qc_insertsize:
  input:
    bam="out/{sample}.sorted.dups.bam"
  output:
    "out/{sample}.metrics.insertsize"
  shell:
    "{config[module_java]} && "
    "{config[module_R]} && "
    "java -jar tools/picard-2.8.2.jar CollectInsertSizeMetrics INPUT={input.bam} OUTPUT={output} HISTOGRAM_FILE={output}.pdf"

rule qc_conpair:
  input:
    reference=config["genome"],
    reference_dict=config["genome_dict"],
    bams=tumour_germline_bams
  output:
    "out/{tumour}.concordance",
    "out/{tumour}.contamination"
  log:
    stderr="log/{tumour}.conpair.stderr",
    stdout="log/{tumour}.conpair.stdout"
  shell:
    "( "
    "{config[module_java]} && "
    "mkdir -p tmp/conpair_$$ && "
    "python tools/Conpair/scripts/run_gatk_pileup_for_sample.py --reference {input.reference} --conpair_dir tools/Conpair --gatk tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -B {input.bams[0]} -O tmp/conpair_$$/tumour.pileup && "
    "python tools/Conpair/scripts/run_gatk_pileup_for_sample.py --reference {input.reference} --conpair_dir tools/Conpair --gatk tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -B {input.bams[1]} -O tmp/conpair_$$/normal.pileup && "
    "PYTHONPATH=tools/Conpair/modules CONPAIR_DIR=tools/Conpair python tools/Conpair/scripts/verify_concordance.py -T tmp/conpair_$$/tumour.pileup -N tmp/conpair_$$/normal.pileup --outfile {output[0]} --normal_homozygous_markers_only && "
    "PYTHONPATH=tools/Conpair/modules CONPAIR_DIR=tools/Conpair python tools/Conpair/scripts/estimate_tumor_normal_contamination.py -T tmp/conpair_$$/tumour.pileup -N tmp/conpair_$$/normal.pileup --outfile {output[1]} && "
    "rm -r tmp/conpair_$$ "
    ") 1>{log.stdout} 2>{log.stderr}"

rule qc_verifybamid_tumour:
  input:
    vcf="out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
    bam="out/{tumour}.sorted.dups.bam",
    bai="out/{tumour}.sorted.dups.bai",
  output:
    "out/{tumour}.verifybamid.somatic.completed"
  log:
    stderr="log/{tumour}.verifybamid.stderr"
  shell:
    "tools/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID --vcf {input.vcf} --bam {input.bam} --bai {input.bai} --out out/{wildcards.tumour}.verifybamid --verbose 2>{log.stderr} && touch {output}"

# verifybamid doesn't like strelka's germline output
#rule qc_verifybamid_germline:
#  input:
#    vcf="out/{germline}.strelka.germline.filter_gt.vcf.gz",
#    bam="out/{germline}.sorted.dups.bam",
#    bai="out/{germline}.sorted.dups.bai",
#  output:
#    "out/{germline}.verifybamid.germline.completed"
#  log:
#    stderr="log/{germline}.verifybamid.stderr"
#  shell:
#    "tools/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID --vcf {input.vcf} --bam {input.bam} --bai {input.bai} --out out/{wildcards.germline}.verifybamid --verbose 2>{log.stderr} && touch {output}"

rule qc_depth_of_coverage:
  input:
    reference=config["genome"],
    bed=config["regions"],
    bam="out/{sample}.sorted.dups.bam"
  output:
    "out/{sample}.depth_of_coverage.sample_summary"
  log:
    "log/{sample}.depth_of_coverage.stderr"
  params:
    prefix="out/{sample}.depth_of_coverage"
  shell:
    "{config[module_java]} && "
    "java -jar tools/GenomeAnalysisTK-3.7.0.jar -T DepthOfCoverage -R {input.reference} -o {params.prefix} -I {input.bam} -L {input.bed} && rm {params.prefix} "
    "2>{log}"

rule multiqc:
  input:
    expand("out/fastqc/{sample}/completed", sample=config['samples']),
    expand("out/{sample}.metrics.alignment", sample=config['samples']),
    expand("out/{sample}.metrics.insertsize", sample=config['samples']),
    expand("out/{sample}.metrics.target", sample=config['samples']),
    expand("out/{tumour}.concordance", tumour=config['tumours']), # TODO
    expand("out/{tumour}.contamination", tumour=config['tumours']), # TODO
    expand("out/{tumour}.verifybamid.somatic.completed", tumour=config['tumours']),
    expand("out/{sample}.depth_of_coverage.sample_summary", sample=config['samples']),
    "out/qc.summary.tsv"
    
  output:
    "out/multiqc.html"
  shell:
    "multiqc --force --filename {output} out"

rule qc_on_target_coverage_hist:
  input:
    bed=config["regions"],
    bam="out/{sample}.sorted.dups.bam"
  output:
    hist="out/{sample}.ontarget.hist",
  shell:
    "{config[module_bedtools]} && "
    #"bedtools sort -g reference/genome.lengths -i {input.bed} | bedtools merge -i - | bedtools coverage -sorted -hist -b {input.bam} -a stdin -g reference/genome.lengths | grep ^all > {output.hist}" # 2.27
    "bedtools sort -faidx reference/genome.lengths -i {input.bed} | bedtools merge -i - | bedtools coverage -sorted -hist -b {input.bam} -a stdin -g reference/genome.lengths | grep ^all > {output.hist}"

rule qc_on_target_coverage:
  input:
    reference=config["genome"],
    bed=config["regions"],
    bam="out/{sample}.sorted.dups.bam"
  output:
    summary="out/{sample}.ontarget.summary"
  shell:
    "{config[module_bedtools]} && "
    #"bedtools sort -g reference/genome.lengths -i {input.bed} | bedtools merge -i - | bedtools coverage -sorted -a stdin -b {input.bam} -d -g reference/genome.lengths | cut -f5 | src/stats.py > {output.summary}" # 2.27
    "bedtools sort -faidx reference/genome.lengths -i {input.bed} | bedtools merge -i - | bedtools coverage -sorted -a stdin -b {input.bam} -d -g reference/genome.lengths | cut -f5 | src/stats.py > {output.summary}"

rule qc_on_target_coverage_plot:
  input:
    expand("out/{sample}.ontarget.hist", sample=config['samples']),
  output:
    "out/ontarget.png"
  shell:
    "src/plot_coverage.py --target {output} --files {input} --max 10000"

rule qc_on_target_coverage_plot_germline:
  input:
    expand("out/{germline}.ontarget.hist", germline=germline_samples()),
  output:
    "out/ontarget_germline.png"
  shell:
    "src/plot_coverage.py --target {output} --files {input} --max 10000"

rule qc_on_target_coverage_plot_tumour:
  input:
    expand("out/{tumour}.ontarget.hist", tumour=config['tumours']),
  output:
    "out/ontarget_tumour.png"
  shell:
    "src/plot_coverage.py --target {output} --files {input} --max 10000"

rule qc_on_target_coverage_combined:
  input:
    expand("out/{sample}.ontarget.summary", sample=config['samples']),
  output:
    "out/ontarget.tsv"
  shell:
    "echo \"Filename    n       Mean    Min     Max     Total\" >{output} && "
    "for f in {input}; do echo \"$f     $(tail -1 $f)\" >> {output}; done"

### alignment ###
rule trim:
  input:
    fastqs=lambda wildcards: config["samples"][wildcards.sample]
  output:
    "out/{sample}_R1.trimmed.paired.fq.gz",
    "out/{sample}_R1.trimmed.unpaired.fq.gz",
    "out/{sample}_R2.trimmed.paired.fq.gz",
    "out/{sample}_R2.trimmed.unpaired.fq.gz"
  log:
    stderr="out/{sample}.trimmomatic.stderr" # has useful output
  params:
    cores=cluster["align"]["n"],
  shell:
    "{config[module_java]} && "
    "java -jar tools/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads {params.cores} -phred33 {input.fastqs} {output} ILLUMINACLIP:tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:1:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>{log.stderr}"

rule align:
  input:
    reference=config["genome"],
    fastq_r1="out/{sample}_R1.trimmed.paired.fq.gz", 
    fastq_r2="out/{sample}_R2.trimmed.paired.fq.gz"
    #fastqs=lambda wildcards: config["samples"][wildcards.sample]

  output:
    "tmp/{sample}.paired.bam"

  log:
    "log/{sample}.paired.bwa.log"

  params:
    cores=cluster["align"]["n"],
    read_group=read_group

  shell:
    "{config[module_bwa]} && {config[module_samtools]} && "
    "(bwa mem -M -t {params.cores} -R \"{params.read_group}\" {input.reference} {input.fastq_r1} {input.fastq_r2} | samtools view -b -h -o {output} -) 2>{log}"

rule align_unpaired:
  input:
    reference=config["genome"],
    fastq_r1="out/{sample}_R1.trimmed.unpaired.fq.gz", 
    fastq_r2="out/{sample}_R2.trimmed.unpaired.fq.gz"

  output:
    r1="tmp/{sample}_R1.unpaired.bam",
    r2="tmp/{sample}_R2.unpaired.bam"

  log:
    "log/{sample}.bwa_unpaired.log"

  params:
    cores=cluster["align_unpaired"]["n"],
    read_group=read_group

  shell:
    "{config[module_bwa]} && {config[module_samtools]} && "
    "(bwa mem -M -t {params.cores} -R \"{params.read_group}\" {input.reference} {input.fastq_r1} | samtools view -b -h -o {output.r1} - && "
    "bwa mem -M -t {params.cores} -R \"{params.read_group}\" {input.reference} {input.fastq_r2} | samtools view -b -h -o {output.r2} -) 2>{log}"

# sort the bam
rule merge_bams:
  input:
    "tmp/{sample}.paired.bam",
    "tmp/{sample}_R1.unpaired.bam",
    "tmp/{sample}_R2.unpaired.bam"
  output:
    "tmp/{sample}.merged.bam"
  log:
    "log/{sample}.merge.log"
  shell:
    "{config[module_samtools]} && "
    "samtools merge {output} {input} 2>{log}"

# sort the bam
rule sort:
  input:
    "tmp/{sample}.merged.bam"

  output:
    bam="tmp/{sample}.sorted.bam",
    bai="tmp/{sample}.sorted.bai"

  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar SortSam INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=2000000 CREATE_INDEX=True"

# duplicates
rule gatk_duplicates:
  input:
    "tmp/{sample}.sorted.bam"
  output:
    "out/{sample}.sorted.dups.bam",
    "out/{sample}.sorted.dups.bai",
    "out/{sample}.markduplicates.metrics"
  log:
    "log/{sample}.markduplicates.stderr"
  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar MarkDuplicates INPUT={input} OUTPUT={output[0]} METRICS_FILE={output[2]} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True CREATE_INDEX=True MAX_RECORDS_IN_RAM=2000000"

### germline variant calling ###

rule gatk_haplotype_caller:
  input:
    bam="out/{germline}.sorted.dups.bam",
    reference=config["genome"]
  output:
    recal="out/{germline}.recal_table",
    bqsr="out/{germline}.sorted.dups.bqsr.bam",
    gvcf="out/{germline}.hc.gvcf.gz"
  log:
    "log/{germline}.hc.log"
  shell:
    "({config[module_java]} && "
    "tools/gatk-4.0.0.0/gatk BaseRecalibrator --input {input.bam} --output {output.recal} -R {input.reference} --known-sites reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz --known-sites reference/gatk-4-bundle-b37/Mills_and_1000G_gold_standard.indels.b37.vcf.bgz --known-sites reference/gatk-4-bundle-b37/1000G_phase1.indels.b37.vcf.bgz && "
    "tools/gatk-4.0.0.0/gatk ApplyBQSR -R {input.reference} -I {input.bam} -bqsr {output.recal} -O {output.bqsr} && "
    "tools/gatk-4.0.0.0/gatk HaplotypeCaller -R {input.reference} -I {output.bqsr} --emit-ref-confidence GVCF --dbsnp reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz -O {output.gvcf}"
    ") 2>{log}"

rule gatk_joint_genotype:
  input:
    gvcfs=expand("out/{germline}.hc.gvcf.gz", germline=germline_samples()),
    reference=config["genome"]
  output:
    "out/germline_joint_{chromosome}.vcf"
  log:
    "log/gatk_joint_{chromosome}.stderr"
  params:
    variant_list=' '.join(['--variant {}'.format(gvcf) for gvcf in expand("out/{germline}.hc.gvcf.gz", germline=germline_samples())])
  shell:
    "({config[module_java]} && "
    "java -jar tools/GenomeAnalysisTK-3.7.0.jar -T CombineGVCFs -R {input.reference} {params.variant_list} -L {wildcards.chromosome} -o tmp/germline_combined_{wildcards.chromosome}.gvcf && "
    "tools/gatk-4.0.0.0/gatk GenotypeGVCFs -R {input.reference} --dbsnp reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz -V tmp/germline_combined_{wildcards.chromosome}.gvcf -L {wildcards.chromosome} --use-new-qual-calculator true --output out/germline_joint_{wildcards.chromosome}.vcf"
    ") 2>{log}"

# call germline variants on the tumour for validation
rule gatk_joint_genotype_tumours:
  input:
    gvcfs=expand("out/{germline}.hc.gvcf.gz", germline=config["tumours"]),
    reference=config["genome"]
  output:
    "out/tumour_joint_{chromosome}.vcf"
  log:
    "log/gatk_joint_tumour_{chromosome}.stderr"
  params:
    variant_list=' '.join(['--variant {}'.format(gvcf) for gvcf in expand("out/{germline}.hc.gvcf.gz", germline=config["tumours"])])
  shell:
    "({config[module_java]} && "
    "java -jar tools/GenomeAnalysisTK-3.7.0.jar -T CombineGVCFs -R {input.reference} {params.variant_list} -L {wildcards.chromosome} -o tmp/tumours_combined_{wildcards.chromosome}.gvcf && "
    "tools/gatk-4.0.0.0/gatk GenotypeGVCFs -R {input.reference} --dbsnp reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz -V tmp/tumours_combined_{wildcards.chromosome}.gvcf -L {wildcards.chromosome} --use-new-qual-calculator true --output out/tumour_joint_{wildcards.chromosome}.vcf"
    ") 2>{log}"

# notes: 
#   VariantRecalibrator removed due to large cohort size requirements
rule gatk_post_genotype:
  input:
    gvcfs=expand("out/germline_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES),
    reference=config["genome"]
  output:
    "out/germline_joint.hc.normalized.vcf"
  log:
    "log/gatk_post.stderr"
  params:
    inputs=' '.join(['--INPUT={}'.format(gvcf) for gvcf in expand("out/germline_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES)])
  shell:
    "({config[module_java]} && {config[module_R]} && {config[module_samtools]} && "
    "{config[module_htslib]} && "
    "tools/gatk-4.0.0.0/gatk GatherVcfs -R {input.reference} --OUTPUT=tmp/germline_joint.vcf {params.inputs} && "
    "bgzip -c < tmp/germline_joint.vcf > tmp/germline_joint.vcf.bgz && tabix -p vcf tmp/germline_joint.vcf.bgz && "
    "tools/gatk-4.0.0.0/gatk CalculateGenotypePosteriors -R {input.reference} --supporting reference/gatk-4-bundle-b37/1000G_phase3_v4_20130502.sites.vcf.bgz -V tmp/germline_joint.vcf.bgz -O tmp/germline_joint.cgp.vcf && "
    "tools/vt-0.577/vt normalize -n -r {input.reference} tmp/germline_joint.cgp.vcf -o tmp/germline_joint.cgp.normalized.vcf && "
    "tools/vt-0.577/vt decompose -s tmp/germline_joint.cgp.normalized.vcf | tools/vt-0.577/vt normalize -r {input.reference} - -o {output}"
    ") 2>{log}"

rule gatk_post_genotype_tumours:
  input:
    gvcfs=expand("out/tumour_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES),
    reference=config["genome"]
  output:
    "out/tumour_joint.hc.normalized.vcf"
  log:
    "log/gatk_post.stderr"
  params:
    inputs=' '.join(['--INPUT={}'.format(gvcf) for gvcf in expand("out/tumour_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES)])
  shell:
    "({config[module_java]} && {config[module_R]} && {config[module_samtools]} && "
    "{config[module_htslib]} && "
    "tools/gatk-4.0.0.0/gatk GatherVcfs -R {input.reference} --OUTPUT=tmp/tumour_joint.vcf {params.inputs} && "
    "bgzip -c < tmp/tumour_joint.vcf > tmp/tumour_joint.vcf.bgz && tabix -p vcf tmp/tumour_joint.vcf.bgz && "
    "tools/gatk-4.0.0.0/gatk CalculateGenotypePosteriors -R {input.reference} --supporting reference/gatk-4-bundle-b37/1000G_phase3_v4_20130502.sites.vcf.bgz -V tmp/tumour_joint.vcf.bgz -O tmp/tumour_joint.cgp.vcf && "
    "tools/vt-0.577/vt normalize -n -r {input.reference} tmp/tumour_joint.cgp.vcf -o tmp/tumour_joint.cgp.normalized.vcf && "
    "tools/vt-0.577/vt decompose -s tmp/tumour_joint.cgp.normalized.vcf | tools/vt-0.577/vt normalize -r {input.reference} - -o {output}"
    ") 2>{log}"

### qc ###
rule qc_max_coverage_combine:
  input:
    expand("out/{sample}.max_coverage", sample=config['samples']),
  output:
    "out/max_coverage.tsv"
  shell:
    ">{output} && "
    "for f in {input}; do echo \"$f $(grep \"Max coverage\" $f)\" >> {output}; done"

rule qc_max_trimmed_coverage_combine:
  input:
    expand("out/{sample}.max_trimmed_coverage", sample=config['samples']),
  output:
    "out/max_trimmed_coverage.tsv"
  shell:
    ">{output} && "
    "for f in {input}; do echo \"$f $(grep \"Max coverage\" $f)\" >> {output}; done"

rule qc_max_coverage:
  input:
    bed=config["regions"],
    fastqs=lambda wildcards: config["samples"][wildcards.sample]
  output:
    "out/{sample}.max_coverage"
  log:
    "log/{sample}.max_coverage.stderr"
  shell:
    "src/max_coverage.py --verbose --bed {input.bed} --fastqs {input.fastqs} >{output} 2>{log}"

rule qc_max_trimmed_coverage:
  input:
    bed=config["regions"],
    fastqs=("out/{sample}_R1.trimmed.paired.fq.gz", "out/{sample}_R1.trimmed.unpaired.fq.gz", "out/{sample}_R2.trimmed.paired.fq.gz", "out/{sample}_R2.trimmed.unpaired.fq.gz")

  output:
    "out/{sample}.max_trimmed_coverage"
  log:
    "log/{sample}.max_trimmed_coverage.stderr"
  shell:
    "src/max_coverage.py --verbose --bed {input.bed} --fastqs {input.fastqs} >{output} 2>{log}"

rule qc_sequencing_artifacts:
  input:
    bam="out/{sample}.sorted.dups.bam",
    reference=config["genome"]

  output:
    "out/{sample}.artifact_metrics.txt.error_summary_metrics"

  params:
    prefix="out/{sample}.artifact_metrics.txt"

  log:
    stderr="log/{sample}.artifact.err",
    stdout="log/{sample}.artifact.out"

  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar CollectSequencingArtifactMetrics I={input.bam} O={params.prefix} R={input.reference} 2>{log.stderr} 1>{log.stdout}"

rule qc_oxidative_artifacts:
  input:
    bam="out/{sample}.sorted.dups.bam",
    reference=config["genome"]

  output:
    "out/{sample}.oxo_metrics.txt"

  log:
    stderr="log/{sample}.oxo.err",
    stdout="log/{sample}.oxo.out"

  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar CollectOxoGMetrics I={input.bam} O={output} R={input.reference} 2>{log.stderr} 1>{log.stdout}"

### somatic variant calling ###
rule strelka_somatic:
  input:
    reference=config["genome"],
    bams=tumour_germline_bams

  output:
    "out/{tumour}.strelka.somatic.snvs.vcf.gz",
    "out/{tumour}.strelka.somatic.indels.vcf.gz",

  log:
    "log/{tumour}.strelka.somatic.log"

  params:
    cores=cluster["strelka_somatic"]["n"]

  shell:
    "(mkdir -p tmp/strelka_{wildcards.tumour}_$$ && "
    "tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py "
    "--ref {input.reference} "
    "--tumorBam {input.bams[0]} "
    "--normalBam {input.bams[1]} "
    "--runDir tmp/strelka_{wildcards.tumour}_$$ "
    "--exome && "
    "tmp/strelka_{wildcards.tumour}_$$/runWorkflow.py -m local -j {params.cores} && "
    "mv tmp/strelka_{wildcards.tumour}_$$/results/variants/somatic.snvs.vcf.gz {output[0]} && " 
    "mv tmp/strelka_{wildcards.tumour}_$$/results/variants/somatic.indels.vcf.gz {output[1]} && " 
    "rm -r tmp/strelka_{wildcards.tumour}_$$ ) 2>{log}"

rule annotate_af_somatic:
  input:
    "out/{tumour}.strelka.somatic.snvs.vcf.gz",
  output:
    "out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
  log:
    stderr="log/{tumour}.annotate_af.stderr"
  shell:
    "{config[module_samtools]} && "
    "{config[module_htslib]} && "
    "src/annotate_af.py TUMOR {input} | bgzip >{output} 2>{log.stderr}"

# tumour only for each germline
rule mutect2_sample_pon:
  input:
    reference=config["genome"],
    bam="out/{germline}.sorted.dups.bam",
    regions="reference/regions.bed"
  output:
    "out/{germline}.mutect2.pon.vcf.gz",
  log:
    stderr="log/{germline}.mutect2.pon.stderr"
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.0.0.0/gatk Mutect2 -R {input.reference} -I {input.bam} --tumor-sample {wildcards.germline} -L {input.regions} -O {output} --interval-padding 1000 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 2>{log.stderr}"

# combine all the samples to make an overall pon
rule mutect2_pon:
  input:
    vcfs=expand("out/{germline}.mutect2.pon.vcf.gz", germline=germline_samples()),
  output:
    "out/mutect2.pon.vcf.gz"
  log:
    stderr="log/mutect2.pon.stderr"
  params:
    vcfs=' '.join(['--vcfs {}'.format(vcf) for vcf in expand("out/{germline}.mutect2.pon.vcf.gz", germline=germline_samples())])
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.0.0.0/gatk CreateSomaticPanelOfNormals {params.vcfs} -O {output} 2>{log.stderr}"

# mutect2 somatic calls
rule mutect2_somatic_chr:
  input:
    reference=config["genome"],
    dbsnp="reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz",
    regions="reference/regions.bed",
    pon="out/mutect2.pon.vcf.gz",
    gnomad="reference/af-only-gnomad.raw.sites.b37.vcf.gz",
    bams=tumour_germline_dup_bams
  output:
    "tmp/{tumour}.{chromosome}.mutect2.vcf.gz",
  log:
    stderr="log/{tumour}.{chromosome}.mutect2.stderr"
  params:
    germline=lambda wildcards: config["tumours"][wildcards.tumour]
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.0.0.0/gatk --java-options '-Xmx30G' Mutect2 -R {input.reference} -I {input.bams[0]} -I {input.bams[1]} --tumor-sample {wildcards.tumour} --normal-sample {params.germline} --output {output} --output-mode EMIT_VARIANTS_ONLY --dbsnp {input.dbsnp} --germline-resource {input.gnomad} --af-of-alleles-not-in-resource 0.0000025 -pon {input.pon} --interval-padding 1000 -L {input.regions} -L {wildcards.chromosome} --interval-set-rule INTERSECTION --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter"

rule mutect2_somatic:
  input:
    vcfs=expand("tmp/{{tumour}}.{chromosome}.mutect2.vcf.gz", chromosome=GATK_CHROMOSOMES)
  output:
    "out/{tumour}.mutect2.vcf.gz"
  log:
    stderr="log/{tumour}.mutect2.mergevcfs.stderr"
  params:
    inputs=' '.join(['I={}'.format(vcf) for vcf in expand("tmp/{{tumour}}.{chromosome}.mutect2.vcf.gz", chromosome=GATK_CHROMOSOMES)])
  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar MergeVcfs {params.inputs} O={output} 2>{log.stderr}"

rule mutect2_filter:
  input:
    vcf="out/{tumour}.mutect2.vcf.gz",
    bam="out/{tumour}.sorted.dups.bam",
    gnomad="reference/af-only-gnomad.raw.sites.b37.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.vcf.gz"
  log:
    stderr="log/{tumour}.mutect2-filter.stderr",
    stdout="log/{tumour}.mutect2-filter.stdout"
  shell:
    "({config[module_java]} && "
    "tools/gatk-4.0.0.0/gatk GetPileupSummaries -I {input.bam} -V {input.gnomad} -O tmp/{wildcards.tumour}.mutect2.pileup.table && "
    "tools/gatk-4.0.0.0/gatk CalculateContamination -I tmp/{wildcards.tumour}.mutect2.pileup.table -O tmp/{wildcards.tumour}.mutect2.contamination.table && "
    "tools/gatk-4.0.0.0/gatk FilterMutectCalls -V {input.vcf} --contamination-table tmp/{wildcards.tumour}.mutect2.contamination.table -O {output}) 1>{log.stdout} 2>{log.stderr}"

### germline variant calling ###
rule filter_strelka_germline:
  input:
    "out/{germline}.strelka.germline.vcf.gz",
  output:
    "out/{germline}.strelka.germline.filter_gt.vcf.gz",
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input} | grep -v NoPassedVariantGTs | bgzip > {output}"
    
rule strelka_germline:
  input:
    reference=config["genome"],
    bam="out/{germline}.sorted.dups.bam"

  output:
    "out/{germline}.strelka.germline.vcf.gz",
    "out/{germline}.strelka.germline.vcf.gz.tbi"

  log:
    "log/{germline}.strelka.germline.log"

  params:
    cores=cluster["strelka_germline"]["n"]

  shell:
    "(mkdir -p tmp/strelka_{wildcards.germline}_$$ && "
    "tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py "
    "--referenceFasta {input.reference} "
    "--bam {input.bam} "
    "--runDir tmp/strelka_{wildcards.germline}_$$ "
    "--exome && "
    "tmp/strelka_{wildcards.germline}_$$/runWorkflow.py -m local -j {params.cores} && "
    "mv tmp/strelka_{wildcards.germline}_$$/results/variants/variants.vcf.gz {output[0]} && " 
    "mv tmp/strelka_{wildcards.germline}_$$/results/variants/variants.vcf.gz.tbi {output[1]} && " 
    "rm -r tmp/strelka_{wildcards.germline}_$$ ) 2>{log}"

### platypus ###
rule platypus_somatic:
  input:
    reference=config["genome"],
    bams=tumour_germline_bams

  output:
    joint="out/{tumour}.platypus.joint.vcf.gz",
    somatic="out/{tumour}.platypus.somatic.vcf.gz"

  log:
    "log/{tumour}.platypus.somatic.log"

  params:
    germline=germline_sample

  shell:
    # platypus has to run from build directory
    "({config[module_python2]} && "
    "{config[module_htslib]} && "
    "{config[module_samtools]} && "
    "tools/Platypus_0.8.1/Platypus.py callVariants --bamFiles={input.bams[0]},{input.bams[1]} --refFile={input.reference} --output=tmp/platypus_{wildcards.tumour}.vcf && "
    "bgzip < tmp/platypus_{wildcards.tumour}.vcf > {output.joint} && "
    "python tools/Platypus/extensions/Cancer/somaticMutationDetector.py --inputVCF tmp/platypus_{wildcards.tumour}.vcf --outputVCF {output.somatic} --tumourSample {wildcards.tumour}.sorted.dups --normalSample {params.germline}.sorted.dups) 2>{log}"

### pindel ###
#rule pindel_somatic:


### annotation ###
rule annotate_vep_somatic_snvs:
  input:
    vcf="out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
    reference=config['genome']
  output:
    "out/{tumour}.strelka.somatic.snvs.af.vep.vcf.gz"
  log:
    "log/{tumour}.vep.log"
  params:
    cores=cluster["annotate_vep_somatic_snvs"]["n"]
  shell:
    "{config[module_samtools]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

rule annotate_vep_somatic_indels:
  input:
    vcf="out/{tumour}.strelka.somatic.indels.vcf.gz",
    reference=config['genome']
  output:
    "out/{tumour}.strelka.somatic.indels.vep.vcf.gz"
  log:
    "log/{tumour}.vep.log"
  params:
    cores=cluster["annotate_vep_somatic_indels"]["n"]
  shell:
    "{config[module_samtools]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

rule annotate_vep_hc:
  input:
    vcf="out/germline_joint.hc.normalized.vcf",
    reference=config['genome']
  output:
    "out/germline_joint.hc.normalized.vep.vcf.gz"
  log:
    "log/hc.vep.log"
  params:
    cores=cluster["annotate_vep_germline"]["n"]
  shell:
    "{config[module_samtools]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

rule annotate_vep_hc_tumours:
  input:
    vcf="out/tumour_joint.hc.normalized.vcf",
    reference=config['genome']
  output:
    "out/tumour_joint.hc.normalized.vep.vcf.gz"
  log:
    "log/hc.vep.log"
  params:
    cores=cluster["annotate_vep_germline"]["n"]
  shell:
    "{config[module_samtools]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

rule annotate_vep_mutect2:
  input:
    vcf="out/{tumour}.mutect2.filter.vcf.gz",
    reference=config['genome']
  output:
    "out/{tumour}.mutect2.filter.vep.vcf.gz"
  log:
    "log/{tumour}.mutect2.vep.log"
  params:
    cores=cluster["annotate_vep_mutect2"]["n"]
  shell:
    "{config[module_samtools]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

rule annotate_vep_germline:
  input:
    vcf="out/{germline}.strelka.germline.filter_gt.vcf.gz",
    reference=config['genome']
  output:
    "out/{germline}.strelka.germline.filter_gt.vep.vcf.gz"
  log:
    "log/{germline}.strelka.vep.log"
  params:
    cores=cluster["annotate_vep_germline"]["n"]
  shell:
    "{config[module_samtools]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

### other post variant calling
# TODO https://github.com/hammerlab/concordance
#rule qc_vcf_concordance:
#  input:
#  output:
#  log:
#  shell:

rule intersect_somatic_callers:
  input:
    reference=config["genome"],
    mutect2="out/{tumour}.mutect2.filter.vep.vcf.gz",
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.vcf.gz" # TODO not used
  output:
    "out/{tumour}.intersect.vcf.gz"
  log:
    stderr="log/{tumour}.intersect.log"
  shell:
    "({config[module_samtools]} && "
    "{config[module_htslib]} && "
    "tools/vt-0.577/vt decompose -s {input.mutect2} | tools/vt-0.577/vt normalize -n -r {input.reference} - -o out/{wildcards.tumour}.mutect2.norm.vcf.gz && "
    "tools/vt-0.577/vt decompose -s {input.strelka_snvs} | tools/vt-0.577/vt normalize -n -r {input.reference} - -o out/{wildcards.tumour}.strelka.somatic.snvs.af.norm.vcf.gz && "
    "src/vcf_intersect.py --inputs out/{wildcards.tumour}.strelka.somatic.snvs.af.norm.vcf.gz out/{wildcards.tumour}.mutect2.norm.vcf.gz | bgzip > {output}) 2>{log.stderr}"

rule filter_intersected_somatic_callers:
  input:
    "out/{tumour}.intersect.vcf.gz"
  output:
    "out/{tumour}.intersect.filter.vcf.gz"
  log:
    stderr="log/{tumour}.intersect_filter.log"
  params:
    af=config["af_threshold"],
    dp=config["dp_threshold"],
    tumour="{tumour}"
  shell:
    "{config[module_htslib]} && "
    "src/filter_af.py --sample {params.tumour} --af {params.af} --dp {params.dp} < {input} 2>{log.stderr} | bgzip > {output}"

#rule annotate_af_germline:
#  input:
#    "out/{germline}.strelka.germline.vcf.gz",
#  output:
#    "out/{germline}.strelka.germline.af.vcf.gz",
#  log:
#    stderr="log/{germline}.annotate_af.stderr"
#  shell:
#    "{config[module_samtools]} && "
#    "src/annotate_af.py {input} | bgzip >{output} 2>{log.stderr}"

rule bias_filter_strelka:
  input:
    reference=config["genome"],
    bam="out/{tumour}.sorted.dups.bam",
    vcf="out/{tumour}.strelka.somatic.snvs.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.snvs.bias.vcf.gz"
  log:
    stderr="log/{tumour}.bias_filter.snvs.bias.err",
    stdout="log/{tumour}.bias_filter.snvs.bias.out"
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input.vcf} | egrep '(^#|PASS)' > tmp/{wildcards.tumour}_bias_filter_strelka.vcf && "
    "python tools/DKFZBiasFilter/scripts/biasFilter.py --tempFolder tmp tmp/{wildcards.tumour}_bias_filter_strelka.vcf {input.bam} {input.reference} tmp/{wildcards.tumour}_bias_filter_out_strelka.vcf 2>{log.stderr} 1>{log.stdout} && "
    "bgzip < tmp/{wildcards.tumour}_bias_filter_out_strelka.vcf > {output} && "
    "rm tmp/{wildcards.tumour}_bias_filter_strelka.vcf tmp/{wildcards.tumour}_bias_filter_out_strelka.vcf"

rule bias_filter_mutect2:
  input:
    reference=config["genome"],
    bam="out/{tumour}.sorted.dups.bam",
    vcf="out/{tumour}.mutect2.filter.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.bias.vcf.gz"
  log:
    stderr="log/{tumour}.bias_filter.mutect2.bias.err",
    stdout="log/{tumour}.bias_filter.mutect2.bias.out"
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input.vcf} | egrep '(^#|PASS)' > tmp/{wildcards.tumour}_bias_filter_mutect2.vcf && "
    "python tools/DKFZBiasFilter/scripts/biasFilter.py --tempFolder tmp tmp/{wildcards.tumour}_bias_filter_mutect2.vcf {input.bam} {input.reference} tmp/{wildcards.tumour}_bias_filter_out_mutect2.vcf 2>{log.stderr} 1>{log.stdout} && "
    "bgzip < tmp/{wildcards.tumour}_bias_filter_out_mutect2.vcf > {output} && "
    "rm tmp/{wildcards.tumour}_bias_filter_mutect2.vcf tmp/{wildcards.tumour}_bias_filter_out_mutect2.vcf"

# filter on genes of interest and convert to tsv
rule filter_genes_of_interest_tumour:
  input:
    vcf="out/{tumour}.mutect2.filter.vep.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.genes_of_interest.tsv"
  log:
    stderr="log/{tumour}.filter_genes_of_interest_tumour.err"
  params:
    gene_list=' '.join(config["genes_of_interest"])
  shell:
    "src/vcf2tsv.py {input.vcf} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK' | "
    "src/filter_tsv.py --column vep_SYMBOL --values {params.gene_list} > {output}"

# TODO
#rule filter_genes_of_interest_germline:
#  input:
#    vcf="out/{tumour}.mutect2.filter.vep.vcf.gz"
#  output:
#    "out/{tumour}.mutect2.filter.genes_of_interest.tsv"
#  log:
#    stderr="log/{tumour}.filter_genes_of_interest_tumour.err"
#  params:
#    gene_list=' '.join(config["genes_of_interest"])
#  shell:
#    "src/vcf2tsv.py {input.vcf} | "
#    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK' | "
#    "src/filter_tsv.py --column vep_SYMBOL --values {params.gene_list} > {output}"

# mutational signatures
rule mutational_signature:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.vcf.gz"
  output:
    "out/{tumour}.mutational_signature.exposures"
  log:
    stderr="out/{tumour}.mutational_signature.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-0.2/count.py --genome {input.reference} --vcf {input.vcf} > out/{wildcards.tumour}.mutational_signature.counts && "
    "python tools/mutational_signature-0.2/plot_counts.py out/{wildcards.tumour}.mutational_signature.png < out/{wildcards.tumour}.mutational_signature.counts && "
    "python tools/mutational_signature-0.2/decompose.py --signatures tools/mutational_signature-0.2/signatures30.txt --counts out/{wildcards.tumour}.mutational_signature.counts > {output}) 2>{log.stderr}"

rule combine_mutational_signatures:
  input:
    expand("out/{tumour}.mutational_signature.exposures", tumour=config['tumours']),
  output:
    "out/mutational_signatures.combined"
  shell:
    "src/combine_tsv.py {input} | sed 's/^out\\/\\(.*\)\\.mutational_signature\\.exposures/\\1/' > {output}"

# mutational signatures with filtered counts
rule mutational_signature_filtered:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.filter.vcf.gz"
  output:
    "out/{tumour}.mutational_signature.filter.exposures"
  log:
    stderr="out/{tumour}.mutational_signature.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-0.2/count.py --genome {input.reference} --vcf {input.vcf} > out/{wildcards.tumour}.mutational_signature.filter.counts && "
    "python tools/mutational_signature-0.2/plot_counts.py out/{wildcards.tumour}.mutational_signature.png < out/{wildcards.tumour}.mutational_signature.filter.counts && "
    "python tools/mutational_signature-0.2/decompose.py --signatures tools/mutational_signature-0.2/signatures30.txt --counts out/{wildcards.tumour}.mutational_signature.filter.counts > {output}) 2>{log.stderr}"

rule combine_mutational_signatures_filtered:
  input:
    expand("out/{tumour}.mutational_signature.filter.exposures", tumour=config['tumours']),
  output:
    "out/mutational_signatures.filter.combined"
  shell:
    "src/combine_tsv.py {input} | sed 's/^out\\/\\(.*\)\\.mutational_signature\\.filter\\.exposures/\\1/' >{output}"

# copy number
rule copy_number_varscan:
  input:
    bams=tumour_germline_bams,
    reference=config["genome"]
  output:
    "out/{tumour}.varscan.copynumber.called"
  params:
    tumour="{tumour}"
  shell:
    "{config[module_java]} && "
    "{config[module_samtools]} && "
    "samtools mpileup -q 1 -f {input.reference} {input.bams[1]} {input.bams[0]} > tmp/{params.tumour}.mpileups && "
    "java -jar tools/VarScan.v2.3.9.jar copynumber tmp/{params.tumour}.mpileups out/{params.tumour}.varscan --mpileup 1 && "
    "java -jar tools/VarScan.v2.3.9.jar copyCaller out/{params.tumour}.varscan.copynumber --output-file {output}"

rule copy_number_varscan_post:
  input:
    "out/{tumour}.varscan.copynumber.called"
  output:
    "out/{tumour}.varscan.copynumber.deletions.bed"
  params:
    tumour="{tumour}"
  shell:
    "{config[module_network]} && {config[module_R]} && "
    "sed '1d' < {input} > tmp/{params.tumour}.varscan.nohead && "
    "src/varscan_cnv_post.R --in tmp/{params.tumour}.varscan.nohead --out out/{params.tumour}.varscan.merged && " # 1       13360   16851   10      0.7773
    "awk -v OFS='\t' '{{ if ($5 < -0.1) {{len=$3-$2; print $1, $2, $3, \"logR=\" $5 \";length=\" len \";markers=\" $4}} }}' < out/{params.tumour}.varscan.merged > {output}"

# burden (currently only exonic snvs)
rule mutation_burden:
  input:
    vcfs=expand("out/{tumour}.intersect.filter.vcf.gz", tumour=config['tumours']),
    regions="reference/regions.bed"
  output:
    "out/mutation_rate.tsv"
  log:
    stderr="log/mutation_rate.stderr"
  shell:
    "src/mutation_rate.py --verbose --vcfs {input.vcfs} --bed {input.regions} >{output} 2>{log.stderr}"

# loh
rule loh:
  input:
    snvs="out/{tumour}.strelka.somatic.snvs.af.vep.vcf.gz", # loh requires strelka for now
    indels="out/{tumour}.strelka.somatic.indels.vep.vcf.gz"
  output:
    "out/{tumour}.loh.bed"
  log:
    stderr="log/{tumour}.loh.stderr"
  params:
    tumour="{tumour}",
    regions=' '.join(config["loh_regions"]),
    region_names=' '.join(config["loh_region_names"]),
    region_padding=' '.join(config["loh_region_padding"])
  shell:
    "(tools/loh_caller-0.1/loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 20 --neutral --min_af 0.1 < {input.snvs} > tmp/{params.tumour}.loh.snvs.tsv && "
    "tools/loh_caller-0.1/loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 20 --neutral --min_af 0.1 < {input.indels} > tmp/{params.tumour}.loh.indels.tsv && "
    "sort -k1,1 -k2,2n tmp/{params.tumour}.loh.snvs.tsv tmp/{params.tumour}.loh.indels.tsv > tmp/{params.tumour}.loh.tsv && "
    "tools/loh_caller-0.1/loh_merge.py --verbose --noheader --min_len 1000 --min_prop 0.1 --plot out/{params.tumour}.loh --regions {params.regions} --region_names {params.region_names} --region_padding {params.region_padding} --plot_chromosomes >{output}) 2>{log.stderr}"
