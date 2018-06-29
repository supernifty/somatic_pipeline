
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

### final outputs ###
rule all:
  input:
    expand("out/{germline}.hc.gvcf.gz", germline=germline_samples()),
    expand("out/{tumour}.strelka.somatic.snvs.af.vep.vcf.gz", tumour=config['tumours']),
    expand("out/{tumour}.strelka.somatic.indels.vep.vcf.gz", tumour=config['tumours']),
    expand("out/{germline}.strelka.germline.filter_gt.vep.vcf.gz", germline=germline_samples()),
    expand("out/{sample}.oxo_metrics.txt", sample=config['samples']),
    expand("out/{sample}.artifact_metrics.txt.error_summary_metrics", sample=config['samples']),
    #expand("out/{tumour}.strelka.somatic.snvs.bias.vcf", tumour=config['tumours']),
    expand("out/fastqc/{sample}/completed", sample=config['samples']),
    expand("out/{sample}.metrics.insertsize", sample=config['samples']),
    expand("out/{sample}.metrics.alignment", sample=config['samples']),
    expand("out/{sample}.metrics.target", sample=config['samples']),
    expand("out/{tumour}.concordance", tumour=config['tumours']),
    expand("out/{tumour}.contamination", tumour=config['tumours']),
    expand("out/{tumour}.verifybamid.somatic.completed", tumour=config['tumours']),
    expand("out/{tumour}.mutect2.filter.vep.vcf.gz", tumour=config['tumours']),
    "out/germline_joint.hc.normalized.vep.vcf.gz",
    "out/qc.summary.tsv",
    "out/multiqc.html"

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
    "module load java/1.8.0_25 && "
    "mkdir -p {output} && "
    "tools/FastQC/fastqc --extract --outdir out/fastqc/{wildcards.sample} {input.fastqs} && "
    "touch {output}"

rule make_sequence_dict:
  input:
    reference=config["genome"]
  output:
    config["genome_dict"]
  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar CreateSequenceDictionary REFERENCE={input.reference} OUTPUT={output}"

rule make_intervals:
  input:
    bed=config["regions"],
    dict=config["genome_dict"]
  output:
    "out/regions.intervals"
  shell:  
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar BedToIntervalList INPUT={input.bed} OUTPUT={output} SEQUENCE_DICTIONARY={input.dict}"

rule qc_target:
  input:
    reference=config["genome"],
    bam="out/{sample}.sorted.dups.bam",
    intervals="out/regions.intervals"
  output:
    "out/{sample}.metrics.target"
  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar CollectHsMetrics REFERENCE_SEQUENCE={input.reference} INPUT={input.bam} OUTPUT={output} BAIT_INTERVALS={input.intervals} TARGET_INTERVALS={input.intervals}"

rule qc_insertsize:
  input:
    reference=config["genome"],
    bam="out/{sample}.sorted.dups.bam"
  output:
    "out/{sample}.metrics.insertsize"
  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE={input.reference} INPUT={input.bam} OUTPUT={output}"

rule qc_alignment:
  input:
    bam="out/{sample}.sorted.dups.bam"
  output:
    "out/{sample}.metrics.alignment"
  shell:
    "module load java/1.8.0_25 && "
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
    "module load java/1.8.0_25 && "
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


rule multiqc:
  input:
    expand("out/fastqc/{sample}/completed", sample=config['samples']),
    expand("out/{sample}.metrics.alignment", sample=config['samples']),
    expand("out/{sample}.metrics.insertsize", sample=config['samples']),
    expand("out/{sample}.metrics.target", sample=config['samples']),
    expand("out/{tumour}.concordance", tumour=config['tumours']),
    expand("out/{tumour}.contamination", tumour=config['tumours']),
    expand("out/{tumour}.verifybamid.somatic.completed", tumour=config['tumours']),
    "out/qc.summary.tsv"
    
  output:
    "out/multiqc.html"
  shell:
    "multiqc --force --filename {output} out"

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
    stderr="log/{sample}.trimmomatic.stderr"
  params:
    cores=cluster["align"]["n"],
  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads {params.cores} -phred33 {input.fastqs} {output} ILLUMINACLIP:tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>{log.stderr}"

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
    "module load bwa-intel/0.7.12 && module load samtools-intel/1.4 && "
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
    "module load bwa-intel/0.7.12 && module load samtools-intel/1.4 && "
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
    "module load samtools-intel/1.4 && "
    "samtools merge {output} {input} 2>{log}"

# sort the bam
rule sort:
  input:
    "tmp/{sample}.merged.bam"

  output:
    bam="tmp/{sample}.sorted.bam",
    bai="tmp/{sample}.sorted.bai"

  shell:
    "module load java/1.8.0_25 && "
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
    "module load java/1.8.0_25 && "
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
    "(module load java/1.8.0_25 && "
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
    "(module load java/1.8.0_25 && "
    "java -jar tools/GenomeAnalysisTK-3.7.0.jar -T CombineGVCFs -R {input.reference} {params.variant_list} -L {wildcards.chromosome} -o tmp/germline_combined_{wildcards.chromosome}.gvcf && "
    "tools/gatk-4.0.0.0/gatk GenotypeGVCFs -R {input.reference} --dbsnp reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz -V tmp/germline_combined_{wildcards.chromosome}.gvcf -L {wildcards.chromosome} --use-new-qual-calculator true --output out/germline_joint_{wildcards.chromosome}.vcf"
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
    "(module load java/1.8.0_25 && module load R-gcc/3.4.0 && module load samtools-intel/1.8 && "
    "tools/gatk-4.0.0.0/gatk GatherVcfs -R {input.reference} --OUTPUT=tmp/germline_joint.vcf {params.inputs} && "
    "bgzip -c < tmp/germline_joint.vcf > tmp/germline_joint.vcf.bgz && tabix -p vcf tmp/germline_joint.vcf.bgz && "
    "tools/gatk-4.0.0.0/gatk CalculateGenotypePosteriors -R {input.reference} --supporting reference/gatk-4-bundle-b37/1000G_phase3_v4_20130502.sites.vcf.bgz -V tmp/germline_joint.vcf.bgz -O tmp/germline_joint.cgp.vcf && "
    "tools/vt-0.577/vt normalize -n -r {input.reference} tmp/germline_joint.cgp.vcf -o tmp/germline_joint.cgp.normalized.vcf && "
    "tools/vt-0.577/vt decompose -s tmp/germline_joint.cgp.normalized.vcf | tools/vt-0.577/vt normalize -r {input.reference} - -o {output}"
    ") 2>{log}"

### qc ###
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
    "module load java/1.8.0_25 && "
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
    "module load java/1.8.0_25 && "
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
    "module load samtools-intel/1.5 && "
    "src/annotate_af.py {input} | bgzip >{output} 2>{log.stderr}"

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
    "tools/gatk-4.0.0.0/gatk Mutect2 -R {input.reference} -I {input.bams[0]} -I {input.bams[1]} --tumor-sample {wildcards.tumour} --normal-sample {params.germline} --output {output} --output-mode EMIT_VARIANTS_ONLY --dbsnp {input.dbsnp} --germline-resource {input.gnomad} --af-of-alleles-not-in-resource 0.0000025 -pon {input.pon} --interval-padding 1000 -L {input.regions} -L {wildcards.chromosome} --interval-set-rule INTERSECTION --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter"

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
    "module load java/1.8.0_25 && "
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
    "(tools/gatk-4.0.0.0/gatk GetPileupSummaries -I {input.bam} -V {input.gnomad} -O tmp/{wildcards.tumour}.mutect2.pileup.table && "
    "tools/gatk-4.0.0.0/gatk CalculateContamination -I tmp/{wildcards.tumour}.mutect2.pileup.table -O tmp/{wildcards.tumour}.mutect2.contamination.table && "
    "tools/gatk-4.0.0.0/gatk FilterMutectCalls -V {input.vcf} --contamination-table tmp/{wildcards.tumour}.mutect2.contamination.table -O {output}) 1>{log.stdout} 2>{log.stderr}"

### germline variant calling ###
rule filter_strelka_germline:
  input:
    "out/{germline}.strelka.germline.vcf.gz",
  output:
    "out/{germline}.strelka.germline.filter_gt.vcf.gz",
  shell:
    "module load samtools-intel/1.5 && "
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
    "module load samtools-intel/1.5 && "
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
    "module load samtools-intel/1.5 && "
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
    "module load samtools-intel/1.5 && "
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
    "module load samtools-intel/1.5 && "
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
    "module load samtools-intel/1.5 && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

#rule annotate_af_germline:
#  input:
#    "out/{germline}.strelka.germline.vcf.gz",
#  output:
#    "out/{germline}.strelka.germline.af.vcf.gz",
#  log:
#    stderr="log/{germline}.annotate_af.stderr"
#  shell:
#    "module load samtools-intel/1.5 && "
#    "src/annotate_af.py {input} | bgzip >{output} 2>{log.stderr}"

#rule bias_filter:
#  input:
#    reference=config["genome"],
#    bam="out/{tumour}.sorted.dups.bam",
#    vcf="out/{tumour}.strelka.somatic.snvs.vcf.gz"
#  output:
#    "out/{tumour}.strelka.somatic.snvs.bias.vcf"
#  log:
#    stderr="log/{tumour}.bias_filter.snvs.bias.err",
#    stdout="log/{tumour}.bias_filter.snvs.bias.out"
#  shell:
#    "gunzip < {input.vcf} > tmp/bias_filter_$$.vcf && python tools/DKFZBiasFilter/scripts/biasFilter.py --tempFolder tmp tmp/bias_filter_$$.vcf {input.bam} {input.reference} {output} && rm tmp/bias_filter_$$.vcf"
