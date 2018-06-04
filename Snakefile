
configfile: "cfg/config.yaml"

cluster = json.load(open("cfg/cluster.json"))

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

def tumour_normal_bams(wildcards):
  tumour_bam = 'out/{}.sorted.bam'.format(wildcards.tumour)
  normal_bam = 'out/{}.sorted.bam'.format(config["tumours"][wildcards.tumour])
  return [tumour_bam, normal_bam]

### final outputs ###
rule all:
  input:
    expand("out/{tumour}.strelka.snvs.vcf.gz", tumour=config['tumours']),
    expand("out/{tumour}.strelka.indels.vcf.gz", tumour=config['tumours']),
    expand("out/{sample}.oxo_metrics.txt", sample=config['samples']),
    expand("out/{sample}.artifact_metrics.txt.error_summary_metrics", sample=config['samples']),
    expand("out/{tumour}.strelka.snvs.bias.vcf", tumour=config['tumours']),
    expand("out/fastqc/{sample}/completed", sample=config['samples']),
    expand("out/{sample}.metrics.insertsize", sample=config['samples']),
    expand("out/{sample}.metrics.alignment", sample=config['samples']),
    "out/qc.summary.tsv",
    "out/multiqc.html"

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

rule qc_insertsize:
  input:
    reference=config["genome"],
    bam="out/{sample}.sorted.bam"
  output:
    "out/{sample}.metrics.insertsize"
  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE={input.reference} INPUT={input.bam} OUTPUT={output}"

rule qc_alignment:
  input:
    bam="out/{sample}.sorted.bam"
  output:
    "out/{sample}.metrics.alignment"
  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar CollectInsertSizeMetrics INPUT={input.bam} OUTPUT={output} HISTOGRAM_FILE={output}.pdf"

rule multiqc:
  input:
    expand("out/fastqc/{sample}/completed", sample=config['samples']),
    expand("out/{sample}.metrics.alignment", sample=config['samples']),
    expand("out/{sample}.metrics.insertsize", sample=config['samples']),
    "out/qc.summary.tsv"
    
  output:
    "out/multiqc.html"
  shell:
    "multiqc --force --filename {output} out"

### alignment ###
rule align:
  input:
    reference=config["genome"],
    fastqs=lambda wildcards: config["samples"][wildcards.sample]

  output:
    "out/{sample}.bam"

  log:
    "log/{sample}.bwa.log"

  params:
    cores=cluster["align"]["n"],
    read_group=read_group

  shell:
    "module load bwa-intel/0.7.12 && module load samtools-intel/1.4 && "
    "(bwa mem -M -t {params.cores} -R \"{params.read_group}\" {input.reference} {input.fastqs} | samtools view -b -h -o {output} -) 2>{log}"

# sort the bam
rule sort:
  input:
    "out/{sample}.bam"

  output:
    "out/{sample}.sorted.bam"

  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar SortSam INPUT={input} OUTPUT={output} VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=2000000 CREATE_INDEX=True"

### qc ###
rule qc_sequencing_artifacts:
  input:
    bam="out/{sample}.sorted.bam",
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
    bam="out/{sample}.sorted.bam",
    reference=config["genome"]

  output:
    "out/{sample}.oxo_metrics.txt"

  log:
    stderr="log/{sample}.oxo.err",
    stdout="log/{sample}.oxo.out"

  shell:
    "module load java/1.8.0_25 && "
    "java -jar tools/picard-2.8.2.jar CollectOxoGMetrics I={input.bam} O={output} R={input.reference} 2>{log.stderr} 1>{log.stdout}"

### variant calling ###
rule strelka:
  input:
    reference=config["genome"],
    bams=tumour_normal_bams

  output:
    "out/{tumour}.strelka.snvs.vcf.gz",
    "out/{tumour}.strelka.indels.vcf.gz",

  log:
    "log/{tumour}.strelka.log"

  params:
    cores=cluster["strelka"]["n"]

  shell:
    "(mkdir -p tmp/strelka_{wildcards.tumour} && "
    "rm -rf tmp/strelka_{wildcards.tumour}/* && "
    "tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py "
    "--ref {input.reference} "
    "--tumorBam {input.bams[0]} "
    "--normalBam {input.bams[1]} "
    "--runDir tmp/strelka_{wildcards.tumour} "
    "--exome && "
    "tmp/strelka_{wildcards.tumour}/runWorkflow.py -m local -j {params.cores} && "
    "mv tmp/strelka_{wildcards.tumour}/results/variants/somatic.snvs.vcf.gz {output[0]} && " 
    "mv tmp/strelka_{wildcards.tumour}/results/variants/somatic.indels.vcf.gz {output[1]} && " 
    "rm -r tmp/strelka_{wildcards.tumour}) 2>{log}"

rule bias_filter:
  input:
    reference=config["genome"],
    bam="out/{tumour}.sorted.bam",
    vcf="out/{tumour}.strelka.snvs.vcf.gz"
  output:
    "out/{tumour}.strelka.snvs.bias.vcf"
  log:
    stderr="log/{tumour}.bias_filter.snvs.bias.err",
    stdout="log/{tumour}.bias_filter.snvs.bias.out"
  shell:
    "gunzip < {input.vcf} > tmp/bias_filter_$$.vcf && python tools/DKFZBiasFilter/scripts/biasFilter.py --tempFolder tmp tmp/bias_filter_$$.vcf {input.bam} {input.reference} {output} && rm tmp/bias_filter_$$.vcf"
