
configfile: "cfg/config.yaml"

cluster = json.load(open("cfg/cluster.json"))

### helper functions ###
def read_group(wildcards):
  '''
    determine read group from sample name
  '''
  fields = config["samples"][wildcards.sample][0].split("_") # in/S1_RG_2_R1.fastq.gz
  readid = fields[1]
  lane = fields[2]
  return "@RG\tID:{readid}\tSM:{sample}_{readid}\tPU:lib1\tLN:{lane}\tPL:Illumina".format(readid=readid, sample=wildcards.sample, lane=lane)

def tumour_normal_bams(wildcards):
  tumour_bam = 'out/{}.sorted.bam'.format(wildcards.tumour)
  normal_bam = 'out/{}.sorted.bam'.format(config["tumours"][wildcards.tumour])
  return [tumour_bam, normal_bam]

### final outputs ###
rule all:
  input:
    expand("out/{tumour}.strelka.snvs.vcf.gz", tumour=config['tumours']),
    expand("out/{tumour}.strelka.indels.vcf.gz", tumour=config['tumours'])

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
    "java -jar tools/picard-2.8.2.jar SortSam INPUT={input} OUTPUT={output} VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True"

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

