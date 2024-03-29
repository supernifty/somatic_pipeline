print("starting")
VERSION="0.3"

configfile: "cfg/config.yaml"

import yaml
cluster = yaml.load(open("cfg/cluster.yaml"))
samples = yaml.load(open("cfg/samples.yaml"))
print("loaded configuration")


# NOTE: no mitochondria MT because they aren't in our exome
GATK_CHROMOSOMES=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')

# split into many more segments to speed up WGS
#GATK_CHROMOSOMES=('1:1-10000100', '1:9999901-20000100', '1:19999901-30000100', '1:29999901-40000100', '1:39999901-50000100', '1:49999901-60000100', '1:59999901-70000100', '1:69999901-80000100', '1:79999901-90000100', '1:89999901-100000100', '1:99999901-110000100', '1:109999901-120000100', '1:119999901-130000100', '1:129999901-140000100', '1:139999901-150000100', '1:149999901-160000100', '1:159999901-170000100', '1:169999901-180000100', '1:179999901-190000100', '1:189999901-200000100', '1:199999901-210000100', '1:209999901-220000100', '1:219999901-230000100', '1:229999901-240000100', '1:239999901-249250621', '2:1-10000100', '2:9999901-20000100', '2:19999901-30000100', '2:29999901-40000100', '2:39999901-50000100', '2:49999901-60000100', '2:59999901-70000100', '2:69999901-80000100', '2:79999901-90000100', '2:89999901-100000100', '2:99999901-110000100', '2:109999901-120000100', '2:119999901-130000100', '2:129999901-140000100', '2:139999901-150000100', '2:149999901-160000100', '2:159999901-170000100', '2:169999901-180000100', '2:179999901-190000100', '2:189999901-200000100', '2:199999901-210000100', '2:209999901-220000100', '2:219999901-230000100', '2:229999901-240000100', '2:239999901-243199373', '3:1-10000100', '3:9999901-20000100', '3:19999901-30000100', '3:29999901-40000100', '3:39999901-50000100', '3:49999901-60000100', '3:59999901-70000100', '3:69999901-80000100', '3:79999901-90000100', '3:89999901-100000100', '3:99999901-110000100', '3:109999901-120000100', '3:119999901-130000100', '3:129999901-140000100', '3:139999901-150000100', '3:149999901-160000100', '3:159999901-170000100', '3:169999901-180000100', '3:179999901-190000100', '3:189999901-198022430', '4:1-10000100', '4:9999901-20000100', '4:19999901-30000100', '4:29999901-40000100', '4:39999901-50000100', '4:49999901-60000100', '4:59999901-70000100', '4:69999901-80000100', '4:79999901-90000100', '4:89999901-100000100', '4:99999901-110000100', '4:109999901-120000100', '4:119999901-130000100', '4:129999901-140000100', '4:139999901-150000100', '4:149999901-160000100', '4:159999901-170000100', '4:169999901-180000100', '4:179999901-190000100', '4:189999901-191154276', '5:1-10000100', '5:9999901-20000100', '5:19999901-30000100', '5:29999901-40000100', '5:39999901-50000100', '5:49999901-60000100', '5:59999901-70000100', '5:69999901-80000100', '5:79999901-90000100', '5:89999901-100000100', '5:99999901-110000100', '5:109999901-120000100', '5:119999901-130000100', '5:129999901-140000100', '5:139999901-150000100', '5:149999901-160000100', '5:159999901-170000100', '5:169999901-180000100', '5:179999901-180915260', '6:1-10000100', '6:9999901-20000100', '6:19999901-30000100', '6:29999901-40000100', '6:39999901-50000100', '6:49999901-60000100', '6:59999901-70000100', '6:69999901-80000100', '6:79999901-90000100', '6:89999901-100000100', '6:99999901-110000100', '6:109999901-120000100', '6:119999901-130000100', '6:129999901-140000100', '6:139999901-150000100', '6:149999901-160000100', '6:159999901-170000100', '6:169999901-171115067', '7:1-10000100', '7:9999901-20000100', '7:19999901-30000100', '7:29999901-40000100', '7:39999901-50000100', '7:49999901-60000100', '7:59999901-70000100', '7:69999901-80000100', '7:79999901-90000100', '7:89999901-100000100', '7:99999901-110000100', '7:109999901-120000100', '7:119999901-130000100', '7:129999901-140000100', '7:139999901-150000100', '7:149999901-159138663', 'X:1-10000100', 'X:9999901-20000100', 'X:19999901-30000100', 'X:29999901-40000100', 'X:39999901-50000100', 'X:49999901-60000100', 'X:59999901-70000100', 'X:69999901-80000100', 'X:79999901-90000100', 'X:89999901-100000100', 'X:99999901-110000100', 'X:109999901-120000100', 'X:119999901-130000100', 'X:129999901-140000100', 'X:139999901-150000100', 'X:149999901-155270560', '8:1-10000100', '8:9999901-20000100', '8:19999901-30000100', '8:29999901-40000100', '8:39999901-50000100', '8:49999901-60000100', '8:59999901-70000100', '8:69999901-80000100', '8:79999901-90000100', '8:89999901-100000100', '8:99999901-110000100', '8:109999901-120000100', '8:119999901-130000100', '8:129999901-140000100', '8:139999901-146364022', '9:1-10000100', '9:9999901-20000100', '9:19999901-30000100', '9:29999901-40000100', '9:39999901-50000100', '9:49999901-60000100', '9:59999901-70000100', '9:69999901-80000100', '9:79999901-90000100', '9:89999901-100000100', '9:99999901-110000100', '9:109999901-120000100', '9:119999901-130000100', '9:129999901-140000100', '9:139999901-141213431', '10:1-10000100', '10:9999901-20000100', '10:19999901-30000100', '10:29999901-40000100', '10:39999901-50000100', '10:49999901-60000100', '10:59999901-70000100', '10:69999901-80000100', '10:79999901-90000100', '10:89999901-100000100', '10:99999901-110000100', '10:109999901-120000100', '10:119999901-130000100', '10:129999901-135534747', '11:1-10000100', '11:9999901-20000100', '11:19999901-30000100', '11:29999901-40000100', '11:39999901-50000100', '11:49999901-60000100', '11:59999901-70000100', '11:69999901-80000100', '11:79999901-90000100', '11:89999901-100000100', '11:99999901-110000100', '11:109999901-120000100', '11:119999901-130000100', '11:129999901-135006516', '12:1-10000100', '12:9999901-20000100', '12:19999901-30000100', '12:29999901-40000100', '12:39999901-50000100', '12:49999901-60000100', '12:59999901-70000100', '12:69999901-80000100', '12:79999901-90000100', '12:89999901-100000100', '12:99999901-110000100', '12:109999901-120000100', '12:119999901-130000100', '12:129999901-133851895', '13:1-10000100', '13:9999901-20000100', '13:19999901-30000100', '13:29999901-40000100', '13:39999901-50000100', '13:49999901-60000100', '13:59999901-70000100', '13:69999901-80000100', '13:79999901-90000100', '13:89999901-100000100', '13:99999901-110000100', '13:109999901-115169878', '14:1-10000100', '14:9999901-20000100', '14:19999901-30000100', '14:29999901-40000100', '14:39999901-50000100', '14:49999901-60000100', '14:59999901-70000100', '14:69999901-80000100', '14:79999901-90000100', '14:89999901-100000100', '14:99999901-107349540', '15:1-10000100', '15:9999901-20000100', '15:19999901-30000100', '15:29999901-40000100', '15:39999901-50000100', '15:49999901-60000100', '15:59999901-70000100', '15:69999901-80000100', '15:79999901-90000100', '15:89999901-100000100', '15:99999901-102531392', '16:1-10000100', '16:9999901-20000100', '16:19999901-30000100', '16:29999901-40000100', '16:39999901-50000100', '16:49999901-60000100', '16:59999901-70000100', '16:69999901-80000100', '16:79999901-90000100', '16:89999901-90354753', '17:1-10000100', '17:9999901-20000100', '17:19999901-30000100', '17:29999901-40000100', '17:39999901-50000100', '17:49999901-60000100', '17:59999901-70000100', '17:69999901-80000100', '17:79999901-81195210', '18:1-10000100', '18:9999901-20000100', '18:19999901-30000100', '18:29999901-40000100', '18:39999901-50000100', '18:49999901-60000100', '18:59999901-70000100', '18:69999901-78077248', '20:1-10000100', '20:9999901-20000100', '20:19999901-30000100', '20:29999901-40000100', '20:39999901-50000100', '20:49999901-60000100', '20:59999901-63025520', 'Y:1-10000100', 'Y:9999901-20000100', 'Y:19999901-30000100', 'Y:29999901-40000100', 'Y:39999901-50000100', 'Y:49999901-59373566', '19:1-10000100', '19:9999901-20000100', '19:19999901-30000100', '19:29999901-40000100', '19:39999901-50000100', '19:49999901-59128983', '22:1-10000100', '22:9999901-20000100', '22:19999901-30000100', '22:29999901-40000100', '22:39999901-50000100', '22:49999901-51304566', '21:1-10000100', '21:9999901-20000100', '21:19999901-30000100', '21:29999901-40000100', '21:39999901-48129895')

### helper functions ###
def read_group(wildcards):
  '''
    determine read group from sample name
  '''
  # in/0656045001_BC_H5CNYDSXX_GTGTTCTA_L004_R1.fastq.gz
  # 0757079003_T_A_H5CNYDSXX_ACGTATCA_L004_R1.fastq.gz
  suffix = samples["samples"][wildcards.sample][0].replace("in/{}_".format(wildcards.sample), "") # in/S1_RG_2_R1.fastq.gz
  fields = suffix.split("_") # H5CNYDSXX_GTGTTCTA_L004_R1 
  flowcell = fields[0]
  barcode = fields[1]
  lane = fields[2]
  return "@RG\tID:{sample}.{flowcell}.{barcode}.{lane}\tSM:{sample}\tPU:{flowcell}.{barcode}.{lane}\tPL:Illumina".format(flowcell=flowcell, sample=wildcards.sample, lane=lane, barcode=barcode)

def tumour_germline_dup_bams(wildcards):
  tumour_bam = 'out/{}.sorted.dups.bam'.format(wildcards.tumour)
  normal_bam = 'out/{}.sorted.dups.bam'.format(samples["tumours"][wildcards.tumour])
  return [tumour_bam, normal_bam]

def tumour_germline_bams(wildcards):
  tumour_bam = 'out/{}.sorted.dups.bam'.format(wildcards.tumour)
  normal_bam = 'out/{}.sorted.dups.bam'.format(samples["tumours"][wildcards.tumour])
  return [tumour_bam, normal_bam]

def tumour_germline_bais(wildcards):
  tumour_bai = 'out/{}.sorted.dups.bai'.format(wildcards.tumour)
  normal_bai = 'out/{}.sorted.dups.bai'.format(samples["tumours"][wildcards.tumour])
  return [tumour_bai, normal_bai]

def germline_samples():
  all_samples = set(samples['samples'])
  tumours = set(samples['tumours'])
  return list(all_samples.difference(tumours))

def tumour_samples():
  return samples['tumours']

def germline_sample(wildcards):
  return samples["tumours"][wildcards.tumour]

def capture_specific(capture, l):
  if samples["capture"] == capture:
    return l
  else:
    return []

### final outputs ###
rule all:
  input:
    expand("out/{tumour}.strelka.somatic.snvs.af.norm.annot.pass.vcf.gz", tumour=samples['tumours']), # somatic snvs strelka
    expand("out/{tumour}.strelka.somatic.indels.norm.annot.pass.vcf.gz", tumour=samples['tumours']), # somatic indels strelka
    expand("out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz", tumour=samples['tumours']), # somatic indels strelka
    expand("out/{germline}.strelka.germline.filter_gt.vep.vcf.gz", germline=germline_samples()), # germline strelka calls
    # -- disabled minimal benefit expand("out/{germline}.strelka.germline.filter_gt.vep.vcf.gz", germline=samples['tumours']), # run the tumours as if they were germline

    expand("out/{sample}.oxo_metrics.txt", sample=samples['samples']),
    expand("out/{sample}.artifact_metrics.txt.error_summary_metrics", sample=samples['samples']),
    # -- disabled too slow expand("out/{tumour}.strelka.somatic.snvs.bias.vcf.gz", tumour=samples['tumours']),
    # -- disabled expand("out/{tumour}.mutect2.filter.bias.vcf.gz", tumour=samples['tumours']), # somatic mutect2 with dkfz bias annotation

    expand("out/fastqc/{sample}/completed", sample=samples['samples']), # fastqc
    expand("out/{sample}.metrics.insertsize", sample=samples['samples']),
    expand("out/{sample}.metrics.alignment", sample=samples['samples']),
    expand("out/{sample}.metrics.target", sample=samples['samples']),

    # -- issues expand("out/{tumour}.verifybamid.somatic.completed", tumour=samples['tumours']),
    expand("out/{tumour}.intersect.pass.filter.capture.vcf.gz", tumour=samples['tumours']), # somatic combination of strelka and mutect2, with filter, intersect with capture
    expand("out/{tumour}.pass_one.vcf.gz", tumour=samples['tumours']), # somatic combination of strelka and mutect2, with filter

    # -- expand("out/{tumour}.varscan.copynumber.deletions.bed", tumour=samples['tumours']), # TODO fails on mini sample
# --    expand("out/{tumour}.platypus.somatic.vcf.gz", tumour=samples['tumours']), # platypus somatic calls
    expand("out/{tumour}.loh.bed", tumour=samples['tumours']), # loh regions
    expand("out/{tumour}.pass.loh.bed", tumour=samples['tumours']), # loh regions
    expand("out/{tumour}.kataegis.png", tumour=samples['tumours']), # kataegis regions
# --    expand("out/{tumour}.cnv.tsv", tumour=samples['tumours']), # cnv regions
# --    "out/cnvkit.reference.cnn", # cnvkit
    expand("out/{tumour}.strelka.somatic.all.af.png", tumour=samples['tumours']), # plot strelka af
    expand("out/{tumour}.strelka.somatic.pass.af.png", tumour=samples['tumours']), # plot strelka af
    expand("out/{tumour}.strelka.somatic.pass.signatures.af.png", tumour=samples['tumours']), # plot strelka af
    expand("out/{tumour}.mutect2.somatic.af.png", tumour=samples['tumours']), # plot mutect2 af
    expand("out/{tumour}.intersect.pass.filter.signatures.vcf.gz", tumour=samples['tumours']), # annotated outputs with signatures
# --    expand("out/{tumour}.strelka.somatic.mutational_signature_v3.2_sbs.png", tumour=samples['tumours']), # signature profile with signature allocation

    expand("out/{sample}.mini.bam", sample=samples['samples']), 

    # tumour purity
    # -- not working expand("out/{tumour}.theta2.done", tumour=samples['tumours']),
    "out/aggregate/purity.tsv",

    # msi
    "out/aggregate/msisensor.tsv",
    "out/aggregate/mantis.tsv",
    "out/aggregate/msiseq.tsv",


    # combined results
    "out/aggregate/concordance_matrix.tsv",
    "out/aggregate/loh.genes.tsv",
    "out/aggregate/mutect2.filter.genes_of_interest.combined.tsv",
    "out/aggregate/mutect2.filter.combined.tsv.gz",
    "out/aggregate/intersect.filter.combined.tsv",
    "out/aggregate/strelka.all.tsv.gz",

    "out/aggregate/mutational_signatures_v2.combined.tsv",
    "out/aggregate/mutational_signatures_v2.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_sbs.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_sbs_capture.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_dbs.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_dbs.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_id.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_id_strelka.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.pass.dp_af.combined.tsv",

    "out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.contexts.tsv", # context counts
    "out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.pass.contexts.tsv",

    "out/aggregate/mutational_signatures_v3.2_sbs.filter.pks.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_id.filter.pks.combined.tsv",
    "out/aggregate/mutational_signatures_v3.2_id_strelka.filter.pks.combined.tsv",
    "out/aggregate/extended_contexts.v2.tsv",
    "out/aggregate/mutational_signatures_v3.2_sbs_alexandrov.filter.pks.combined.tsv",

    "out/aggregate/mutational_signatures_v3.2_sbs.combined.tsv",

    "out/aggregate/max_coverage.tsv",
    #"out/aggregate/max_trimmed_coverage.tsv",
    "out/aggregate/ontarget.tsv",
    "out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.genes_of_interest.tsv.gz", # gatk calls for all germline samples (disabled for now kinda slow)
    "out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.cadd.tsv.gz",
    # -- disabled "out/aggregate/tumour_joint.hc.normalized.annot.vcf.gz", # gatk calls for all tumour samples (as germline)
    "out/aggregate/qc.summary.tsv",
    "out/aggregate/multiqc.html", # overall general qc
    "out/aggregate/ontarget.png", # combined ontarget coverage plots
    "out/aggregate/ontarget_tumour.png", # somatic ontarget coverage plots
    "out/aggregate/ontarget_germline.png", # germline ontarget coverage plots
    "out/aggregate/ontarget_tumour.dedup.png",
    "out/aggregate/mutation_rate.tsv",
# --    "out/aggregate/mutation_rate.artefact_filter.tsv",
    "out/aggregate/mutation_rate.high_af.tsv",
    "out/aggregate/mutation_rate.low_af.tsv",
    "out/aggregate/msi_burden.tsv",
    # -- disabled "out/aggregate/targetted_gene_summary.somatic.tsv",
    # -- disabled "out/aggregate/targetted_gene_summary.germline.tsv",
    "out/aggregate/gene_coverage.exons.tsv",
    "out/aggregate/gene_coverage_min.tumour.exons.tsv",
    # flagged exons in regions of interest
    "out/aggregate/gene_coverage_min.tumour.exons.qc.flagged.tsv",
    "out/aggregate/gene_coverage_min.germline.exons.qc.flagged.tsv",
    capture_specific("panel", ["out/aggregate/gene_coverage_min.tumour.regions.tsv", "out/aggregate/gene_coverage_min.germline.regions.tsv"]),
# --    "out/aggregate/final.html",

    "out/aggregate/somalier.done",
    "out/aggregate/msmutect.combined.tsv"

# -- disabled for wgs    "out/aggregate/optitype.tsv",

### aggregate ###

# write out all tool versions (TODO)
rule make_versions:
  output:
    versions="out/aggregate/versions.txt"
  shell:
    "src/make_tsv.py --columns Tool Version --rows "
    "Pipeline,{VERSION} "
    ">{output.versions}"

# -- all_germline_variants="out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.tsv.gz"
rule report_md:
  input:
    versions="out/aggregate/versions.txt",
    signatures="out/aggregate/mutational_signatures_v3.2_sbs.filter.combined.tsv",
    burden="out/aggregate/mutation_rate.tsv",
    qc="out/aggregate/qc.summary.tsv",
    selected_somatic_variants="out/aggregate/mutect2.filter.genes_of_interest.combined.tsv",
    all_somatic_variants="out/aggregate/mutect2.filter.combined.tsv.gz"

  output:
    md="out/aggregate/final.md",
    html="out/aggregate/final.html"

  log:
    stderr="log/make_report.stderr"

  shell:
    "src/make_report.py {config[make_report_params]} --versions {input.versions} --signatures {input.signatures} --burden {input.burden} --qc {input.qc} --selected_somatic_variants {input.selected_somatic_variants} --all_somatic_variants {input.all_somatic_variants} > {output.md} 2>{log.stderr} && "
    "{config[module_pandoc]} && "
    "pandoc {output.md} | src/style_report.py > {output.html}"

### QC ###

rule qc_summary:
  input:
    expand("out/{sample}.artifact_metrics.txt.error_summary_metrics", sample=samples['samples'])
  output:
    "out/aggregate/qc.summary.tsv"
  log:
    stderr="log/make_summary.stderr"
  shell:
    "python src/make_summary.py --verbose --samples {input} > {output} 2>{log.stderr}"

rule fastqc:
  input:
    fastqs=lambda wildcards: samples["samples"][wildcards.sample]
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

rule qc_conpair_matrix:
  input:
    reference=config["genome"],
    tumours=expand("out/{tumour}.sorted.dups.bam", tumour=samples['tumours']),
    germlines=expand("out/{germline}.sorted.dups.bam", germline=germline_samples())
  output:
    "out/aggregate/concordance_matrix.tsv"
  log:
    stdout="log/conpair.stdout",
    stderr="log/conpair.stderr"
  shell:
    "( "
    "{config[module_java]} && "
    "mkdir -p tmp/conpair && "
    "python src/conpair_matrix.py --tumours {input.tumours} --germlines {input.germlines} --reference {input.reference} --working tmp/conpair --output {output} && "
    "python tools/plotme-{config[plotme_version]}/plotme/heatmap.py --x Tumour --y Normal --z Concordance --fontsize 6 --x_rotation vertical --target {output}.png < {output}"
    ") 1>{log.stdout} 2>{log.stderr}"


#rule qc_verifybamid_tumour:
#  input:
#    vcf="out/{tumour}.strelka.somatic.snvs.vcf.gz",
#    bam="out/{tumour}.sorted.dups.bam",
#    bai="out/{tumour}.sorted.dups.bai",
#  output:
#    "out/{tumour}.verifybamid.somatic.completed"
#  log:
#    stderr="log/{tumour}.verifybamid.stderr"
#  shell:
#    "tools/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID --noPhoneHome --vcf {input.vcf} --bam {input.bam} --bai {input.bai} --out out/{wildcards.tumour}.verifybamid --verbose 2>{log.stderr} && touch {output}"

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
    expand("out/fastqc/{sample}/completed", sample=samples['samples']),
    expand("out/{sample}.metrics.alignment", sample=samples['samples']),
    expand("out/{sample}.metrics.insertsize", sample=samples['samples']),
    expand("out/{sample}.metrics.target", sample=samples['samples']),
    expand("out/{tumour}.concordance", tumour=samples['tumours']), # TODO
    expand("out/{tumour}.contamination", tumour=samples['tumours']), # TODO
    # -- expand("out/{tumour}.verifybamid.somatic.completed", tumour=samples['tumours']),
    expand("out/{sample}.depth_of_coverage.sample_summary", sample=samples['samples']),
    "out/aggregate/qc.summary.tsv"
    
  output:
    "out/aggregate/multiqc.html"
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
    "bedtools sort -faidx reference/genome.lengths -i {input.bed} | bedtools merge -i - | bedtools coverage -sorted -hist -b {input.bam} -a stdin -g reference/genome.lengths | grep ^all > {output.hist}"

rule qc_on_target_coverage_hist_dedup:
  input:
    bed=config["regions"],
    bam="out/{sample}.sorted.dups.bam"
  output:
    hist="out/{sample}.ontarget.dedup.hist",
  shell:
    "{config[module_samtools]} && "
    "samtools view -uF 0x400 {input.bam} > tmp/{wildcards.sample}.dedup.bam && "
    "{config[module_bedtools]} && "
    "bedtools sort -faidx reference/genome.lengths -i {input.bed} | bedtools merge -i - > tmp/{wildcards.sample}.qc_on_target_coverage_hist_dedup.bed && "
    "bedtools coverage -sorted -hist -b tmp/{wildcards.sample}.dedup.bam -a tmp/{wildcards.sample}.qc_on_target_coverage_hist_dedup.bed -g reference/genome.lengths > tmp/{wildcards.sample}.qc_on_target_coverage_hist_dedup.hist && "
    "grep ^all tmp/{wildcards.sample}.qc_on_target_coverage_hist_dedup.hist > {output.hist} && "
    "rm tmp/{wildcards.sample}.qc_on_target_coverage_hist_dedup.bed tmp/{wildcards.sample}.dedup.bam"

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

rule qc_on_target_coverage_plot_dedup:
  input:
    samples=expand("out/{sample}.ontarget.dedup.hist", sample=samples['samples']),
    germline=expand("out/{germline}.ontarget.dedup.hist", germline=germline_samples()),
    tumour=expand("out/{tumour}.ontarget.dedup.hist", tumour=samples['tumours'])
  output:
    samples="out/aggregate/ontarget.dedup.png",
    germline="out/aggregate/ontarget_germline.dedup.png",
    tumour="out/aggregate/ontarget_tumour.dedup.png"
  shell:
    "src/plot_coverage.py --target {output.samples} --files {input.samples} --max 10000 && "
    "src/plot_coverage.py --target {output.germline} --files {input.germline} --max 10000 && "
    "src/plot_coverage.py --target {output.tumour} --files {input.tumour} --max 10000"

rule qc_on_target_coverage_plot:
  input:
    expand("out/{sample}.ontarget.hist", sample=samples['samples']),
  output:
    "out/aggregate/ontarget.png"
  shell:
    "src/plot_coverage.py --target {output} --files {input} --max 10000"

rule qc_on_target_coverage_plot_germline:
  input:
    expand("out/{germline}.ontarget.hist", germline=germline_samples()),
  output:
    "out/aggregate/ontarget_germline.png"
  shell:
    "src/plot_coverage.py --target {output} --files {input} --max 10000"

rule qc_on_target_coverage_plot_tumour:
  input:
    expand("out/{tumour}.ontarget.hist", tumour=samples['tumours']),
  output:
    "out/aggregate/ontarget_tumour.png"
  shell:
    "src/plot_coverage.py --target {output} --files {input} --max 10000"

rule qc_on_target_coverage_combined:
  input:
    expand("out/{sample}.ontarget.summary", sample=samples['samples']),
  output:
    "out/aggregate/ontarget.tsv"
  shell:
    "echo \"Filename	n	Mean	Min	Max	Total\" >{output} && "
    "for f in {input}; do echo \"$f	$(tail -1 $f)\" >> {output}; done"

### coverage ###

# exon coverage in genes of interest
rule qc_papillon_exons_min:
  input:
    bams=expand("out/{sample}.sorted.dups.bam", sample=samples['samples']),
    tumours=expand("out/{tumour}.sorted.dups.bam", tumour=samples['tumours']),
    germlines=expand("out/{germline}.sorted.dups.bam", germline=germline_samples())
  output:
    tumour="out/aggregate/gene_coverage_min.tumour.exons.tsv",
    germline="out/aggregate/gene_coverage_min.germline.exons.tsv"
  log:
    stderr="log/papillon_min.exons.stderr"
  params:
    gene_list=' '.join(config["genes_of_interest"])
  shell:
    "(tools/papillon-{config[papillon_version]}/papillon/papillon.py --max_coverage {config[papillon_max_coverage_tumour]} --stat min --bed {config[refseq_bed]} --bams {input.tumours} --genes {params.gene_list} --plot out/aggregate/gene_coverage.tumour.min --padding {config[papillon_padding]} --exon_plots >{output.tumour} && "
    "tools/papillon-{config[papillon_version]}/papillon/papillon.py --max_coverage {config[papillon_max_coverage_germline]} --stat min --bed {config[refseq_bed]} --bams {input.germlines} --genes {params.gene_list} --plot out/aggregate/gene_coverage.germline.min --padding {config[papillon_padding]} --exon_plots >{output.germline}) 2>{log.stderr}"

rule qc_papillon_exons:
  input:
    bams=expand("out/{sample}.sorted.dups.bam", sample=samples['samples']),
    tumours=expand("out/{tumour}.sorted.dups.bam", tumour=samples['tumours']),
    germlines=expand("out/{germline}.sorted.dups.bam", germline=germline_samples())
  output:
    all="out/aggregate/gene_coverage.exons.tsv",
    tumour="out/aggregate/gene_coverage.tumour.exons.tsv",
    germline="out/aggregate/gene_coverage.germline.exons.tsv"
  log:
    stderr="log/papillon.exons.stderr"
  params:
    gene_list=' '.join(config["genes_of_interest"])
  shell:
    "(tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat median --bed {config[refseq_bed]} --bams {input.bams} --genes {params.gene_list} --plot out/aggregate/gene_coverage.exons.median --exon_plots >{output.all} && "
    "tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat median --bed {config[refseq_bed]} --bams {input.tumours} --genes {params.gene_list} --plot out/aggregate/gene_coverage.tumour.exons.median >{output.tumour} && "
    "tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat median --bed {config[refseq_bed]} --bams {input.germlines} --genes {params.gene_list} --plot out/aggregate/gene_coverage.germline.exons.median >{output.germline}) 2>{log.stderr}"

rule qc_papillon_exons_flag:
  input:
    tumour="out/aggregate/gene_coverage_min.tumour.exons.tsv",
    germline="out/aggregate/gene_coverage_min.germline.exons.tsv"
  output:
    tumour="out/aggregate/gene_coverage_min.tumour.exons.qc.flagged.tsv",
    germline="out/aggregate/gene_coverage_min.germline.exons.qc.flagged.tsv"
  shell:
    "python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filter 'Min<{config[papillon_min_coverage_tumour_flag]}' < {input.tumour} > {output.tumour} && "
    "python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filter 'Min<{config[papillon_min_coverage_germline_flag]}' < {input.germline} > {output.germline}"

# region coverage in genes of interest - panel only
rule qc_papillon_regions_min:
  input:
    bams=expand("out/{sample}.sorted.dups.bam", sample=samples['samples']),
    tumours=expand("out/{tumour}.sorted.dups.bam", tumour=samples['tumours']),
    germlines=expand("out/{germline}.sorted.dups.bam", germline=germline_samples())
  output:
    tumour="out/aggregate/gene_coverage_min.tumour.regions.tsv",
    germline="out/aggregate/gene_coverage_min.germline.regions.tsv"
  log:
    stderr="log/papillon_min.regions.stderr"
  shell:
    "(tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat min --bed {config[regions]} --bams {input.tumours} --padding 0 >{output.tumour} && "
    "tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat min --bed {config[regions]} --bams {input.germlines} --padding 0 >{output.germline}) 2>{log.stderr}"

# region coverage in genes of interest - panel only
rule qc_papillon_regions_flag:
  input:
    tumour="out/aggregate/gene_coverage_min.tumour.regions.tsv",
    germline="out/aggregate/gene_coverage_min.germline.regions.tsv"
  output:
    tumour="out/aggregate/gene_coverage_min.tumour.regions.qc.flagged.tsv",
    germline="out/aggregate/gene_coverage_min.germline.regions.qc.flagged.tsv"
  shell:
    "python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filter 'Min<{config[papillon_min_coverage_tumour_flag]}' < {input.tumour} > {output.tumour} && "
    "python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filter 'Min<{config[papillon_min_coverage_germline_flag]}' < {input.germline} > {output.germline}"

### alignment ###
rule trim:
  input:
    fastqs=lambda wildcards: samples["samples"][wildcards.sample]
  output:
    "tmp/{sample}_R1.trimmed.paired.fq.gz",
    "tmp/{sample}_R1.trimmed.unpaired.fq.gz",
    "tmp/{sample}_R2.trimmed.paired.fq.gz",
    "tmp/{sample}_R2.trimmed.unpaired.fq.gz"
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
    fastq_r1="tmp/{sample}_R1.trimmed.paired.fq.gz", 
    fastq_r2="tmp/{sample}_R2.trimmed.paired.fq.gz"
    #fastqs=lambda wildcards: samples["samples"][wildcards.sample]

  output:
    bam="tmp/{sample}.paired.bam",
    sentinel="out/{sample}.paired.bam.completed"

  log:
    "log/{sample}.paired.bwa.log"

  params:
    cores=cluster["align"]["n"],
    read_group=read_group

  shell:
    "{config[module_bwa]} && {config[module_samtools]} && "
    "(bwa mem -M -t {params.cores} -R \"{params.read_group}\" {input.reference} {input.fastq_r1} {input.fastq_r2} | samtools view -b -h -o {output.bam} - && "
    "rm {input.fastq_r1} {input.fastq_r2}) 2>{log} && touch {output.sentinel}"

rule align_unpaired:
  input:
    reference=config["genome"],
    fastq_r1="tmp/{sample}_R1.trimmed.unpaired.fq.gz", 
    fastq_r2="tmp/{sample}_R2.trimmed.unpaired.fq.gz"

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
    "bwa mem -M -t {params.cores} -R \"{params.read_group}\" {input.reference} {input.fastq_r2} | samtools view -b -h -o {output.r2} - && "
    "rm {input.fastq_r1} {input.fastq_r2}) 2>{log}"

# sort the bam
rule merge_bams:
  input:
    "tmp/{sample}.paired.bam",
    "tmp/{sample}_R1.unpaired.bam",
    "tmp/{sample}_R2.unpaired.bam"
  output:
    bam="tmp/{sample}.merged.bam",
    sentinel="out/{sample}.merged.bam.complete"
  log:
    "log/{sample}.merge.log"
  shell:
    "{config[module_samtools]} && "
    "(samtools merge {output.bam} {input} 2>{log}) && "
    "rm {input} && touch {output.sentinel}"

# sort the bam
rule sort:
  input:
    "tmp/{sample}.merged.bam"

  output:
    bam="tmp/{sample}.sorted.bam",
    bai="tmp/{sample}.sorted.bai",
    sentinel="out/{sample}.sorted.bam.sentinel"

  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar SortSam INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=2000000 CREATE_INDEX=True && "
    "rm {input} && touch {output.sentinel}"

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
    "java -jar tools/picard-2.8.2.jar MarkDuplicates INPUT={input} OUTPUT={output[0]} METRICS_FILE={output[2]} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=True CREATE_INDEX=True MAX_RECORDS_IN_RAM=2000000 && "
    "rm {input}"

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
    "tmp/germline_joint_{chromosome}.vcf"
  log:
    "log/gatk_joint_{chromosome}.stderr"
  params:
    variant_list=' '.join(['--variant {}'.format(gvcf) for gvcf in expand("out/{germline}.hc.gvcf.gz", germline=germline_samples())])
  shell:
    "({config[module_java]} && "
    "java -jar tools/GenomeAnalysisTK-3.7.0.jar -T CombineGVCFs -R {input.reference} {params.variant_list} -L {wildcards.chromosome} -o tmp/germline_combined_{wildcards.chromosome}.gvcf && "
    "tools/gatk-4.0.0.0/gatk GenotypeGVCFs -R {input.reference} --dbsnp reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz -V tmp/germline_combined_{wildcards.chromosome}.gvcf -L {wildcards.chromosome} --use-new-qual-calculator true --output tmp/germline_joint_{wildcards.chromosome}.vcf"
    ") 2>{log}"

# call germline variants on the tumour for validation
rule gatk_joint_genotype_tumours:
  input:
    gvcfs=expand("out/{germline}.hc.gvcf.gz", germline=samples["tumours"]),
    reference=config["genome"]
  output:
    "tmp/tumour_joint_{chromosome}.vcf"
  log:
    "log/gatk_joint_tumour_{chromosome}.stderr"
  params:
    variant_list=' '.join(['--variant {}'.format(gvcf) for gvcf in expand("out/{germline}.hc.gvcf.gz", germline=samples["tumours"])])
  shell:
    "({config[module_java]} && "
    "java -jar tools/GenomeAnalysisTK-3.7.0.jar -T CombineGVCFs -R {input.reference} {params.variant_list} -L {wildcards.chromosome} -o tmp/tumours_combined_{wildcards.chromosome}.gvcf && "
    "tools/gatk-4.0.0.0/gatk GenotypeGVCFs -R {input.reference} --dbsnp reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz -V tmp/tumours_combined_{wildcards.chromosome}.gvcf -L {wildcards.chromosome} --use-new-qual-calculator true --output out/tumour_joint_{wildcards.chromosome}.vcf"
    ") 2>{log}"

# notes: 
#   VariantRecalibrator removed due to large cohort size requirements
rule gatk_post_genotype:
  input:
    gvcfs=expand("tmp/germline_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES),
    reference=config["genome"]
  output:
    "out/germline_joint.hc.normalized.vcf"
  log:
    "log/gatk_post.stderr"
  params:
    inputs=' '.join(['--INPUT={}'.format(gvcf) for gvcf in expand("tmp/germline_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES)])
  shell:
    "({config[module_java]} && {config[module_R]} && {config[module_samtools]} && "
    "{config[module_htslib]} && "
    "tools/gatk-4.0.0.0/gatk GatherVcfs -R {input.reference} --OUTPUT=tmp/germline_joint.vcf {params.inputs} && "
    "bgzip -c < tmp/germline_joint.vcf > tmp/germline_joint.vcf.bgz && tabix -p vcf tmp/germline_joint.vcf.bgz && "
    "tools/gatk-4.0.0.0/gatk CalculateGenotypePosteriors -R {input.reference} --supporting reference/gatk-4-bundle-b37/1000G_phase3_v4_20130502.sites.vcf.bgz -V tmp/germline_joint.vcf.bgz -O tmp/germline_joint.cgp.vcf && "
    "tools/vt-{config[vt_version]}/vt normalize -n -r {input.reference} tmp/germline_joint.cgp.vcf -o tmp/germline_joint.cgp.normalized.vcf && "
    "tools/vt-{config[vt_version]}/vt decompose {config[vt_decompose_params]} tmp/germline_joint.cgp.normalized.vcf | tools/vt-{config[vt_version]}/vt normalize -r {input.reference} - -o {output}"
    ") 2>{log}"

rule gatk_post_genotype_tumours:
  input:
    gvcfs=expand("tmp/tumour_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES),
    reference=config["genome"]
  output:
    "out/tumour_joint.hc.normalized.vcf"
  log:
    "log/gatk_post.stderr"
  params:
    inputs=' '.join(['--INPUT={}'.format(gvcf) for gvcf in expand("tmp/tumour_joint_{chromosome}.vcf", chromosome=GATK_CHROMOSOMES)])
  shell:
    "({config[module_java]} && {config[module_R]} && {config[module_samtools]} && "
    "{config[module_htslib]} && "
    "tools/gatk-4.0.0.0/gatk GatherVcfs -R {input.reference} --OUTPUT=tmp/tumour_joint.vcf {params.inputs} && "
    "bgzip -c < tmp/tumour_joint.vcf > tmp/tumour_joint.vcf.bgz && tabix -p vcf tmp/tumour_joint.vcf.bgz && "
    "tools/gatk-4.0.0.0/gatk CalculateGenotypePosteriors -R {input.reference} --supporting reference/gatk-4-bundle-b37/1000G_phase3_v4_20130502.sites.vcf.bgz -V tmp/tumour_joint.vcf.bgz -O tmp/tumour_joint.cgp.vcf && "
    "tools/vt-{config[vt_version]}/vt normalize -n -r {input.reference} tmp/tumour_joint.cgp.vcf -o tmp/tumour_joint.cgp.normalized.vcf && "
    "tools/vt-{config[vt_version]}/vt decompose {config[vt_decompose_params]} tmp/tumour_joint.cgp.normalized.vcf | tools/vt-{config[vt_version]}/vt normalize -r {input.reference} - -o {output}"
    ") 2>{log}"

### qc ###
rule qc_max_coverage_combine:
  input:
    expand("out/{sample}.max_coverage", sample=samples['samples']),
  output:
    "out/aggregate/max_coverage.tsv"
  shell:
    ">{output} && "
    "for f in {input}; do echo \"$f $(grep \"Max coverage\" $f)\" >> {output}; done"

#rule qc_max_trimmed_coverage_combine:
#  input:
#    expand("out/{sample}.max_trimmed_coverage", sample=samples['samples']),
#  output:
#    "out/aggregate/max_trimmed_coverage.tsv"
#  shell:
#    ">{output} && "
#    "for f in {input}; do echo \"$f $(grep \"Max coverage\" $f)\" >> {output}; done"

rule qc_max_coverage:
  input:
    bed=config["regions"],
    fastqs=lambda wildcards: samples["samples"][wildcards.sample]
  output:
    "out/{sample}.max_coverage"
  log:
    "log/{sample}.max_coverage.stderr"
  shell:
    "src/max_coverage.py --verbose --bed {input.bed} --fastqs {input.fastqs} >{output} 2>{log}"

#rule qc_max_trimmed_coverage:
#  input:
#    bed=config["regions"],
#    fastqs=("out/{sample}_R1.trimmed.paired.fq.gz", "out/{sample}_R1.trimmed.unpaired.fq.gz", "out/{sample}_R2.trimmed.paired.fq.gz", "out/{sample}_R2.trimmed.unpaired.fq.gz")
#
#  output:
#    "out/{sample}.max_trimmed_coverage"
#  log:
#    "log/{sample}.max_trimmed_coverage.stderr"
#  shell:
#    "src/max_coverage.py --verbose --bed {input.bed} --fastqs {input.fastqs} >{output} 2>{log}"

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
    "{config[strelka_params]} && "
    "tmp/strelka_{wildcards.tumour}_$$/runWorkflow.py -m local -j {params.cores} && "
    "mv tmp/strelka_{wildcards.tumour}_$$/results/variants/somatic.snvs.vcf.gz {output[0]} && " 
    "mv tmp/strelka_{wildcards.tumour}_$$/results/variants/somatic.indels.vcf.gz {output[1]} && " 
    "rm -r tmp/strelka_{wildcards.tumour}_$$ ) 2>{log}"

rule annotate_af_somatic:
  input:
    vcf="out/{tumour}.strelka.somatic.snvs.vcf.gz",
    bam="out/{tumour}.sorted.dups.bam"
  output:
    "out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
  log:
    stderr="log/{tumour}.annotate_af.stderr"
  shell:
    "({config[module_samtools]} && "
    "{config[module_htslib]} && "
    "{config[module_java]} && "
    "([ ! -e {input.vcf}.tbi ] && tabix -p vcf {input.vcf} || true) && "
    "tools/gatk-4.0.0.0/gatk AnnotateVcfWithBamDepth -V {input.vcf} -I {input.bam} -O tmp/{wildcards.tumour}.strelka.somatic.snvs.af.vcf.gz --lenient && "
    "src/annotate_af.py TUMOR tmp/{wildcards.tumour}.strelka.somatic.snvs.af.vcf.gz | bgzip >{output} ) 2>{log.stderr}"

# tumour only for each germline
rule mutect2_sample_pon_chr:
  input:
    reference=config["genome"],
    bam="out/{germline}.sorted.dups.bam",
    regions=config["regions"]
  output:
    vcf="tmp/{germline}.{chromosome}.mutect2.pon.vcf.gz"
  log:
    stderr="log/{germline}.{chromosome}.mutect2.pon.stderr"
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.0.0.0/gatk Mutect2 -R {input.reference} -I {input.bam} --tumor-sample {wildcards.germline} -L {input.regions} -L {wildcards.chromosome} --interval-set-rule INTERSECTION -O {output.vcf} --interval-padding 1000 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 2>{log.stderr}"

rule mutect2_sample_pon:
  input:
    vcfs=expand("tmp/{{germline}}.{chromosome}.mutect2.pon.vcf.gz", chromosome=GATK_CHROMOSOMES)
  output:
    vcf="out/{germline}.mutect2.pon.vcf.gz"
  log:
    stderr="log/{germline}.mutect2.pon.stderr"
  params:
    inputs=' '.join(['I={}'.format(vcf) for vcf in expand("tmp/{{germline}}.{chromosome}.mutect2.pon.vcf.gz", chromosome=GATK_CHROMOSOMES)])
  shell:
    "{config[module_java]} && "
    "java -jar tools/picard-2.8.2.jar MergeVcfs {params.inputs} O={output.vcf} 2>{log.stderr}"

#rule mutect2_sample_pon:
#  input:
#    reference=config["genome"],
#    bam="out/{germline}.sorted.dups.bam",
#    regions=config["regions"]
#  output:
#    vcf="out/{germline}.mutect2.pon.vcf.gz"
#  log:
#    stderr="log/{germline}.mutect2.pon.stderr"
#  shell:
#    "{config[module_java]} && "
#    "tools/gatk-4.0.0.0/gatk Mutect2 -R {input.reference} -I {input.bam} --tumor-sample {wildcards.germline} -L {input.regions} -O tmp/{wildcards.tumour}.mutect2.pon.vcf.gz --interval-padding 1000 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 2>{log.stderr} && mv tmp/{wildcards.tumour}.mutect2.pon.vcf.gz {output.vcf}"

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
    dbsnp="reference/gatk-4-bundle-b37/dbsnp_138.b37.vcf.bgz",
    pon="out/mutect2.pon.vcf.gz",
    gnomad="reference/af-only-gnomad.raw.sites.b37.vcf.gz",
    bams=tumour_germline_dup_bams
  output:
    "tmp/{tumour}.{chromosome}.mutect2.vcf.gz",
  log:
    stderr="log/{tumour}.{chromosome}.mutect2.stderr"
  params:
    germline=lambda wildcards: samples["tumours"][wildcards.tumour]
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.0.0.0/gatk --java-options '-Xmx30G' Mutect2 -R {config[genome]} -I {input.bams[0]} -I {input.bams[1]} --tumor-sample {wildcards.tumour} --normal-sample {params.germline} --output {output} --output-mode EMIT_VARIANTS_ONLY --dbsnp {input.dbsnp} --germline-resource {input.gnomad} --af-of-alleles-not-in-resource 0.0000025 -pon {input.pon} --interval-padding 1000 -L {config[regions]} -L {wildcards.chromosome} --interval-set-rule INTERSECTION --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter"

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
    "{config[module_htslib]} && "
    "tools/gatk-4.0.0.0/gatk GetPileupSummaries -I {input.bam} -V {input.gnomad} -O tmp/{wildcards.tumour}.mutect2.pileup.table && "
    "tools/gatk-4.0.0.0/gatk CalculateContamination -I tmp/{wildcards.tumour}.mutect2.pileup.table -O tmp/{wildcards.tumour}.mutect2.contamination.table && "
    "tools/gatk-4.0.0.0/gatk FilterMutectCalls -V {input.vcf} --contamination-table tmp/{wildcards.tumour}.mutect2.contamination.table -O tmp/{wildcards.tumour}.mutect2.filter.vcf.gz && "
    "([ ! -e tmp/{wildcards.tumour}.mutect2.filter.vcf.gz.tbi ] && tabix -p vcf tmp/{wildcards.tumour}.mutect2.filter.vcf.gz || true) && "
    "tools/gatk-4.0.0.0/gatk AnnotateVcfWithBamDepth --lenient -I {input.bam} -V tmp/{wildcards.tumour}.mutect2.filter.vcf.gz -O {output}) 1>{log.stdout} 2>{log.stderr}"

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
    "{config[strelka_params]} && "
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
    vcf="out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz",
    reference=config['genome']
  output:
    "out/{tumour}.strelka.somatic.snvs.af.norm.vep.vcf.gz"
  log:
    "log/{tumour}.vep.log"
  params:
    cores=cluster["annotate_vep_somatic_snvs"]["n"]
  shell:
    "{config[module_samtools]} && "
    "{config[module_htslib]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} . 2>{log}"

rule annotate_vep_somatic_indels:
  input:
    vcf="out/{tumour}.strelka.somatic.indels.norm.vcf.gz",
    bam="out/{tumour}.sorted.dups.bam",
    reference=config['genome']
  output:
    "out/{tumour}.strelka.somatic.indels.norm.vep.vcf.gz"
  log:
    "log/{tumour}.vep.log"
  params:
    cores=cluster["annotate_vep_somatic_indels"]["n"]
  shell:
    "({config[module_samtools]} && "
    "{config[module_java]} && "
    "{config[module_htslib]} && "
    "([ ! -e {input.vcf}.tbi ] && tabix -p vcf {input.vcf} || true) && "
    "tools/gatk-4.0.0.0/gatk AnnotateVcfWithBamDepth --lenient -O tmp/{wildcards.tumour}.strelka.somatic.indels.norm.vcf.gz -I {input.bam} -V {input.vcf} && "
    "src/annotate.sh tmp/{wildcards.tumour}.strelka.somatic.indels.norm.vcf.gz {output} {input.reference} {params.cores} .) 2>{log}"

rule annotate_vep_hc:
  input:
    vcf="out/germline_joint.hc.normalized.vcf",
    reference=config['genome']
  output:
    result="out/aggregate/germline_joint.hc.normalized.annot.vcf.gz"
  log:
    "log/hc.vep.log"
  params:
    cores=cluster["annotate_vep_germline"]["n"]
  shell:
    "{config[module_samtools]} && "
    "{config[module_htslib]} && "
    "src/annotate.sh {input.vcf} tmp/germline_joint.hc.normalized.vep.vcf.gz {input.reference} {params.cores} . 2>{log} && mv tmp/germline_joint.hc.normalized.vep.vcf.gz {output.result}"

rule annotate_vep_mutect2:
  input:
    vcf="out/{tumour}.mutect2.filter.vcf.gz",
    reference=config['genome']
  output:
    "out/{tumour}.mutect2.filter.norm.vep.vcf.gz"
  log:
    "log/{tumour}.mutect2.vep.log"
  params:
    cores=cluster["annotate_vep_mutect2"]["n"]
  shell:
    "{config[module_samtools]} && "
    "{config[module_htslib]} && "
    "{config[module_bedtools]} && "
    "tools/vt-{config[vt_version]}/vt decompose {config[vt_decompose_params]} {input.vcf} | tools/vt-{config[vt_version]}/vt normalize -n -r {input.reference} - -o out/{wildcards.tumour}.mutect2.filter.norm.vcf.gz && "
    "src/annotate.sh out/{wildcards.tumour}.mutect2.filter.norm.vcf.gz {output} {input.reference} {params.cores} . 2>{log}"

##### additional variant annotation #####
# cosmic
# clinvar
# revel

rule annotate_cosmic_mutect2:
  input:
    vcf="out/{tumour}.mutect2.filter.norm.vep.vcf.gz"
  output:
    vcf="out/{tumour}.mutect2.filter.norm.annot.vcf.gz"
  log:
    "log/{tumour}.cosmic.log"
  shell:
    "{config[module_htslib]} && "
    "python src/annotate_cosmic.py --cosmic {config[cosmic_counts]} < {input.vcf} | "
    "bgzip > {output.vcf} 2>{log}"

#    "tools/vcfanno_linux64 -lua cfg/vcfanno.lua cfg/vcfanno.cfg {input.vcf} | "

rule annotate_cosmic_strelka_snvs:
  input:
    vcf="out/{tumour}.strelka.somatic.snvs.af.norm.vep.vcf.gz"
  output:
    vcf="out/{tumour}.strelka.somatic.snvs.af.norm.annot.vcf.gz"
  log:
    "log/{tumour}.cosmic.log"
  shell:
    "{config[module_htslib]} && "
    "python src/annotate_cosmic.py --cosmic {config[cosmic_counts]} < {input.vcf} | "
    "bgzip > {output.vcf} 2>{log}"

rule annotate_cosmic_strelka_indels:
  input:
    vcf="out/{tumour}.strelka.somatic.indels.norm.vep.vcf.gz"
  output:
    vcf="out/{tumour}.strelka.somatic.indels.norm.annot.vcf.gz"
  log:
    "log/{tumour}.cosmic.log"
  shell:
    "{config[module_htslib]} && "
    "python src/annotate_cosmic.py --cosmic {config[cosmic_counts]} < {input.vcf} | "
    "bgzip > {output.vcf} 2>{log}"

rule annotate_clinvar:
  input:
    vcfs_mutect2=expand("out/{tumour}.mutect2.filter.norm.annot.revel.vcf.gz", tumour=samples['tumours']),
    vcfs_strelka=expand("out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.vcf.gz", tumour=samples['tumours']),
    vcfs_indels=expand("out/{tumour}.strelka.somatic.indels.norm.annot.revel.vcf.gz", tumour=samples['tumours']),
    vcfs_gl="out/aggregate/germline_joint.hc.normalized.annot.revel.vcf.gz"
  output:
    vcfs_mutect2=expand("out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.vcf.gz", tumour=samples['tumours']),
    vcfs_strelka=expand("out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.clinvar.vcf.gz", tumour=samples['tumours']),
    vcfs_indels=expand("out/{tumour}.strelka.somatic.indels.norm.annot.revel.clinvar.vcf.gz", tumour=samples['tumours']),
    vcfs_gl="out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.vcf.gz"
  log:
    "log/clinvar.log"
  shell:
    "src/annotate_vcf.py "
      "--vcf {config[clinvar]} "
      "--fields CLNDN CLNSIG "
      "--vcfs {input.vcfs_mutect2} {input.vcfs_strelka} {input.vcfs_indels} {input.vcfs_gl} "
      "--suffix clinvar 2>{log}"

rule annotate_revel:
  input:
    vcfs_mutect2=expand("out/{tumour}.mutect2.filter.norm.annot.vcf.gz", tumour=samples['tumours']),
    vcfs_strelka=expand("out/{tumour}.strelka.somatic.snvs.af.norm.annot.vcf.gz", tumour=samples['tumours']),
    vcfs_indels=expand("out/{tumour}.strelka.somatic.indels.norm.annot.vcf.gz", tumour=samples['tumours']),
    vcfs_gl="out/aggregate/germline_joint.hc.normalized.annot.vcf.gz"
  output:
    vcfs_mutect2=expand("out/{tumour}.mutect2.filter.norm.annot.revel.vcf.gz", tumour=samples['tumours']),
    vcfs_strelka=expand("out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.vcf.gz", tumour=samples['tumours']),
    vcfs_indels=expand("out/{tumour}.strelka.somatic.indels.norm.annot.revel.vcf.gz", tumour=samples['tumours']),
    vcfs_gl="out/aggregate/germline_joint.hc.normalized.annot.revel.vcf.gz"
  log:
    "log/revel.log"
  shell:
    "src/annotate_vcf.py "
      "--vcf {config[revel]} "
      "--is_tsv --tsv_zipped --tsv_chrom_column chr --tsv_pos_column hg19_pos --tsv_ref_column ref --tsv_alt_column alt --fields REVEL "
      "--vcfs {input.vcfs_mutect2} {input.vcfs_strelka} {input.vcfs_indels} {input.vcfs_gl} "
      "--suffix revel 2>{log}"

rule annotate_cadd:
  input:
    vcfs_mutect2=expand("out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.vcf.gz", tumour=samples['tumours']),
    vcfs_strelka=expand("out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.clinvar.vcf.gz", tumour=samples['tumours']),
    vcfs_indels=expand("out/{tumour}.strelka.somatic.indels.norm.annot.revel.clinvar.vcf.gz", tumour=samples['tumours']),
    vcfs_gl="out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.vcf.gz",
    cadd="out/cadd.vcf.gz"
  output:
    vcfs_mutect2=expand("out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.cadd.vcf.gz", tumour=samples['tumours']),
    vcfs_strelka=expand("out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.clinvar.cadd.vcf.gz", tumour=samples['tumours']),
    vcfs_indels=expand("out/{tumour}.strelka.somatic.indels.norm.annot.revel.clinvar.cadd.vcf.gz", tumour=samples['tumours']),
    vcfs_gl="out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.cadd.vcf.gz"
  log:
    "log/cadd.log"
  shell:
    "src/annotate_vcf.py "
      "--vcf {input.cadd} "
      "--fields phred --rename 'phred=cadd_phred' "
      "--vcfs {input.vcfs_mutect2} {input.vcfs_strelka} {input.vcfs_indels} {input.vcfs_gl} "
      "--suffix cadd 2>{log}"

rule filter_germline:
  input:
    "out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.cadd.tsv.gz"
  output:
    "out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.genes_of_interest.tsv.gz"
  params:
    region_names=' '.join(['vep_SYMBOL={}'.format(x) for x in config["genes_of_interest"]])
  shell:
    "gunzip < {input} | "
    "csvfilter.py --delimiter '	' --filter {params.region_names} | "
    "csvcols.py --delimiter '	' --cols VCF_SAMPLE_ID CHROM POS ID REF ALT vep_SYMBOL vep_HGVSc vep_HGVSp AD DP CLNDN CLNSIG vep_Consequence vep_IMPACT vep_gnomAD_AF vep_SIFT vep_PolyPhen REVEL | "
    "csvfilter.py --delimiter '	' --filter 'DP!DP' | "
    "src/ad_to_af.py | "
    "csvmap.py --delimiter '	' --map 'vep_gnomAD_AF,,0' | "
    "csvfilter.py --delimiter '	' --filter 'DP>30' 'VAF>0.05' 'vep_gnomAD_AF<0.01' | "
    "gzip > tmp/germline_joint.hc.normalized.annot.revel.clinvar.genes_of_interest.tsv.gz && mv tmp/germline_joint.hc.normalized.annot.revel.clinvar.genes_of_interest.tsv.gz {output}"


#####

rule pass_mutect2:
  input:
    "out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.cadd.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.norm.annot.pass.vcf.gz"
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input} | egrep '(^#|PASS|	str_contraction	)' | bgzip > {output}"

rule dp_af_filter_mutect2:
  input:
    "out/{tumour}.mutect2.filter.norm.annot.pass.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.norm.annot.pass.dp_af.vcf.gz"
  log:
    stderr="log/{tumour}.dp_af_filter_mutect2.log"
  params:
    af=config["af_threshold"],
    dp=config["dp_threshold"],
    tumour="{tumour}"
  shell:
    "{config[module_htslib]} && "
    "src/filter_af.py --af {params.af} --dp {params.dp} --dp_field BAM_DEPTH --sample {params.tumour} < {input} 2>{log.stderr} | bgzip > {output}"

rule pass_strelka_indels:
  input:
    "out/{tumour}.strelka.somatic.indels.norm.annot.revel.clinvar.cadd.vcf.gz",
  output:
    "out/{tumour}.strelka.somatic.indels.norm.annot.pass.vcf.gz",
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input} | egrep '(^#|PASS|	LowDepth	)' | bgzip > {output}"

rule pass_strelka_snvs:
  input:
    "out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.clinvar.cadd.vcf.gz",
  output:
    "out/{tumour}.strelka.somatic.snvs.af.norm.annot.pass.vcf.gz",
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input} | egrep '(^#|PASS|	LowDepth	)' | bgzip > {output}"

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
    "{config[module_htslib]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} . 2>{log}"

rule strelka_normalise:
  input:
    reference=config["genome"],
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.vcf.gz" 
  output:
    snvs_norm="out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz",
    indels_norm="out/{tumour}.strelka.somatic.indels.norm.vcf.gz"
  shell:
    "{config[module_samtools]} && "
    "{config[module_htslib]} && "
    "tools/vt-{config[vt_version]}/vt decompose {config[vt_decompose_params]} {input.strelka_snvs} | tools/vt-{config[vt_version]}/vt normalize -n -r {input.reference} - -o {output.snvs_norm} && "
    "tools/vt-{config[vt_version]}/vt decompose {config[vt_decompose_params]} {input.strelka_indels} | tools/vt-{config[vt_version]}/vt normalize -n -r {input.reference} - -o {output.indels_norm}"
 
rule intersect_somatic_callers:
  input:
    reference=config["genome"],
    mutect2="out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.vcf.gz",
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.norm.vcf.gz" 
  output:
    "out/{tumour}.intersect.vcf.gz"
  log:
    stderr="log/{tumour}.intersect.log"
  shell:
    "({config[module_samtools]} && "
    "{config[module_bedtools]} && "
    "{config[module_htslib]} && "
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_snvs} {input.mutect2} > tmp/{wildcards.tumour}.intersect.unsorted.vcf && "
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_indels} {input.mutect2} | sed -n '/^#/!p' >> tmp/{wildcards.tumour}.intersect.unsorted.vcf && "
    "grep '^#' tmp/{wildcards.tumour}.intersect.unsorted.vcf > tmp/{wildcards.tumour}.intersect.vcf && "
    "bedtools sort -faidx reference/genome.lengths < tmp/{wildcards.tumour}.intersect.unsorted.vcf >> tmp/{wildcards.tumour}.intersect.vcf && "
    "bgzip < tmp/{wildcards.tumour}.intersect.vcf > {output}"
    ") 2>{log.stderr}"
 
rule pass_one_somatic_callers:
  input:
    reference=config["genome"],
    mutect2="out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.vcf.gz",
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.norm.vcf.gz" 
  output:
    "out/{tumour}.pass_one.vcf.gz"
  log:
    stderr="log/{tumour}.pass_one.log"
  shell:
    "({config[module_samtools]} && "
    "{config[module_bedtools]} && "
    "{config[module_htslib]} && "
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_snvs} {input.mutect2} --pass_one > tmp/{wildcards.tumour}.pass_one.unsorted.vcf && " # snvs vs mutect
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_indels} {input.mutect2} --pass_one | sed -n '/^#/!p' >> tmp/{wildcards.tumour}.pass_one.unsorted.vcf && " # indels vs mutect
    "grep '^#' tmp/{wildcards.tumour}.pass_one.unsorted.vcf > tmp/{wildcards.tumour}.pass_one.vcf && " # header
    "bedtools sort -faidx reference/genome.lengths < tmp/{wildcards.tumour}.pass_one.unsorted.vcf >> tmp/{wildcards.tumour}.pass_one.vcf && "
    "bgzip < tmp/{wildcards.tumour}.pass_one.vcf > {output}"
    ") 2>{log.stderr}"
 
rule intersect_pass_somatic_callers:
  input:
    reference=config["genome"],
    mutect2="out/{tumour}.mutect2.filter.norm.annot.pass.vcf.gz",
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.annot.pass.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.norm.annot.pass.vcf.gz" 
  output:
    "out/{tumour}.intersect.pass.vcf.gz"
  log:
    stderr="log/{tumour}.intersect.log"
  shell:
    "({config[module_samtools]} && "
    "{config[module_bedtools]} && "
    "{config[module_htslib]} && "
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_snvs} {input.mutect2} > tmp/{wildcards.tumour}.intersect.unsorted.pass.vcf && " # snvs vs mutect
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_indels} {input.mutect2} | sed -n '/^#/!p' >> tmp/{wildcards.tumour}.intersect.unsorted.pass.vcf && " # indels vs mutect
    "grep '^#' tmp/{wildcards.tumour}.intersect.unsorted.pass.vcf > tmp/{wildcards.tumour}.intersect.pass.vcf && " # header
    "bedtools sort -faidx reference/genome.lengths < tmp/{wildcards.tumour}.intersect.unsorted.pass.vcf >> tmp/{wildcards.tumour}.intersect.pass.vcf && "
    "bgzip < tmp/{wildcards.tumour}.intersect.pass.vcf > {output}"
    ") 2>{log.stderr}"

rule filter_intersected_somatic_callers:
  input:
    nopass="out/{tumour}.intersect.vcf.gz",
    passed="out/{tumour}.intersect.pass.vcf.gz"
  output:
    nopass="out/{tumour}.intersect.filter.vcf.gz",
    passed="out/{tumour}.intersect.pass.filter.vcf.gz"
  log:
    stderr="log/{tumour}.intersect_filter.log"
  params:
    af=config["af_threshold"],
    dp=config["dp_threshold"],
    tumour="{tumour}"
  shell:
    "{config[module_htslib]} && "
    "src/filter_af.py --sample {params.tumour} --af {params.af} --dp {params.dp} < {input.nopass} 2>{log.stderr} | bgzip > {output.nopass} && "
    "src/filter_af.py --sample {params.tumour} --af {params.af} --dp {params.dp} < {input.passed} 2>{log.stderr} | bgzip > {output.passed}"

rule filter_capture:
  input:
    vcf_intersect="out/{tumour}.intersect.pass.filter.vcf.gz",
    vcf_indels="out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz",
    regions=config["regions"]
  output:
    vcf_intersect="out/{tumour}.intersect.pass.filter.capture.vcf.gz",
    vcf_indels="out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.capture.vcf.gz"
  shell:
    "{config[module_bedtools]} && "
    "{config[module_htslib]} && "
    "bedtools intersect -a {input.vcf_intersect} -b {input.regions} -header | bgzip > {output.vcf_intersect} && "
    "bedtools intersect -a {input.vcf_indels} -b {input.regions} -header | bgzip > {output.vcf_indels}"

# mutational signatures with filtered counts
rule mutational_signature_capture_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.filter.capture.vcf.gz",
    vcf_indels="out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz"
  output:
    counts="out/{tumour}.mutational_signature_v3.2_capture.filter.counts",
    counts_indels="out/{tumour}.mutational_signature_v3.2_id_capture.filter.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs_capture.filter.exposures",
    sbs_pks="out/{tumour}.mutational_signature_v3.2_sbs_pks_capture.filter.exposures",
    dbs="out/{tumour}.mutational_signature_v3.2_dbs_capture.filter.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id_capture.filter.exposures",
    id_strelka="out/{tumour}.mutational_signature_v3.2_id_strelka_capture.filter.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_capture.stderr" # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --doublets --genome {input.reference} --vcf {input.vcf} > {output.counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --doublets --genome {input.reference} --vcf {input.vcf_indels} > {output.counts_indels} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures /home/peter/crc/peter/analyses/yocrc/ids_with_pks.tsv --counts {output.counts_indels} > {output.id_strelka} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts {output.counts} > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures /home/peter/crc/peter/analyses/yocrc/sbs_with_pks.tsv --counts {output.counts} > {output.sbs_pks} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts {output.counts} > {output.id} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.1_db.txt --counts {output.counts} > {output.dbs}) 2>{log.stderr}"



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


rule annotate_af_indels:
  input:
    "out/{tumour}.strelka.somatic.indels.norm.annot.pass.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.vcf.gz"
  log:
    stderr="log/{tumour}.strelka.indels.annotate_af.stderr"
  shell:
    "src/annotate_indel_af.py --verbose --sample TUMOR < {input} > {output} 2>{log.stderr}"

rule filter_strelka_snvs:
  input:
    "out/{tumour}.strelka.somatic.snvs.af.norm.annot.pass.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.snvs.af.norm.annot.pass.filter.vcf.gz"
  log:
    stderr="log/{tumour}.strelka.snvs.filter.stderr"
  shell:
    "src/filter_af.py --verbose --sample TUMOR --info_af --af {config[af_threshold]} --dp {config[dp_threshold]} < {input} > {output} 2>{log.stderr}"

rule filter_indels:
  input:
    "out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz"
  log:
    stderr="log/{tumour}.strelka.indels.filter.stderr"
  shell:
    "src/filter_af.py --verbose --sample TUMOR --info_af --af {config[af_threshold]} --dp {config[dp_threshold]} < {input} > {output} 2>{log.stderr}"

rule bias_filter_strelka_indels:
  input:
    reference=config["genome"],
    bam="out/{tumour}.sorted.dups.bam",
    vcf="out/{tumour}.strelka.somatic.indels.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.indels.bias.vcf.gz"
  log:
    stderr="log/{tumour}.bias_filter.indels.bias.err",
    stdout="log/{tumour}.bias_filter.indels.bias.out"
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input.vcf} | egrep '(^#|PASS)' > tmp/{wildcards.tumour}_bias_filter_strelka.vcf && "
    "python tools/DKFZBiasFilter/scripts/biasFilter.py --tempFolder tmp tmp/{wildcards.tumour}_bias_filter_strelka_indels.vcf {input.bam} {input.reference} tmp/{wildcards.tumour}_bias_filter_out_strelka_indels.vcf 2>{log.stderr} 1>{log.stdout} && "
    "bgzip < tmp/{wildcards.tumour}_bias_filter_out_strelka_indels.vcf > {output} && "
    "rm tmp/{wildcards.tumour}_bias_filter_strelka_indels.vcf tmp/{wildcards.tumour}_bias_filter_out_strelka_indels.vcf"

# TODO deprecated
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

rule combine_intersect_tsv:
  input:
    expand("out/{tumour}.intersect.pass.filter.annot.tsv", tumour=samples['tumours'])
  output:
    "out/aggregate/intersect.filter.combined.tsv"
  shell:
    "src/combine_tsv_raw.py {input} | sed 's/^out\\/\\([^.]*\\)\\.[^\\t]*/\\1/' > {output}"

rule intersect_tsv:
  input:
    vcf="out/{tumour}.intersect.pass.filter.signatures.vcf.gz"
  output:
    "out/{tumour}.intersect.pass.filter.annot.tsv"
  shell:
    "src/vcf2tsv.py {input.vcf} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --transcript CANONICAL=YES --override 'POLD1=Feature|NM_002691.4' 'BRAF=Feature|NM_004333.6' >{output}"

rule combine_mutect2_tsv:
  input:
    expand("out/{tumour}.mutect2.filter.annot.tsv", tumour=samples['tumours'])
  output:
    "out/aggregate/mutect2.filter.combined.tsv.gz"
  shell:
    "src/combine_tsv_raw.py {input} | sed 's/^out\\/\\([^.]*\\)\\.[^\\t]*/\\1/' | gzip > {output}"

rule combine_strelka_tsv:
  input:
    snvs=expand("out/{tumour}.strelka.snvs.annot.tsv", tumour=samples['tumours']),
    indels=expand("out/{tumour}.strelka.indels.annot.tsv", tumour=samples['tumours'])
  output:
    "out/aggregate/strelka.all.tsv.gz"
  shell:
    "csvmerge.py --delimiter '	' --union {input.snvs} {input.indels} | gzip>{output}"
    # TODO can remove inputs

# .strelka.somatic.snvs.af.norm.annot.revel.clinvar.cadd
rule strelka_tsv:
  input:
    snvs="out/{tumour}.strelka.somatic.snvs.af.norm.annot.pass.vcf.gz",
    indels="out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.vcf.gz"
  output:
    snvs="out/{tumour}.strelka.snvs.annot.tsv",
    indels="out/{tumour}.strelka.indels.annot.tsv"
  shell:
    "src/vcf2tsv.py --keep_rejected_calls {input.snvs} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --transcript CANONICAL=YES --override 'POLD1=Feature|NM_002691.4' 'BRAF=Feature|NM_004333.6' "
    "| csvfilter.py --delimiter '	' --filter VCF_SAMPLE_ID=TUMOR | csvmap.py --delimiter '	' --map VCF_SAMPLE_ID,TUMOR,{wildcards.tumour} "
    ">{output.snvs} "
    "&& src/vcf2tsv.py --keep_rejected_calls {input.indels} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --transcript CANONICAL=YES --override 'POLD1=Feature|NM_002691.4' 'BRAF=Feature|NM_004333.6' "
    "| csvfilter.py --delimiter '	' --filter VCF_SAMPLE_ID=TUMOR | csvmap.py --delimiter '	' --map VCF_SAMPLE_ID,TUMOR,{wildcards.tumour} "
    ">{output.indels}"

rule mutect2_tsv:
  input:
    vcf="out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.cadd.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.annot.tsv"
  shell:
    "src/vcf2tsv.py {input.vcf} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --transcript CANONICAL=YES --override 'POLD1=Feature|NM_002691.4' 'BRAF=Feature|NM_004333.6' "
    "| src/ad_to_af.py "
    ">{output}"

rule germline_tsv:
  input:
    vcf="out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.cadd.vcf.gz"
  output:
    "out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.cadd.tsv.gz"
  shell:
    "src/vcf2tsv.py {input.vcf} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --transcript CANONICAL=YES --override 'POLD1=Feature|NM_002691.4' 'BRAF=Feature|NM_004333.6' | python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filters 'GT!0/0' 'GT!./.' | gzip >{output}"

rule combine_genes_of_interest:
  input:
    expand("out/{tumour}.mutect2.filter.genes_of_interest.tsv", tumour=samples['tumours']),
  output:
    "out/aggregate/mutect2.filter.genes_of_interest.combined.tsv"
  shell:
    "src/combine_tsv_raw.py {input} | sed 's/^out\\/\\([^.]*\\)\\.[^\\t]*/\\1/' > {output}"

# filter on genes of interest and convert to tsv
rule filter_genes_of_interest_tumour:
  input:
    vcf="out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.genes_of_interest.tsv"
  log:
    stderr="log/{tumour}.filter_genes_of_interest_tumour.err"
  params:
    gene_list=' '.join(config["genes_of_interest"])
  shell:
    "src/vcf2tsv.py {input.vcf} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --transcript CANONICAL=YES --override 'POLD1=Feature|NM_002691.4' 'BRAF=Feature|NM_004333.6' | "
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
#    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --transcript CANONICAL=YES --override 'POLD1=Feature|NM_002691.4' 'BRAF=Feature|NM_004333.6' | "
#    "src/filter_tsv.py --column vep_SYMBOL --values {params.gene_list} > {output}"

# ----- mutational signatures
# mutational signatures with filtered counts
rule mutational_signature_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.vcf.gz"
  output:
    counts="out/{tumour}.mutational_signature_v3.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs.exposures",
    dbs="out/{tumour}.mutational_signature_v3.2_dbs.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.stderr" # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --doublets --genome {input.reference} --vcf {input.vcf} > {output.counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts {output.counts} > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts {output.counts} > {output.id} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.1_db.txt --counts {output.counts} > {output.dbs}) 2>{log.stderr}"

# mutational signatures with filtered counts
rule mutational_signature_filtered_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    counts="out/{tumour}.mutational_signature_v3.filter.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs.filter.exposures",
    dbs="out/{tumour}.mutational_signature_v3.2_dbs.filter.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id.filter.exposures",
    plot_sbs="out/{tumour}.mutational_signature_v3.2_sbs.filter.png",
    plot_id="out/{tumour}.mutational_signature_v3.2_id.filter.png"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_filter.stderr"
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --doublets --genome {input.reference} --vcf {input.vcf} > {output.counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts {output.counts} > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/plot_counts.py --target {output.plot_sbs} --name {wildcards.tumour} --type sbs < {output.counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/plot_counts.py --target {output.plot_id} --type id --name {wildcards.tumour} < {output.counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts {output.counts} > {output.id} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.1_db.txt --counts {output.counts} > {output.dbs}) 2>{log.stderr}"

# mutational signatures with filtered counts
rule mutational_signature_filtered_v3_2_id:
  input:
    reference=config["genome"],
    vcf_indels="out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz",
    vcf_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.annot.pass.filter.vcf.gz"
  output:
    indel_counts="out/{tumour}.mutational_signature_v3.2_strelka_indels.filter.counts",
    snv_counts="out/{tumour}.mutational_signature_v3.2_strelka_snvs.filter.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs_strelka.filter.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id_strelka.filter.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_strelka.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf_indels} > {output.indel_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf_snvs} > {output.snv_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts {output.snv_counts} > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts {output.indel_counts} > {output.id}) 2>{log.stderr}"

# mutational signatures with mutect
rule mutational_signature_mutect2_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.vcf.gz"
  output:
    indel_counts="out/{tumour}.mutational_signature_v3.2_mutect2_indels.filter.counts",
    snv_counts="out/{tumour}.mutational_signature_v3.2_mutect2_snvs.filter.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs_mutect2.filter.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id_mutect2.filter.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_mutect2.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf} > {output.indel_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} > {output.snv_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3.2_mutect2_snvs.filter.counts > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts out/{wildcards.tumour}.mutational_signature_v3.2_mutect2_indels.filter.counts > {output.id}) 2>{log.stderr}"

rule mutational_signature_mutect2_pass_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.mutect2.filter.norm.annot.pass.vcf.gz",
  output:
    indel_counts="out/{tumour}.mutational_signature_v3.2_mutect2_indels.filter.pass.counts",
    snv_counts="out/{tumour}.mutational_signature_v3.2_mutect2_snvs.filter.pass.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs_mutect2.filter.pass.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id_mutect2.filter.pass.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_mutect2.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf} > {output.indel_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} > {output.snv_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3.2_mutect2_snvs.filter.pass.counts > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts out/{wildcards.tumour}.mutational_signature_v3.2_mutect2_indels.filter.pass.counts > {output.id}) 2>{log.stderr}"

rule mutational_signature_mutect2_pass_dp_af_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.mutect2.filter.norm.annot.pass.dp_af.vcf.gz",
  output:
    indel_counts="out/{tumour}.mutational_signature_v3.2_mutect2_indels.filter.pass.dp_af.counts",
    snv_counts="out/{tumour}.mutational_signature_v3.2_mutect2_snvs.filter.pass.dp_af.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs_mutect2.filter.pass.dp_af.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id_mutect2.filter.pass.dp_af.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_mutect2_dp_af.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf} > {output.indel_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} > {output.snv_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3.2_mutect2_snvs.filter.pass.dp_af.counts > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts out/{wildcards.tumour}.mutational_signature_v3.2_mutect2_indels.filter.pass.dp_af.counts > {output.id}) 2>{log.stderr}"

# mutational signatures with filtered counts
rule mutational_signature_filtered_pks_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz",
    vcf_strelka="out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz"
  output:
    sbs_counts="out/{tumour}.mutational_signature_v3.2_sbs.filter.pks.counts",
    id_counts="out/{tumour}.mutational_signature_v3.2_id.filter.pks.counts",
    id_strelka_counts="out/{tumour}.mutational_signature_v3.2_id_strelka.filter.pks.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs.filter.pks.exposures",
    id="out/{tumour}.mutational_signature_v3.2_id.filter.pks.exposures",
    id_strelka="out/{tumour}.mutational_signature_v3.2_id_strelka.filter.pks.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_filter.pks.stderr"
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} > {output.sbs_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf} > {output.id_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf_strelka} > {output.id_strelka_counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures /home/peter/crc/peter/analyses/yocrc/sbs_with_pks.tsv --counts {output.sbs_counts} > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures /home/peter/crc/peter/analyses/yocrc/ids_with_pks.tsv --counts {output.id_strelka_counts} > {output.id_strelka} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures /home/peter/crc/peter/analyses/yocrc/ids_with_pks.tsv --counts {output.id_counts} > {output.id}) 2>{log.stderr}"

# mutational signatures with filtered counts
rule mutational_signature_v3_2:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.filter.capture.vcf.gz",
    vcf_indels="out/{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz"
  output:
    counts="out/{tumour}.mutational_signature_v3.2.filter.counts",
    counts_indels="out/{tumour}.mutational_signature_v3.2.strelka.counts",
    sbs="out/{tumour}.mutational_signature_v3.2_sbs.filter.capture.exposures",
    dbs="out/{tumour}.mutational_signature_v3.2_db.filter.capture.exposures",
    id_strelka="out/{tumour}.mutational_signature_v3.2_id.strelka.exposures",
    id_capture="out/{tumour}.mutational_signature_v3.2_id.capture.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2.stderr" # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} --indels > {output.counts} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf_indels} > {output.counts_indels} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --counts {output.counts} > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts {output.counts_indels} > {output.id_strelka} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_id.txt --counts {output.counts} > {output.id_capture} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.1_db.txt --counts {output.counts} > {output.dbs}) 2>{log.stderr}"

# mutational signatures with filtered counts
rule mutational_signature_filtered_pks_alexandrov_v3:
  input:
    reference=config["genome"],
    sbs_counts="out/{tumour}.mutational_signature_v3.2_sbs.filter.pks.counts",
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    sbs="out/{tumour}.mutational_signature_v3.2_sbs_alexandrov.filter.pks.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.2_filter.pks.stderr"
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures /home/peter/crc/peter/analyses/yocrc/sbs_with_pks_alexandrov.tsv --counts {input.sbs_counts} > {output.sbs}) 2>{log.stderr}"

# signatures from all samples
#rule combine_mutational_signatures_v3_2: -> combine_mutational_signatures_filtered_v3
#  input:
#    sbs=expand("out/{tumour}.mutational_signature_v3.2_sbs.exposures", tumour=samples['tumours']),
#    dbs=expand("out/{tumour}.mutational_signature_v3.2_db.exposures", tumour=samples['tumours']),
#    id=expand("out/{tumour}.mutational_signature_v3.2_id.strelka.exposures", tumour=samples['tumours'])
#  output:
#    sbs="out/aggregate/mutational_signatures_v3.2_sbs.combined.tsv",
#    dbs="out/aggregate/mutational_signatures_v3.2_db.combined.tsv",
#    id="out/aggregate/mutational_signatures_v3.2_id.combined.tsv"
#  shell:
#    "src/combine_tsv.py {input.sbs} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs\\.exposures/\\1/' >{output.sbs} && "
#    "src/combine_tsv.py {input.dbs} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_db\\.exposures/\\1/' >{output.dbs} && "
#    "src/combine_tsv.py {input.id} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id\\.strelka\\.exposures/\\1/' >{output.id}"

# signatures from all samples
rule combine_mutational_signatures_filtered_v3:
  input:
    sbs=expand("out/{tumour}.mutational_signature_v3.2_sbs.filter.exposures", tumour=samples['tumours']),
    sbs_capture=expand("out/{tumour}.mutational_signature_v3.2_sbs_capture.filter.exposures", tumour=samples['tumours']),
    sbs2=expand("out/{tumour}.mutational_signature_v3.2_sbs_strelka.filter.exposures", tumour=samples['tumours']),
    sbs_mutect2=expand("out/{tumour}.mutational_signature_v3.2_sbs_mutect2.filter.exposures", tumour=samples['tumours']),
    sbs_mutect2_pass=expand("out/{tumour}.mutational_signature_v3.2_sbs_mutect2.filter.pass.exposures", tumour=samples['tumours']),
    sbs_mutect2_pass_dp_af=expand("out/{tumour}.mutational_signature_v3.2_sbs_mutect2.filter.pass.dp_af.exposures", tumour=samples['tumours']),
    sbs_pks=expand("out/{tumour}.mutational_signature_v3.2_sbs.filter.pks.exposures", tumour=samples['tumours']),
    sbs_pks_alexandrov=expand("out/{tumour}.mutational_signature_v3.2_sbs_alexandrov.filter.pks.exposures", tumour=samples['tumours']),
    sbs_pks_capture=expand("out/{tumour}.mutational_signature_v3.2_sbs_pks_capture.filter.exposures", tumour=samples['tumours']),
    dbs=expand("out/{tumour}.mutational_signature_v3.2_dbs.filter.exposures", tumour=samples['tumours']),
    dbs2=expand("out/{tumour}.mutational_signature_v3.2_dbs.exposures", tumour=samples['tumours']),
    id=expand("out/{tumour}.mutational_signature_v3.2_id.exposures", tumour=samples['tumours']),
    id_capture=expand("out/{tumour}.mutational_signature_v3.2_id_capture.filter.exposures", tumour=samples['tumours']),
    id2=expand("out/{tumour}.mutational_signature_v3.2_id_strelka.filter.exposures", tumour=samples['tumours']),
    id_mutect2=expand("out/{tumour}.mutational_signature_v3.2_id_mutect2.filter.exposures", tumour=samples['tumours']),
    id_mutect2_pass=expand("out/{tumour}.mutational_signature_v3.2_id_mutect2.filter.pass.exposures", tumour=samples['tumours']),
    id_mutect2_pass_dp_af=expand("out/{tumour}.mutational_signature_v3.2_id_mutect2.filter.pass.dp_af.exposures", tumour=samples['tumours']),
    id_pks=expand("out/{tumour}.mutational_signature_v3.2_id.filter.pks.exposures", tumour=samples['tumours']),
    id_strelka_pks=expand("out/{tumour}.mutational_signature_v3.2_id_strelka.filter.pks.exposures", tumour=samples['tumours']),
    id_strelka_capture=expand("out/{tumour}.mutational_signature_v3.2_id_strelka_capture.filter.exposures", tumour=samples['tumours'])
  output:
    sbs="out/aggregate/mutational_signatures_v3.2_sbs.filter.combined.tsv",
    sbs_capture="out/aggregate/mutational_signatures_v3.2_sbs_capture.filter.combined.tsv",
    sbs2="out/aggregate/mutational_signatures_v3.2_sbs_strelka.filter.combined.tsv",
    sbs_mutect2="out/aggregate/mutational_signatures_v3.2_sbs_mutect2.filter.combined.tsv",
    sbs_mutect2_pass="out/aggregate/mutational_signatures_v3.2_sbs_mutect2.filter.pass.combined.tsv",
    sbs_mutect2_pass_dp_af="out/aggregate/mutational_signatures_v3.2_sbs_mutect2.filter.pass.dp_af.combined.tsv",
    sbs_pks="out/aggregate/mutational_signatures_v3.2_sbs.filter.pks.combined.tsv",
    sbs_pks_alexandrov="out/aggregate/mutational_signatures_v3.2_sbs_alexandrov.filter.pks.combined.tsv",
    sbs_pks_capture="out/aggregate/mutational_signatures_v3.2_sbs_capture.filter.pks.combined.tsv",
    dbs="out/aggregate/mutational_signatures_v3.2_dbs.filter.combined.tsv",
    dbs2="out/aggregate/mutational_signatures_v3.2_dbs.combined.tsv",
    id="out/aggregate/mutational_signatures_v3.2_id.combined.tsv",
    id_capture="out/aggregate/mutational_signatures_v3.2_id_capture.filter.combined.tsv",
    id2="out/aggregate/mutational_signatures_v3.2_id_strelka.filter.combined.tsv",
    id_mutect2="out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.combined.tsv",
    id_mutect2_pass="out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.pass.combined.tsv",
    id_mutect2_pass_dp_af="out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.pass.dp_af.combined.tsv",
    id_pks="out/aggregate/mutational_signatures_v3.2_id.filter.pks.combined.tsv",
    id_strelka_pks="out/aggregate/mutational_signatures_v3.2_id_strelka.filter.pks.combined.tsv",
    id_strelka_capture="out/aggregate/mutational_signatures_v3.2_id_strelka_capture.filter.pks.combined.tsv"
  shell:
    "src/combine_tsv.py {input.sbs} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs\\.filter\\.exposures/\\1/' >{output.sbs} && "
    "src/combine_tsv.py {input.sbs_capture} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs_capture\\.filter\\.exposures/\\1/' >{output.sbs_capture} && "
    "src/combine_tsv.py {input.sbs2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs_strelka\\.filter\\.exposures/\\1/' >{output.sbs2} && "
    "src/combine_tsv.py {input.sbs_mutect2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs_mutect2.filter\\.exposures/\\1/' >{output.sbs_mutect2} && "
    "src/combine_tsv.py {input.sbs_mutect2_pass} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs_mutect2.filter.pass\\.exposures/\\1/' >{output.sbs_mutect2_pass} && "
    "src/combine_tsv.py {input.sbs_mutect2_pass_dp_af} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs_mutect2.filter.pass.dp_af\\.exposures/\\1/' >{output.sbs_mutect2_pass_dp_af} && "
    "src/combine_tsv.py {input.sbs_pks} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs\\.filter\\.pks\\.exposures/\\1/' >{output.sbs_pks} && "
    "src/combine_tsv.py {input.sbs_pks_alexandrov} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs_alexandrov\\.filter\\.pks\\.exposures/\\1/' >{output.sbs_pks_alexandrov} && "
    "src/combine_tsv.py {input.sbs_pks_capture} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs_capture\\.filter\\.pks\\.exposures/\\1/' >{output.sbs_pks_capture} && "
    "src/combine_tsv.py {input.dbs} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_dbs\\.filter\\.exposures/\\1/' >{output.dbs} && "
    "src/combine_tsv.py {input.dbs2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_dbs\\.exposures/\\1/' >{output.dbs2} && "
    "src/combine_tsv.py {input.id} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id\\.exposures/\\1/' >{output.id} && "
    "src/combine_tsv.py {input.id_capture} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id_capture\\.filter\\.exposures/\\1/' >{output.id_capture} && "
    "src/combine_tsv.py {input.id2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id_strelka\\.filter\\.exposures/\\1/' >{output.id2} && "
    "src/combine_tsv.py {input.id_mutect2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id_mutect2\\.filter\\.exposures/\\1/' >{output.id_mutect2} && "
    "src/combine_tsv.py {input.id_mutect2_pass} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id_mutect2\\.filter\\.pass\\.exposures/\\1/' >{output.id_mutect2_pass} && "
    "src/combine_tsv.py {input.id_mutect2_pass_dp_af} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id_mutect2\\.filter\\.pass.dp_af\\.exposures/\\1/' >{output.id_mutect2_pass_dp_af} && "
    "src/combine_tsv.py {input.id_pks} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id\\.filter\\.pks\\.exposures/\\1/' >{output.id_pks} && "
    "src/combine_tsv.py {input.id_strelka_pks} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id_strelka\\.filter\\.pks\\.exposures/\\1/' >{output.id_strelka_pks} && "
    "src/combine_tsv.py {input.id_strelka_capture} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_id_strelka_capture\\.filter\\.pks\\.exposures/\\1/' >{output.id_strelka_capture}"

# mutational signature contexts from samples
rule combine_mutational_signatures_contexts_v3:
  input:
    sbs=expand("out/{tumour}.mutational_signature_v3.filter.counts", tumour=samples['tumours']),
    sbs2=expand("out/{tumour}.mutational_signature_v3.2_strelka_snvs.filter.counts", tumour=samples['tumours']),
    sbs_mutect2=expand("out/{tumour}.mutational_signature_v3.2_mutect2_snvs.filter.counts", tumour=samples['tumours']),
    id=expand("out/{tumour}.mutational_signature_v3.2_strelka_indels.filter.counts", tumour=samples['tumours']),
    id_mutect2=expand("out/{tumour}.mutational_signature_v3.2_mutect2_indels.filter.counts", tumour=samples['tumours']),
    id_mutect2_pass=expand("out/{tumour}.mutational_signature_v3.2_mutect2_indels.filter.pass.counts", tumour=samples['tumours']),
    sbs_pks=expand("out/{tumour}.mutational_signature_v3.2_sbs.filter.pks.counts", tumour=samples['tumours'])
  output:
    sbs="out/aggregate/mutational_signatures_v3.2_sbs.filter.contexts.tsv",
    sbs2="out/aggregate/mutational_signatures_v3.2_sbs_strelka.filter.contexts.tsv",
    sbs_mutect2="out/aggregate/mutational_signatures_v3.2_sbs_mutect2.filter.contexts.tsv",
    id="out/aggregate/mutational_signatures_v3.2_id.contexts.tsv",
    id_mutect2="out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.contexts.tsv",
    id_mutect2_pass="out/aggregate/mutational_signatures_v3.2_id_mutect2.filter.pass.contexts.tsv",
    sbs_pks="out/aggregate/mutational_signatures_v3.2_sbs.filter.contexts.pks.tsv"
  shell:
    "src/combine_tsv.py {input.sbs} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs\\.filter\\.counts/\\1/' >{output.sbs} && "
    "src/combine_tsv.py {input.sbs2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_strelka_snvs\\.filter\\.counts/\\1/' >{output.sbs2} && "
    "src/combine_tsv.py {input.sbs_mutect2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_mutect2_snvs.filter\\.counts/\\1/' >{output.sbs_mutect2} && "
    "src/combine_tsv.py {input.id} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_strelka_indels\\.filter\\.counts/\\1/' >{output.id} && "
    "src/combine_tsv.py {input.id_mutect2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_mutect2_indels\\.filter\\.counts/\\1/' >{output.id_mutect2} && "
    "src/combine_tsv.py {input.id_mutect2_pass} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_mutect2_indels\\.filter\\.pass\\.counts/\\1/' >{output.id_mutect2_pass} && "
    "src/combine_tsv.py {input.sbs_pks} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3.2_sbs\\.filter\\.pks.counts/\\1/' >{output.sbs_pks}"

######### extended contexts
rule combine_extended_contexts:
  input:
    flattened=expand("tmp/{tumour}.extended_contexts.v2.flattened", tumour=samples['tumours']),
  output:
    combined="out/aggregate/extended_contexts.v2.tsv"
  shell:
    "src/combine_tsv.py {input.flattened} | sed 's/^tmp\\/\\(.*\)\\.extended_contexts.v2.flattened/\\1/' >{output.combined}"

rule extended_contexts:
  input:
    reference=config["genome"],
    transcripts=config["transcripts"],
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    contexts="out/{tumour}.extended_contexts.v2.tsv",
    detail="out/{tumour}.extended_contexts.detail.tsv",
    flattened="tmp/{tumour}.extended_contexts.v2.flattened"
  log:
    stderr="log/{tumour}.extended_contexts.log"
  shell:
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/extended_context.py --genome {input.reference} --vcf {input.vcf} --transcripts {input.transcripts} --output {output.detail} "
      "--rules 'T>*,-1~AT,+1~AT,-3=A' 'A>*,-1~AT,+1~AT,+3=T' "
      "'T>,-1=A,-2=A,-3=A,-4=A,+1!T' 'A>,+4=T,+3=T,+2=T,+1=T,-1!A' "
      "'T>,-1=A,-2=A,-3=A,-4!A,+1=T,+2!T' 'A>,+4!T,+3=T,+2=T,+1=T,-1=A,-2!A' "
      "'T>,-1=A,-2=A,-3!A,+1=T,+2=T,+3!T' 'A>,+3!T,+2=T,+1=T,-1=A,-2=A,-3!A' "
      "'T>,-1=A,+1=T,+2=T,+3=T' 'A>,+1=T,-1=A,-2=A,-3=A' "
      ">{output.contexts} 2>{log.stderr} && "
    "python tools/csvtools-{config[csvtools_version]}/csvtools/csvflatten.py --delimiter '	' --key Rule < {output.contexts} > {output.flattened}"

# ----- copy number
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

# burden snvs
rule mutation_burden_low_af:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.vcf.gz", tumour=samples['tumours']),
    regions=config["regions"]
  output:
    "out/aggregate/mutation_rate.low_af.tsv"
  log:
    stderr="log/mutation_rate.stderr"
  shell:
    "src/mutation_rate.py --vcfs {input.vcfs} --bed {input.regions} --min_af 0.01 >{output} 2>{log.stderr}"

# burden snvs
rule mutation_burden_high_af:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']),
    regions=config["regions"],
  output:
    "out/aggregate/mutation_rate.high_af.tsv"
  log:
    stderr="log/mutation_rate.stderr"
  shell:
    "src/mutation_rate.py --vcfs {input.vcfs} --bed {input.regions} --min_af 0.2 >{output} 2>{log.stderr}"

# burden snvs
rule mutation_burden:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']),
    regions=config["regions"],
  output:
    "out/aggregate/mutation_rate.tsv"
  log:
    stderr="log/mutation_rate.stderr"
  shell:
    "src/mutation_rate.py --vcfs {input.vcfs} --bed {input.regions} >{output} 2>{log.stderr}"

# burden indels 
rule msi_burden:
  input:
    vcfs=expand("out/{tumour}.strelka.somatic.indels.norm.annot.pass.vcf.gz", tumour=samples['tumours']),
    regions="reference/msi.regions.bed"
  output:
    "out/aggregate/msi_burden.tsv"
  log:
    stderr="log/msi_burden.stderr"
  shell:
    "src/mutation_rate.py --vcfs {input.vcfs} --bed {input.regions} --indels_only --sample_name TUMOR >{output} 2>{log.stderr}"

# af distribution
rule plot_af_strelka:
  input:
    "out/{tumour}.strelka.somatic.snvs.af.vcf.gz"
  output:
    all="out/{tumour}.strelka.somatic.all.af.png",
    just_pass="out/{tumour}.strelka.somatic.pass.af.png"
  shell:
    "src/plot_af.py --log --sample TUMOR --target {output.all} --dp {config[dp_threshold]} --info_af --title 'Variant count as a function of somatic allele fraction for {wildcards.tumour}' < {input} && "
    "src/plot_af.py --just_pass --log --sample TUMOR --target {output.just_pass} --dp {config[dp_threshold]} --info_af --title 'Variant count as a function of somatic allele fraction for {wildcards.tumour}' < {input}"

# af distribution

rule annotate_strelka_signatures:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
    exposures="out/{tumour}.mutational_signature_v3.2_sbs.filter.exposures"
  output:
    vcf="out/{tumour}.strelka.somatic.snvs.af.signatures.vcf.gz",
    plot="out/{tumour}.strelka.somatic.mutational_signature_v3.2_sbs.png"
  shell:
    "{config[module_htslib]} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/annotate_context.py --genome {input.reference} --vcf {input.vcf} > tmp/{wildcards.tumour}.somatic.snvs.af.contexts.vcf.gz && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/assign_signatures.py --definition tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --signatures {input.exposures} --artefacts tools/mutational_signature-{config[signature_version]}/data/signature_summary.tsv --vcf tmp/{wildcards.tumour}.somatic.snvs.af.contexts.vcf.gz --threshold 0.01 --plot {output.plot} > {output.vcf}"
  
rule plot_af_strelka_signatures:
  input:
    vcf="out/{tumour}.strelka.somatic.snvs.af.signatures.vcf.gz",
  output:
    nopass="out/{tumour}.strelka.somatic.signatures.af.png",
    justpass="out/{tumour}.strelka.somatic.pass.signatures.af.png",
    nopass_percent="out/{tumour}.strelka.somatic.signatures.percent.af.png",
    justpass_percent="out/{tumour}.strelka.somatic.pass.signatures.percent.af.png"
  shell:
    "src/plot_af.py --signature_likelihoods --sample TUMOR --target {output.nopass} --dp {config[dp_threshold]} --info_af --title 'Signature breakdown as a function of somatic allele fraction for {wildcards.tumour}' < {input.vcf} && "
    "src/plot_af.py --just_pass --signature_likelihoods --sample TUMOR --target {output.justpass} --dp {config[dp_threshold]} --info_af --title 'Signature breakdown as a function of somatic allele fraction for {wildcards.tumour} PASS variants' < {input.vcf} && "
    "src/plot_af.py --signature_likelihoods --sample TUMOR --target {output.nopass_percent} --dp {config[dp_threshold]} --info_af --percent --title 'Signature proportion as a function of somatic allele fraction for {wildcards.tumour}' < {input.vcf} && "
    "src/plot_af.py --just_pass --signature_likelihoods --sample TUMOR --target {output.justpass_percent} --dp {config[dp_threshold]} --info_af --percent --title 'Signature proportion as a function of somatic allele fraction for {wildcards.tumour} PASS variants' < {input.vcf}"

rule plot_af_mutect2:
  input:
    "out/{tumour}.mutect2.filter.norm.annot.revel.clinvar.vcf.gz"
  output:
    "out/{tumour}.mutect2.somatic.af.png"
  shell:
    "src/plot_af.py --verbose --log --sample {wildcards.tumour} --target {output} --dp {config[dp_threshold]} --vep_format 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL' --genes {config[af_genes]} --impacts HIGH --gene_colors {config[af_gene_colors]} < {input}"

#"src/plot_af.py --log --sample {wildcards.tumour} --target {output} --dp {config[dp_threshold]} --vep_format "Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK|CANONICAL" --genes {config[af_genes]} --consequences {config[af_consequences]} < {input}"

#####################
# loh
# --germline_vcf $OUT/$gl*.strelka.germline.norm.annot.revel.cadd.gnomad.dbnsfp.vcf.gz
#    "out/{germline}.strelka.germline.filter_gt.vcf.gz",
rule loh:
  input:
    snvs="out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.clinvar.vcf.gz", # loh requires strelka for now
    indels="out/{tumour}.strelka.somatic.indels.norm.annot.revel.clinvar.vcf.gz",
    purity="out/aggregate/purity.tsv"
  output:
    tsv="out/{tumour}.loh.tsv",
    bed="out/{tumour}.loh.bed"
  log:
    stderr="log/{tumour}.loh.stderr"
  params:
    tumour="{tumour}",
    regions=' '.join(config["regions_of_interest"]),
    region_names=' '.join(config["genes_of_interest"]),
    region_padding=' '.join([str(x) for x in config["loh_region_padding"]]),
    germline=germline_sample # TODO assumes this file is complete

  shell:
    "(purity=$(python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filter Sample={params.tumour} < {input.purity} | python tools/csvtools-{config[csvtools_version]}/csvtools/csvcols.py --delimiter '	' --cols Best | sed '1d' | tr '\\r\\n' ' ') && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 30 --neutral --min_af 0.1 --tumour_cellularity $purity < {input.snvs} > tmp/{params.tumour}.loh.snvs.tsv && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 30 --neutral --min_af 0.1 --tumour_cellularity $purity < {input.indels} > tmp/{params.tumour}.loh.indels.tsv && "
    "sort -k1,1 -k2,2n tmp/{params.tumour}.loh.snvs.tsv tmp/{params.tumour}.loh.indels.tsv > {output.tsv} && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh_merge.py --verbose --noheader --min_len 1000 --min_prop 0.1 --min_count 2 --max_loh_without_evidence 10000000 --plot out/{params.tumour}.loh --regions {params.regions} --region_names {params.region_names} --region_padding {params.region_padding} --plot_chromosomes --title \"Evidence for LOH in {params.tumour} using estimated purity $purity\" <{output.tsv} >{output.bed} && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh_merge.py --verbose --noheader --min_len 1000 --min_prop 0.1 --min_count 2 --max_loh_without_evidence 10000000 --plot out/{params.tumour}.wide.loh --regions {params.regions} --region_names {params.region_names} --region_padding 50000000 --title \"Evidence for LOH in {params.tumour} using estimated purity $purity\" <{output.tsv} >/dev/null) 2>{log.stderr}"

#####################
rule loh_pass:
  input:
    snvs="out/{tumour}.strelka.somatic.snvs.af.norm.annot.revel.clinvar.vcf.gz", # loh requires strelka for now
    indels="out/{tumour}.strelka.somatic.indels.norm.annot.revel.clinvar.vcf.gz",
    purity="out/aggregate/purity.tsv"
  output:
    tsv="out/{tumour}.pass.loh.tsv",
    bed="out/{tumour}.pass.loh.bed"
  log:
    stderr="log/{tumour}.pass.loh.stderr"
  params:
    tumour="{tumour}",
    regions=' '.join(config["regions_of_interest"]),
    region_names=' '.join(config["genes_of_interest"]),
    region_padding=' '.join([str(x) for x in config["loh_region_padding"]]),
    germline=germline_sample # TODO assumes this file is complete
  shell:
    "(purity=$(python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filter Sample={params.tumour} < {input.purity} | python tools/csvtools-{config[csvtools_version]}/csvtools/csvcols.py --delimiter '	' --cols Best | sed '1d' | tr '\\r\\n' ' ') && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh.py --germline NORMAL --tumour TUMOR --min_dp_germline 10 --min_dp_tumour 30 --neutral --min_af 0.1 --tumour_cellularity $purity < {input.snvs} > tmp/{params.tumour}.pass.loh.snvs.tsv && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh.py --germline NORMAL --tumour TUMOR --min_dp_germline 10 --min_dp_tumour 30 --neutral --min_af 0.1 --tumour_cellularity $purity < {input.indels} > tmp/{params.tumour}.pass.loh.indels.tsv && "
    "sort -k1,1 -k2,2n tmp/{params.tumour}.pass.loh.snvs.tsv tmp/{params.tumour}.pass.loh.indels.tsv > {output.tsv} && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh_merge.py --noheader --min_len 1000 --min_prop 0.1 --min_count 2 --max_loh_without_evidence 10000000 --plot out/{params.tumour}.pass.loh --regions {params.regions} --region_names {params.region_names} --region_padding {params.region_padding} --plot_chromosomes --title \"Evidence for LOH in {params.tumour} using estimated purity $purity and pass variants\" <{output.tsv} >{output.bed} && "
    "tools/LOHdeTerminator-{config[loh_version]}/loh_merge.py --noheader --min_len 1000 --min_prop 0.1 --min_count 2 --max_loh_without_evidence 10000000 --plot out/{params.tumour}.pass.wide.loh --regions {params.regions} --region_names {params.region_names} --region_padding 50000000 --title \"Evidence for LOH in {params.tumour} using estimated purity $purity and pass variants\" <{output.tsv} >/dev/null) 2>{log.stderr}"


rule loh_summary:
  input:
    beds=expand("out/{tumour}.loh.bed", tumour=samples['tumours']),
    transcripts=config["transcripts"],
  output:
    tsv="out/aggregate/loh.genes.tsv"
  log:
    stderr="log/loh_summary.stderr"
  shell:
    "src/combine_loh.py --min_accept 2 --lohs {input.beds} --transcripts {input.transcripts} >{output} 2>{log.stderr}"

#####################
# cnv not working yet

# cnvkit
#rule cnvkit_prep:
#  input:
#    tumours=expand("out/{sample}.sorted.dups.bam", sample=tumour_samples()),
#    germlines=expand("out/{sample}.sorted.dups.bam", sample=germline_samples()),
#    bed=config["regions"],
#    reference=config["genome"]
#  output:
#    "out/cnvkit.reference.cnn"
#  shell:
#    "cnvkit.py access {input.reference} -o tmp/access.hg19.bed && "
#    "cnvkit.py autobin {input.tumours} {input.germlines} -t {input.bed} -g tmp/access.hg19.bed --annotate {config[refflat]} --short-names && "
#    "(for bam in {input.tumours} {input.germlines}; do "
#      "cnvkit.py coverage $bam {input.bed} -o ${bam}.targetcoverage.cnn && "
#      cnvkit.py coverage $bam baits.antitarget.bed -o Sample.antitargetcoverage.cnn

#rule cnvkit_coverage:
# With all normal samples...
#cnvkit.py reference *Normal.{,anti}targetcoverage.cnn --fasta hg19.fa -o my_reference.cnn

# For each tumor sample...
#cnvkit.py fix Sample.targetcoverage.cnn Sample.antitargetcoverage.cnn my_reference.cnn -o Sample.cnr
#cnvkit.py segment Sample.cnr -o Sample.cns

# Optionally, with --scatter and --diagram
#cnvkit.py scatter Sample.cnr -s Sample.cns -o Sample-scatter.pdf
#cnvkit.py diagram Sample.cnr -s Sample.cns -o Sample-diagram.pdf


#rule cnvkit:
#  input:
#    tumours=expand("out/{sample}.sorted.dups.bam", sample=tumour_samples()),
#    germlines=expand("out/{sample}.sorted.dups.bam", sample=germline_samples()),
#    bed=config["regions"],
#    reference=config["genome"]
#  output:
#    "out/cnvkit.reference.cnn"
#  shell:
#    "{config[module_R]} && "
#    "cnvkit.py batch {input.tumours} "
#      "--normal {input.germlines} "
#      "--targets {input.bed} "
#      "--annotate {config[refflat]} "
#      "--fasta {input.reference} "
#      "--access {config[cnvkit_access]} "
#      "--output-reference out/cnvkit.reference.cnn --output-dir out/ "
#      "--segment-method hmm "
#      "--diagram --scatter"

# this is an experimental in-house CNV caller
rule cnv_caller:
  input:
    bams=tumour_germline_bams,
    bed=config["regions"]

  output:
    "out/{tumour}.cnv.calls.tsv"

  log:
    "log/{tumour}.cnv.calls.log"

  shell:
    "tools/cnv_caller-{config[cnv_version]}/cnv_caller/call.py --verbose --tumour {input.bams[0]} --normal {input.bams[1]} --bed {input.bed} > {output}"

# cnv
rule cnv_summary:
  input:
    "out/{tumour}.cnv.calls.tsv"

  output:
    "out/{tumour}.cnv.tsv"

  log:
    "log/{tumour}.cnv.log"

  params:
    regions=' '.join(config["regions_of_interest"]),
    region_names=' '.join(config["genes_of_interest"]),
    region_padding=' '.join([str(x) for x in config["cnv_region_padding"]])

  shell:
    "tools/cnv_caller-{config[cnv_version]}/cnv_caller/group.py --plot out/{wildcards.tumour}.cnv --regions {params.regions} --region_names {params.region_names} --region_padding {params.region_padding} <{input} >{output} 2>{log}"

######################
# msi
rule msisensor_prep:
  input:
    reference=config["genome"]
  output:
    "out/msisensor.list"
  log:
    stderr="log/msisensor.list.log"
  shell:
    "tools/msisensor-{config[msisensor_version]}/binary/msisensor.linux scan -d {input.reference} -o {output}"

rule msisensor:
  input:
    microsatellites="out/msisensor.list",
    bed=config["regions"],
    bams=tumour_germline_bams
  output:
    "out/{tumour}.msisensor.tsv"
  log:
    stderr="log/{tumour}.msisensor.stderr"
  params:
    tumour="{tumour}"
  shell:
    "tools/msisensor-{config[msisensor_version]}/binary/msisensor.linux msi -d {input.microsatellites} -n {input.bams[1]} -t {input.bams[0]} -e {input.bed} -o tmp/{params.tumour}.msisensor && "
    "mv tmp/{params.tumour}.msisensor out/{params.tumour}.msisensor.tsv"

rule msisensor_combine:
  input:
    expand("out/{tumour}.msisensor.tsv", tumour=samples['tumours']),
  output:
    "out/aggregate/msisensor.tsv"
  shell:
    "src/combine_msisensor.py {input} > {output}"

## mutational signature context based filter
rule filter_context_artefact:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz",
    exposures="out/{tumour}.mutational_signature_v3.2_sbs.filter.exposures"
  output:
    vcf="out/{tumour}.intersect.pass.filter.signatures.vcf.gz",
    png="out/{tumour}.signatures_cosmic_v3.2_sbs.png"
  shell:
    "{config[module_htslib]} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/annotate_context.py --genome {input.reference} --vcf {input.vcf} > tmp/{wildcards.tumour}.intersect.pass.filter.context.vcf.gz && " # annotates
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/assign_signatures.py --definition tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3.2_sbs.txt --signatures {input.exposures} --artefacts tools/mutational_signature-{config[signature_version]}/data/signature_summary.tsv --vcf tmp/{wildcards.tumour}.intersect.pass.filter.context.vcf.gz --threshold 0.01 --plot {output.png} > {output.vcf}"

# burden exonic snvs
rule mutation_burden_artefact_filter:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.filter.signatures.vcf.gz", tumour=samples['tumours']),
    regions=config["regions"],
  output:
    "out/aggregate/mutation_rate.artefact_filter.tsv"
  log:
    stderr="log/mutation_rate.artefact_filter.stderr"
  shell:
    "src/mutation_rate.py --vcfs {input.vcfs} --bed {input.regions} --signature_artefacts tools/mutational_signature-{config[signature_version]}/data/signature_summary.tsv --signature_artefact_penalty 0.5 >{output} 2>{log.stderr}"

#rule mutation_rate_artefact:
#  input:
#    vcfs=expand("out/{tumour}.intersect.pass.filter.denoised.vcf.gz", tumour=samples['tumours']),
#    regions=config["regions"]
#  output:
#    "out/aggregate/mutation_rate.denoised.tsv"
#  log:
#    stderr="log/mutation_rate.denoised.stderr"
#  shell:
#    "src/mutation_rate.py --vcfs {input.vcfs} --bed {input.regions} >{output} 2>{log.stderr}"

##### mantis #####
rule mantis_prep:
  input:
    mantis="reference/mantis.bed",
    capture=config["regions"]
  output:
    "out/mantis.intersect.bed"
  shell:
    "{config[module_bedtools]} && "
    "bedtools intersect -a {input.mantis} -b {input.capture} > {output}"
    
rule mantis:
  input:
    genome=config["genome"],
    mantis="out/mantis.intersect.bed",
    bams=tumour_germline_bams,
    bais=tumour_germline_bais
  output:
    "out/{tumour}.mantis.tsv"
  params:
    tumour="{tumour}"
  shell:
    "([[ ! -e {input.bams[0]}.bai ]] && ln -sf $(pwd)/{input.bais[0]} {input.bams[0]}.bai || true) && "
    "([[ ! -e {input.bams[1]}.bai ]] && ln -sf $(pwd)/{input.bais[1]} {input.bams[1]}.bai || true) && "
    "{config[module_samtools]} && "
    "python tools/MANTIS-{config[mantis_version]}/mantis.py --tumor {input.bams[0]} --normal {input.bams[1]} --bedfile {input.mantis} --genome {input.genome} --output tmp/{params.tumour}.mantis -mrq 20.0 -mlq 25.0 -mlc 20 -mrr 1 && "
    "mv tmp/{params.tumour}.mantis {output}"

# take the Average row Difference column
rule mantis_combine:
  input:
    expand("out/{tumour}.mantis.tsv", tumour=samples['tumours'])
  output:
    "out/aggregate/mantis.tsv"
  shell:
    "src/combine_mantis.py {input} > {output}"

##### hla typing #####
# split fastqs because the memory usage of razers3 is ludicrous
#    "fastqsplitter -i {input.fastqs[0]} -o tmp/optitype.{wildcards.sample}/i1.1.fq.gz -o tmp/optitype.{wildcards.sample}/i1.2.fq.gz -o tmp/optitype.{wildcards.sample}/i1.3.fq.gz -o tmp/optitype.{wildcards.sample}/i1.4.fq.gz -o tmp/optitype.{wildcards.sample}/i1.5.fq.gz
#    "fastqsplitter -i {input.fastqs[1]} -o tmp/optitype.{wildcards.sample}/i2.1.fq.gz -o tmp/optitype.{wildcards.sample}/i2.2.fq.gz -o tmp/optitype.{wildcards.sample}/i1.3.fq.gz -o tmp/optitype.{wildcards.sample}/i1.4.fq.gz -o tmp/optitype.{wildcards.sample}/i1.5.fq.gz
rule optitype:
  input:
    fastqs=lambda wildcards: samples["samples"][wildcards.sample]
  output:
    result="out/{sample}.optitype.tsv",
    report="out/{sample}.optitype.pdf"
  shell:
    "{config[module_cplex]} && "
    "{config[module_gcc]} && "
    "mkdir -p tmp/optitype.{wildcards.sample} && "
    "tools/OptiType-1.3.5/OptiTypePipeline.py -c cfg/optitype.cfg -i {input.fastqs[0]} {input.fastqs[1]} --dna --verbose --outdir tmp/optitype.{wildcards.sample} --enumerate 10 && "
    "mv tmp/optitype.{wildcards.sample}/*/*.tsv {output.result} && "
    "mv tmp/optitype.{wildcards.sample}/*/*.pdf {output.report}"

rule optitype_combine:
  input:
    expand("out/{sample}.optitype.tsv", sample=samples['samples']),
  output:
    "out/aggregate/optitype.tsv"
  shell:
    "src/combine_optitype.py {input} > {output}"

# V=$W/AGRF_CAGRF18627_H7MG5DSXX/out/0350358019_T.intersect.pass.filter.vcf.gz
# B=$W/AGRF_CAGRF18627_H7MG5DSXX/out/0350358019_T.sorted.dups.bam
#gunzip < $V > v.vcf
#cp $V v.vcf
#neoepiscope swap -i v.vcf -o vs.vcf
#$S/HapCUT2-v.1.3.3/build/extractHAIRS --indels 1 --bam $B --VCF vs.vcf --out fragment_file
#$S/HapCUT2-v.1.3.3/build/HAPCUT2 --fragments fragment_file --VCF vs.vcf --output haplotype_output_file
#neoepiscope prep -v vs.vcf -c haplotype_output_file -o adjusted.vcf
#module load CUDA
#neoepiscope call -b hg19 -c adjusted.vcf -o result3.tsv -d ./neoepiscope.data --alleles 'A*02:01,A*02:01,B*56:01,B*15:01,C*03:04,C*01:02'
#rule neoepiscope:
#  input:
#    expand("out/{sample}.optitype.tsv", sample=samples['samples']),
#  output:
#    "out/aggregate/optitype.tsv"
#  shell:
#    "src/combine_optitype.py {input} > {output}"

##### msiseq style #####
rule msiseq:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']),
    repeats="reference/hg19repeats.msiseq.tsv",
    capture=config["regions"]
  output:
    result="out/aggregate/msiseq.tsv"
  shell:
    "src/msiseq.py --verbose --vcfs {input.vcfs} --repeats {input.repeats} --capture {input.capture} --threshold 0.1 > {output.result}"

##### cadd #####
rule cadd_prepare:
  input:
    regions=config["regions"],
    cadd="/data/scratch/projects/punim0567/peter/cadd_v1.6.vcf.gz" #config["cadd"]
  output:
    cadd="out/cadd.vcf.gz"
  shell:
    "{config[module_htslib]} && "
    "{config[module_bedtools]} && "
    "bedtools slop -b 100 -g reference/genome.lengths -i {input.regions} "
    "| sort -k1,1 -k2,2n "
    "| bedtools merge -i - > tmp/cadd_prepare.bed && "
    "bedtools intersect -a {input.cadd} -b tmp/cadd_prepare.bed -header "
    "| bgzip > {output.cadd} && "
    "tabix -p vcf {output.cadd}"

##### kataegis #####
#    "tools/kataegis-{config[kataegis_version]}/kataegis/annotate.py --plot_genome {output} --plot_prefix out/${tumour}.kataegis.zoomed. --just_kataegis < $f | bgzip > out/${t}.kataegis.vcf.gz"
# count 3 is more sensitive
rule kataegis:
  input:
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    png="out/{tumour}.kataegis.png",
    vcf="out/{tumour}.kataegis.vcf.gz"
  shell:
    "{config[module_htslib]} && "
    "tools/kataegis-{config[kataegis_version]}/kataegis/annotate.py --count 3 --plot_genome {output.png} --just_kataegis < {input.vcf} | bgzip > {output.vcf}"

# tumour purity 
rule purity:
  input:
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    "out/{tumour}.purity.tsv"
  params:
    tumour="{tumour}"
  shell:
    "python tools/purity-{config[purity_version]}/purity.py --tumour {params.tumour} < {input} | python tools/csvtools-{config[csvtools_version]}/csvtools/csvadd.py --delimiter '	' --name Sample --value {params.tumour} > {output}"

rule purity_combine:
  input:
    expand("out/{tumour}.purity.tsv", tumour=samples['tumours'])
  output:
    "out/aggregate/purity.tsv"
  shell:
    "python tools/csvtools-{config[csvtools_version]}/csvtools/csvmerge.py --delimiter '	' {input} > {output}"


# tumour purity - currently not working
rule purity_theta:
  input:
    reference=config["genome"],
    bed=config["regions"],
    bams=tumour_germline_bams
  output:
    "out/{tumour}.theta2.done"
  log:
    stderr="log/{tumour}.purity.stderr"
  params:
    tumour="{tumour}",
  shell:
    #"cnvkit.py batch {input.bams[0]} --normal {input.bams[1]} --targets {input.bed} --fasta {input.reference} --output-reference tmp/{params.tumour}.cnn --output-dir tmp/{params.tumour}.cnvkit --diagram && "
    #"cnvkit.py export theta tmp/{params.tumour}.cnvkit/{params.tumour}.sorted.dups.cns -r tmp/{params.tumour}.cnn -o tmp/{params.tumour}.cnvkit.interval_count && "
    "{config[module_python2]} && {config[module_samtools]} && "
    "mkdir -p ./tmp/{params.tumour}.theta2 && "
    "tools/theta2/bin/CreateExomeInput -s tmp/{params.tumour}.cnvkit.interval_count -t {input.bams[0]} -n {input.bams[1]} --FA {input.reference} --EXON_FILE tools/theta2/data/hg19.exons.bed --DIR ./tmp/{params.tumour}.theta2 --QUALITY 1 && "
    "touch {output}"

rule minibam:
  input:
    bam="out/{sample}.sorted.dups.bam"
  output:
    "out/{sample}.mini.bam"
  shell:
    "{config[module_samtools]} && "
    "samtools view -h -b {input.bam} {config[regions_of_interest_minibam]} > {output} && "
    "samtools index {output}"

# somatic
rule somatic_variant_summary:
  input:
    vcfs1=expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']),
    vcfs2=expand("out/{tumour}.strelka.somatic.indels.norm.annot.pass.vcf.gz", tumour=samples['tumours']),
    lohs=expand("out/{tumour}.loh.bed", tumour=samples['tumours'])

  output:
    tsv="out/aggregate/targetted_gene_summary.somatic.tsv",
    png="out/aggregate/targetted_gene_summary.somatic.png"

  shell:
    "src/targetted_gene_summary.py --vcfs {input.vcfs1} {input.vcfs2} --lohs {input.lohs} --genes {config[genes_of_interest]} --loci {config[regions_of_interest]} > {output.tsv} && "
    "src/targetted_gene_plot.py --target {output.png} < {output.tsv}"


# germline
rule germline_variant_summary:
  input:
    vcf="out/aggregate/germline_joint.hc.normalized.annot.revel.clinvar.cadd.vcf.gz"

  output:
    tsv="out/aggregate/targetted_gene_summary.germline.tsv",
    png="out/aggregate/targetted_gene_summary.germline.png"

  shell:
    "src/targetted_gene_summary.py --vcfs {input.vcf} --genes {config[genes_of_interest]} --loci {config[regions_of_interest]} --multisample > {output.tsv} && "
    "src/targetted_gene_plot.py --target {output.png} < {output.tsv}"

##### msmutect
rule msmutect_somatic:
  input:
    regions_msmutect=config["regions_msmutect"],
    bams=tumour_germline_dup_bams
  output:
    tsv="out/msmutect/{tumour}.msmutect.tsv"
  params:
    cores=cluster["msmutect_somatic"]["n"]
  shell:
    "msmutect -T {input.bams[0]} -N {input.bams[1]} -l {input.regions_msmutect} -O tmp/{wildcards.tumour} -c {params.cores} -m -A && "
    "mv tmp/{wildcards.tumour}.full.mut.tsv {output}"

rule combine_msmutect_somatic:
  input:
    expand("out/msmutect/{tumour}.msmutect.tsv", tumour=samples['tumours'])
  output:
    "out/aggregate/msmutect.combined.tsv"
  shell:
    "src/combine_msmutect.py {input} > {output}"


##### ethnicity with somalier
rule somalier_extract:
  input:
    reference=config["genome"],
    bam="out/{sample}.sorted.dups.bam",
  output:
    "out/{sample}.somalier"
  shell:
    "tools/somalier-v{config[somalier_version]}/somalier extract -d out --sites tools/somalier-v{config[somalier_version]}/sites.GRCh37.vcf.gz -f {input.reference} {input.bam}"

rule somalier:
  input:
    samples=expand("out/{sample}.somalier", sample=samples['samples'])
  output:
    "out/aggregate/somalier.done"
  shell:
    "tools/somalier-v{config[somalier_version]}/somalier relate --output-prefix=out/aggregate/somalier-relate {input.samples} && "
    "tools/somalier-v{config[somalier_version]}/somalier ancestry --output-prefix=out/aggregate/somalier-ancestry --labels tools/somalier-v{config[somalier_version]}/ancestry-labels-1kg.tsv tools/somalier-v{config[somalier_version]}/1kg-somalier/*.somalier ++ {input.samples} && "
    "touch out/aggregate/somalier.done"

