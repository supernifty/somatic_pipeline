print("starting")
VERSION="0.3"

configfile: "cfg/config.yaml"

import yaml
cluster = yaml.load(open("cfg/cluster.yaml"))
samples = yaml.load(open("cfg/samples.yaml"))
print("loaded configuration")


# NOTE: no mitochondria MT because they aren't in our exome
GATK_CHROMOSOMES=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')

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

def germline_samples():
  all_samples = set(samples['samples'])
  tumours = set(samples['tumours'])
  return list(all_samples.difference(tumours))

def germline_sample(wildcards):
  return samples["tumours"][wildcards.tumour]

### final outputs ###
rule all:
  input:
    expand("out/{tumour}.strelka.somatic.snvs.af.norm.vep.pass.vcf.gz", tumour=samples['tumours']), # somatic snvs strelka
    expand("out/{tumour}.strelka.somatic.indels.norm.vep.pass.vcf.gz", tumour=samples['tumours']), # somatic indels strelka
    expand("out/{tumour}.strelka.somatic.indels.norm.vep.pass.af.filter.vcf.gz", tumour=samples['tumours']), # somatic indels strelka
    expand("out/{germline}.strelka.germline.filter_gt.vep.vcf.gz", germline=germline_samples()), # germline strelka calls
    # -- disabled for now expand("out/{germline}.strelka.germline.filter_gt.vep.vcf.gz", germline=samples['tumours']), # run the tumours as if they were germline

    expand("out/{sample}.oxo_metrics.txt", sample=samples['samples']),
    expand("out/{sample}.artifact_metrics.txt.error_summary_metrics", sample=samples['samples']),
    # -- disabled too slow expand("out/{tumour}.strelka.somatic.snvs.bias.vcf.gz", tumour=samples['tumours']),
    # -- disabled expand("out/{tumour}.mutect2.filter.bias.vcf.gz", tumour=samples['tumours']), # somatic mutect2 with dkfz bias annotation

    expand("out/{tumour}.mutect2.filter.norm.vep.pass.vcf.gz", tumour=samples['tumours']), # somatic mutect2
    expand("out/fastqc/{sample}/completed", sample=samples['samples']), # fastqc
    expand("out/{sample}.metrics.insertsize", sample=samples['samples']),
    expand("out/{sample}.metrics.alignment", sample=samples['samples']),
    expand("out/{sample}.metrics.target", sample=samples['samples']),

    # -- issues expand("out/{tumour}.verifybamid.somatic.completed", tumour=samples['tumours']),
    expand("out/{tumour}.mutect2.filter.norm.annot.vcf.gz", tumour=samples['tumours']),
    expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']), # somatic combination of strelka and mutect2, with filter
    expand("out/{tumour}.pass_one.vcf.gz", tumour=samples['tumours']), # somatic combination of strelka and mutect2, with filter

    # -- expand("out/{tumour}.varscan.copynumber.deletions.bed", tumour=samples['tumours']), # TODO fails on mini sample
    expand("out/{tumour}.platypus.somatic.vcf.gz", tumour=samples['tumours']), # platypus somatic calls
    expand("out/{tumour}.loh.bed", tumour=samples['tumours']), # loh regions
    expand("out/{tumour}.cnv.tsv", tumour=samples['tumours']), # cnv regions
    expand("out/{tumour}.strelka.somatic.af.png", tumour=samples['tumours']), # plot strelka af
    expand("out/{tumour}.strelka.somatic.pass.signatures.af.png", tumour=samples['tumours']), # plot strelka af
    expand("out/{tumour}.mutect2.somatic.af.png", tumour=samples['tumours']), # plot mutect2 af
    expand("out/{tumour}.intersect.pass.filter.signatures.vcf.gz", tumour=samples['tumours']), # annotated outputs with signatures
    expand("out/{tumour}.strelka.somatic.mutational_signature_v3_sbs.png", tumour=samples['tumours']), # signature profile with signature allocation

    expand("out/{sample}.mini.bam", sample=samples['samples']), 

    # tumour purity
    # -- not working expand("out/{tumour}.theta2.done", tumour=samples['tumours']),

    # msi
    "out/aggregate/msisensor.tsv",

    # combined results
    "out/aggregate/mutect2.filter.genes_of_interest.combined.tsv",
    "out/aggregate/mutect2.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v2.combined.tsv",
    "out/aggregate/mutational_signatures_v2.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3_sbs.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3_dbs.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3_sbs_strand.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3_sbs_strelka_strand.filter.combined.tsv",
    "out/aggregate/mutational_signatures_v3_dbs.combined.tsv",
    "out/aggregate/mutational_signatures_v3_id.combined.tsv",
    "out/aggregate/mutational_signatures_v3_id_strelka.filter.combined.tsv",

    "out/aggregate/max_coverage.tsv",
    "out/aggregate/max_trimmed_coverage.tsv",
    "out/aggregate/ontarget.tsv",
    "out/aggregate/germline_joint.hc.normalized.annot.tsv", # gatk calls for all germline samples
    # -- disabled "out/aggregate/tumour_joint.hc.normalized.annot.vcf.gz", # gatk calls for all tumour samples (as germline)
    "out/aggregate/qc.summary.tsv",
    "out/aggregate/multiqc.html", # overall general qc
    "out/aggregate/ontarget.png", # combined ontarget coverage plots
    "out/aggregate/ontarget_tumour.png", # somatic ontarget coverage plots
    "out/aggregate/ontarget_germline.png", # germline ontarget coverage plots
    "out/aggregate/ontarget_tumour.dedup.png",
    "out/aggregate/mutation_rate.tsv",
    "out/aggregate/mutation_rate.artefact_filter.tsv",
    "out/aggregate/mutation_rate.high_af.tsv",
    "out/aggregate/mutation_rate.low_af.tsv",
    "out/aggregate/msi_burden.tsv",
    "out/aggregate/targetted_gene_summary.somatic.tsv",
    "out/aggregate/targetted_gene_summary.germline.tsv",
    "out/aggregate/gene_coverage.tsv",
    "out/aggregate/gene_coverage_min.tumour.tsv",
    "out/aggregate/final.html",

### aggregate ###

# write out all tool versions (TODO)
rule make_versions:
  output:
    versions="out/aggregate/versions.txt"
  shell:
    "src/make_tsv.py --columns Tool Version --rows "
    "Pipeline,{VERSION} "
    ">{output.versions}"

rule report_md:
  input:
    versions="out/aggregate/versions.txt",
    signatures="out/aggregate/mutational_signatures_v3_sbs.filter.combined.tsv",
    burden="out/aggregate/mutation_rate.tsv",
    qc="out/aggregate/qc.summary.tsv",
    selected_somatic_variants="out/aggregate/mutect2.filter.genes_of_interest.combined.tsv",
    all_somatic_variants="out/aggregate/mutect2.filter.combined.tsv",
    all_germline_variants="out/aggregate/germline_joint.hc.normalized.annot.tsv"
  
  output:
    md="out/aggregate/final.md",
    html="out/aggregate/final.html"

  log:
    stderr="log/make_report.stderr"

  shell:
    "src/make_report.py --versions {input.versions} --signatures {input.signatures} --burden {input.burden} --qc {input.qc} --selected_somatic_variants {input.selected_somatic_variants} --all_somatic_variants {input.all_somatic_variants} --all_germline_variants {input.all_germline_variants} > {output.md} 2>{log.stderr} && "
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
    "{config[module_bedtools]} && "
    "samtools view -uF 0x400 {input.bam} > tmp/{wildcards.sample}.dedup.bam && "
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

rule qc_papillon_min:
  input:
    bams=expand("out/{sample}.sorted.dups.bam", sample=samples['samples']),
    tumours=expand("out/{tumour}.sorted.dups.bam", tumour=samples['tumours']),
    germlines=expand("out/{germline}.sorted.dups.bam", germline=germline_samples())
  output:
    tumour="out/aggregate/gene_coverage_min.tumour.tsv",
    germline="out/aggregate/gene_coverage_min.germline.tsv"
  log:
    stderr="log/papillon_min.stderr"
  params:
    gene_list=' '.join(config["genes_of_interest"])
  shell:
    "(tools/papillon-{config[papillon_version]}/papillon/papillon.py --max_coverage {config[papillon_max_coverage_tumour]} --stat min --bed {config[refseq_bed]} --bams {input.tumours} --genes {params.gene_list} --plot out/aggregate/gene_coverage.tumour.min --padding {config[papillon_padding]} --exon_plots >{output.tumour} && "
    "tools/papillon-{config[papillon_version]}/papillon/papillon.py --max_coverage {config[papillon_max_coverage_germline]} --stat min --bed {config[refseq_bed]} --bams {input.germlines} --genes {params.gene_list} --plot out/aggregate/gene_coverage.germline.min --padding {config[papillon_padding]} --exon_plots >{output.germline}) 2>{log.stderr}"

rule qc_papillon:
  input:
    bams=expand("out/{sample}.sorted.dups.bam", sample=samples['samples']),
    tumours=expand("out/{tumour}.sorted.dups.bam", tumour=samples['tumours']),
    germlines=expand("out/{germline}.sorted.dups.bam", germline=germline_samples())
  output:
    all="out/aggregate/gene_coverage.tsv",
    tumour="out/aggregate/gene_coverage.tumour.tsv",
    germline="out/aggregate/gene_coverage.germline.tsv"
  log:
    stderr="log/papillon.stderr"
  params:
    gene_list=' '.join(config["genes_of_interest"])
  shell:
    "(tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat median --bed {config[refseq_bed]} --bams {input.bams} --genes {params.gene_list} --plot out/aggregate/gene_coverage.median --exon_plots >{output.all} && "
    "tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat median --bed {config[refseq_bed]} --bams {input.tumours} --genes {params.gene_list} --plot out/aggregate/gene_coverage.tumour.median >{output.tumour} && "
    "tools/papillon-{config[papillon_version]}/papillon/papillon.py --stat median --bed {config[refseq_bed]} --bams {input.germlines} --genes {params.gene_list} --plot out/aggregate/gene_coverage.germline.median >{output.germline}) 2>{log.stderr}"

### alignment ###
rule trim:
  input:
    fastqs=lambda wildcards: samples["samples"][wildcards.sample]
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
    #fastqs=lambda wildcards: samples["samples"][wildcards.sample]

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
    gvcfs=expand("out/{germline}.hc.gvcf.gz", germline=samples["tumours"]),
    reference=config["genome"]
  output:
    "out/tumour_joint_{chromosome}.vcf"
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
    expand("out/{sample}.max_coverage", sample=samples['samples']),
  output:
    "out/aggregate/max_coverage.tsv"
  shell:
    ">{output} && "
    "for f in {input}; do echo \"$f $(grep \"Max coverage\" $f)\" >> {output}; done"

rule qc_max_trimmed_coverage_combine:
  input:
    expand("out/{sample}.max_trimmed_coverage", sample=samples['samples']),
  output:
    "out/aggregate/max_trimmed_coverage.tsv"
  shell:
    ">{output} && "
    "for f in {input}; do echo \"$f $(grep \"Max coverage\" $f)\" >> {output}; done"

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
rule mutect2_sample_pon:
  input:
    reference=config["genome"],
    bam="out/{germline}.sorted.dups.bam",
    regions=config["regions"]
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
    regions=config["regions"],
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
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

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
    "src/annotate.sh tmp/{wildcards.tumour}.strelka.somatic.indels.norm.vcf.gz {output} {input.reference} {params.cores}) 2>{log}"

rule annotate_vep_hc:
  input:
    vcf="out/germline_joint.hc.normalized.vcf",
    reference=config['genome']
  output:
    "tmp/germline_joint.hc.normalized.vep.vcf.gz"
  log:
    "log/hc.vep.log"
  params:
    cores=cluster["annotate_vep_germline"]["n"]
  shell:
    "{config[module_samtools]} && "
    "{config[module_htslib]} && "
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

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
    "tools/vt-0.577/vt decompose -s {input.vcf} | tools/vt-0.577/vt normalize -n -r {input.reference} - -o out/{wildcards.tumour}.mutect2.filter.norm.vcf.gz && "
    "src/annotate.sh out/{wildcards.tumour}.mutect2.filter.norm.vcf.gz {output} {input.reference} {params.cores} 2>{log}"

rule annotate_clinvar_mutect2:
  input:
    vcf="out/{tumour}.mutect2.filter.norm.vep.vcf.gz"
  output:
    vcf="out/{tumour}.mutect2.filter.norm.annot.vcf.gz"
  log:
    "log/{tumour}.clinvar.log"
  shell:
    "{config[module_htslib]} && "
    "tools/vcfanno_linux64 -lua cfg/vcfanno.lua cfg/vcfanno.cfg {input.vcf} | bgzip > {output.vcf} 2>{log}"

rule annotate_clinvar_strelka_snvs:
  input:
    vcf="out/{tumour}.strelka.somatic.snvs.af.norm.vep.vcf.gz"
  output:
    vcf="out/{tumour}.strelka.somatic.snvs.af.norm.annot.vcf.gz"
  log:
    "log/{tumour}.clinvar.log"
  shell:
    "{config[module_htslib]} && "
    "tools/vcfanno_linux64 -lua cfg/vcfanno.lua cfg/vcfanno.cfg {input.vcf} | bgzip > {output.vcf} 2>{log}"

rule annotate_clinvar_strelka_indels:
  input:
    vcf="out/{tumour}.strelka.somatic.indels.norm.vep.vcf.gz"
  output:
    vcf="out/{tumour}.strelka.somatic.indels.norm.annot.vcf.gz"
  log:
    "log/{tumour}.clinvar.log"
  shell:
    "{config[module_htslib]} && "
    "tools/vcfanno_linux64 -lua cfg/vcfanno.lua cfg/vcfanno.cfg {input.vcf} | bgzip > {output.vcf} 2>{log}"

rule annotate_clinvar_hc:
  input:
    vcf="tmp/germline_joint.hc.normalized.vep.vcf.gz"
  output:
    vcf="out/aggregate/germline_joint.hc.normalized.annot.vcf.gz"
  log:
    "log/hc.clinvar.log"
  shell:
    "{config[module_htslib]} && "
    "tools/vcfanno_linux64 -lua cfg/vcfanno.lua cfg/vcfanno.cfg {input.vcf} | bgzip > {output.vcf} 2>{log}"

rule pass_mutect2:
  input:
    "out/{tumour}.mutect2.filter.norm.annot.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.norm.annot.pass.vcf.gz"
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input} | egrep '(^#|PASS)' | bgzip > {output}"

# TODO LowDepth fix
rule pass_strelka_indels:
  input:
    "out/{tumour}.strelka.somatic.indels.norm.annot.vcf.gz",
  output:
    "out/{tumour}.strelka.somatic.indels.norm.annot.pass.vcf.gz",
  shell:
    "{config[module_htslib]} && "
    "gunzip < {input} | egrep '(^#|PASS|	LowDepth	)' | bgzip > {output}"

# TODO LowDepth fix
rule pass_strelka_snvs:
  input:
    "out/{tumour}.strelka.somatic.snvs.af.norm.annot.vcf.gz",
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
    "src/annotate.sh {input.vcf} {output} {input.reference} {params.cores} 2>{log}"

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
    "tools/vt-0.577/vt decompose -s {input.strelka_snvs} | tools/vt-0.577/vt normalize -n -r {input.reference} - -o {output.snvs_norm} && "
    "tools/vt-0.577/vt decompose -s {input.strelka_indels} | tools/vt-0.577/vt normalize -n -r {input.reference} - -o {output.indels_norm}"
 
rule intersect_somatic_callers:
  input:
    reference=config["genome"],
    mutect2="out/{tumour}.mutect2.filter.norm.annot.vcf.gz",
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.norm.vcf.gz" 
  output:
    "out/{tumour}.intersect.vcf.gz"
  log:
    stderr="log/{tumour}.intersect.log"
  shell:
    "({config[module_samtools]} && "
    "{config[module_htslib]} && "
    "{config[module_bedtools]} && "
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_snvs} {input.mutect2} > tmp/{wildcards.tumour}.intersect.unsorted.vcf && "
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_indels} {input.mutect2} | sed -n '/^#/!p' >> tmp/{wildcards.tumour}.intersect.unsorted.vcf && "
    "grep '^#' tmp/{wildcards.tumour}.intersect.unsorted.vcf > tmp/{wildcards.tumour}.intersect.vcf && "
    "bedtools sort -faidx reference/genome.lengths < tmp/{wildcards.tumour}.intersect.unsorted.vcf >> tmp/{wildcards.tumour}.intersect.vcf && "
    "bgzip < tmp/{wildcards.tumour}.intersect.vcf > {output}"
    ") 2>{log.stderr}"
 
rule pass_one_somatic_callers:
  input:
    reference=config["genome"],
    mutect2="out/{tumour}.mutect2.filter.norm.annot.vcf.gz",
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.norm.vcf.gz" 
  output:
    "out/{tumour}.pass_one.vcf.gz"
  log:
    stderr="log/{tumour}.pass_one.log"
  shell:
    "({config[module_samtools]} && "
    "{config[module_htslib]} && "
    "{config[module_bedtools]} && "
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_snvs} {input.mutect2} --pass_one > tmp/{wildcards.tumour}.pass_one.unsorted.vcf && " # snvs vs mutect
    "src/vcf_intersect.py --allowed_filters str_contraction LowDepth --inputs {input.strelka_indels} {input.mutect2} --pass_one | sed -n '/^#/!p' >> tmp/{wildcards.tumour}.pass_one.unsorted.vcf && " # indels vs mutect
    "grep '^#' tmp/{wildcards.tumour}.pass_one.unsorted.vcf > tmp/{wildcards.tumour}.pass_one.vcf && " # header
    "bedtools sort -faidx reference/genome.lengths < tmp/{wildcards.tumour}.pass_one.unsorted.vcf >> tmp/{wildcards.tumour}.pass_one.vcf && "
    "bgzip < tmp/{wildcards.tumour}.pass_one.vcf > {output}"
    ") 2>{log.stderr}"
 
rule intersect_pass_somatic_callers:
  input:
    reference=config["genome"],
    mutect2="out/{tumour}.mutect2.filter.norm.vep.pass.vcf.gz",
    strelka_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.vep.pass.vcf.gz",
    strelka_indels="out/{tumour}.strelka.somatic.indels.norm.vep.pass.vcf.gz" 
  output:
    "out/{tumour}.intersect.pass.vcf.gz"
  log:
    stderr="log/{tumour}.intersect.log"
  shell:
    "({config[module_samtools]} && "
    "{config[module_htslib]} && "
    "{config[module_bedtools]} && "
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
    "out/{tumour}.strelka.somatic.indels.norm.vep.pass.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.indels.norm.vep.pass.af.vcf.gz"
  log:
    stderr="log/{tumour}.strelka.indels.annotate_af.stderr"
  shell:
    "src/annotate_indel_af.py --verbose --sample TUMOR < {input} > {output} 2>{log.stderr}"

rule filter_strelka_snvs:
  input:
    "out/{tumour}.strelka.somatic.snvs.af.norm.vep.pass.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.snvs.af.norm.vep.pass.filter.vcf.gz"
  log:
    stderr="log/{tumour}.strelka.snvs.filter.stderr"
  shell:
    "src/filter_af.py --verbose --sample TUMOR --info_af --af {config[af_threshold]} --dp {config[dp_threshold]} < {input} > {output} 2>{log.stderr}"

rule filter_indels:
  input:
    "out/{tumour}.strelka.somatic.indels.norm.vep.pass.af.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.indels.norm.vep.pass.af.filter.vcf.gz"
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

rule combine_mutect2_tsv:
  input:
    expand("out/{tumour}.mutect2.filter.vep.tsv", tumour=samples['tumours'])
  output:
    "out/aggregate/mutect2.filter.combined.tsv"
  shell:
    "src/combine_tsv_raw.py {input} | sed 's/^out\\/\\([^.]*\\)\\.[^\\t]*/\\1/' > {output}"

rule mutect2_tsv:
  input:
    vcf="out/{tumour}.mutect2.filter.norm.annot.vcf.gz"
  output:
    "out/{tumour}.mutect2.filter.vep.tsv"
  shell:
    "src/vcf2tsv.py {input.vcf} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK' >{output}"

rule germline_tsv:
  input:
    vcf="out/aggregate/germline_joint.hc.normalized.annot.vcf.gz"
  output:
    "out/aggregate/germline_joint.hc.normalized.annot.tsv"
  shell:
    "src/vcf2tsv.py {input.vcf} | "
    "src/extract_vep.py --header 'Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK' | python tools/csvtools-{config[csvtools_version]}/csvtools/csvfilter.py --delimiter '	' --filters 'GT!0/0' 'GT!./.' >{output}"

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
    vcf="out/{tumour}.mutect2.filter.norm.annot.vcf.gz"
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

# ----- mutational signatures
rule mutational_signature_v2:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.vcf.gz"
  output:
    "out/{tumour}.mutational_signature_v2.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v2.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} > out/{wildcards.tumour}.mutational_signature_v2.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v2.txt --counts out/{wildcards.tumour}.mutational_signature_v2.counts > {output}) 2>{log.stderr}"
#    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/plot_counts.py out/{wildcards.tumour}.mutational_signature_v2.png < out/{wildcards.tumour}.mutational_signature_v2.counts && "

rule combine_mutational_signatures_v2:
  input:
    expand("out/{tumour}.mutational_signature_v2.exposures", tumour=samples['tumours']),
  output:
    "out/aggregate/mutational_signatures_v2.combined.tsv"
  shell:
    "src/combine_tsv.py {input} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v2\\.exposures/\\1/' > {output}"

# mutational signatures with filtered counts
rule mutational_signature_filtered_v2:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    "out/{tumour}.mutational_signature_v2.filter.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v2.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} > out/{wildcards.tumour}.mutational_signature_v2.filter.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v2.txt --counts out/{wildcards.tumour}.mutational_signature_v2.filter.counts > {output}) 2>{log.stderr}"
#    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/plot_counts.py out/{wildcards.tumour}.mutational_signature_v2.filter.png < out/{wildcards.tumour}.mutational_signature_v2.filter.counts && "

rule combine_mutational_signatures_filtered_v2:
  input:
    expand("out/{tumour}.mutational_signature_v2.filter.exposures", tumour=samples['tumours']),
  output:
    "out/aggregate/mutational_signatures_v2.filter.combined.tsv"
  shell:
    "src/combine_tsv.py {input} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v2\\.filter\\.exposures/\\1/' >{output}"


# disabled for now   "python tools/mutational_signature-{config[signature_version]}/mutational_signature/plot_components.py --threshold 0.05 --show_signature --target {output}.png --order 'Signature.1' 'Signature.2' 'Signature.3' 'Signature.4' 'Signature.5' 'Signature.6' 'Signature.7' 'Signature.8' 'Signature.9' 'Signature.10' 'Signature.11' 'Signature.12' 'Signature.13' 'Signature.14' 'Signature.15' 'Signature.16' 'Signature.17' 'Signature.18' 'Signature.19' 'Signature.20' 'Signature.21' 'Signature.22' 'Signature.23' 'Signature.24' 'Signature.25' 'Signature.26' 'Signature.27' 'Signature.28' 'Signature.29' 'Signature.30' --descriptions '5-methylcytosine deamination' 'APOBEC' 'double-strand break-repair failure' 'tobacco mutagens' '' 'defective mismatch repair' 'ultraviolet light exposure' '' '' 'POLE mutations' 'alkylating agents' '' 'APOBEC' '' 'defective mismatch repair' '' '' '' '' 'defective mismatch repair' '' '' '' 'aflatoxin exposure' '' 'defective mismatch repair' '' '' 'tobacco' '' < {output}"

# mutational signatures with filtered counts
rule mutational_signature_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.vcf.gz"
  output:
    sbs="out/{tumour}.mutational_signature_v3_sbs.exposures",
    dbs="out/{tumour}.mutational_signature_v3_dbs.exposures",
    id="out/{tumour}.mutational_signature_v3_id.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3.stderr" # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --doublets --genome {input.reference} --vcf {input.vcf} > out/{wildcards.tumour}.mutational_signature_v3.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_sbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3.counts > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_id.txt --counts out/{wildcards.tumour}.mutational_signature_v3.counts > {output.id} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_dbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3.counts > {output.dbs}) 2>{log.stderr}"

# mutational signatures with filtered counts
rule mutational_signature_filtered_v3:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    sbs="out/{tumour}.mutational_signature_v3_sbs.filter.exposures",
    dbs="out/{tumour}.mutational_signature_v3_dbs.filter.exposures",
    id="out/{tumour}.mutational_signature_v3_id.filter.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3_filter.stderr"
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --doublets --genome {input.reference} --vcf {input.vcf} > out/{wildcards.tumour}.mutational_signature_v3.filter.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_sbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3.filter.counts > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_id.txt --counts out/{wildcards.tumour}.mutational_signature_v3.filter.counts > {output.id} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_dbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3.filter.counts > {output.dbs}) 2>{log.stderr}"

# mutational signatures with filtered counts
rule mutational_signature_filtered_v3_id:
  input:
    reference=config["genome"],
    vcf_indels="out/{tumour}.strelka.somatic.indels.norm.vep.pass.af.filter.vcf.gz",
    vcf_snvs="out/{tumour}.strelka.somatic.snvs.af.norm.vep.pass.filter.vcf.gz"
  output:
    sbs="out/{tumour}.mutational_signature_v3_sbs_strelka.filter.exposures",
    id="out/{tumour}.mutational_signature_v3_id_strelka.filter.exposures"
  log:
    stderr="out/{tumour}.mutational_signature_v3_strelka.stderr", # keep for now
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --indels --just_indels --genome {input.reference} --vcf {input.vcf_indels} > out/{wildcards.tumour}.mutational_signature_v3_strelka_indels.filter.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf_snvs} > out/{wildcards.tumour}.mutational_signature_v3_strelka_snvs.filter.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_sbs.txt --counts out/{wildcards.tumour}.mutational_signature_v3_strelka_snvs.filter.counts > {output.sbs} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_id.txt --counts out/{wildcards.tumour}.mutational_signature_v3_strelka_indels.filter.counts > {output.id}) 2>{log.stderr}"

# strand bias included (TODO some redundancy in calculation)
rule mutational_signature_filtered_v3_sbs_strelka_strand:
  input:
    reference=config["genome"],
    transcripts=config["transcripts"],
    vcf="out/{tumour}.strelka.somatic.snvs.af.norm.vep.pass.filter.vcf.gz"
  output:
    sbs="out/{tumour}.mutational_signature_v3_sbs_strelka_strand.filter.exposures",
  log:
    stderr="out/{tumour}.mutational_signature_v3_sbs_strelka_strand_filter.stderr"
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} --transcripts {input.transcripts} > out/{wildcards.tumour}.mutational_signature_v3_sbs_strelka_strand.filter.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_sbs_tx.txt --counts out/{wildcards.tumour}.mutational_signature_v3_sbs_strelka_strand.filter.counts --strand > {output.sbs}) 2>{log.stderr}"

# strand bias included (TODO some redundancy in calculation)
rule mutational_signature_filtered_v3_sbs_strand:
  input:
    reference=config["genome"],
    transcripts=config["transcripts"],
    vcf="out/{tumour}.intersect.pass.filter.vcf.gz"
  output:
    sbs="out/{tumour}.mutational_signature_v3_sbs_strand.filter.exposures",
  log:
    stderr="out/{tumour}.mutational_signature_v3_sbs_strand_filter.stderr"
  shell:
    "(python tools/mutational_signature-{config[signature_version]}/mutational_signature/count.py --genome {input.reference} --vcf {input.vcf} --transcripts {input.transcripts} > out/{wildcards.tumour}.mutational_signature_v3_sbs_strand.filter.counts && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/decompose.py --signatures tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_sbs_tx.txt --counts out/{wildcards.tumour}.mutational_signature_v3_sbs_strand.filter.counts --strand > {output.sbs}) 2>{log.stderr}"

# signatures from all samples
rule combine_mutational_signatures_filtered_v3:
  input:
    sbs=expand("out/{tumour}.mutational_signature_v3_sbs.filter.exposures", tumour=samples['tumours']),
    sbs2=expand("out/{tumour}.mutational_signature_v3_sbs_strelka.filter.exposures", tumour=samples['tumours']),
    sbs_strand=expand("out/{tumour}.mutational_signature_v3_sbs_strand.filter.exposures", tumour=samples['tumours']),
    sbs2_strand=expand("out/{tumour}.mutational_signature_v3_sbs_strelka_strand.filter.exposures", tumour=samples['tumours']),
    dbs=expand("out/{tumour}.mutational_signature_v3_dbs.filter.exposures", tumour=samples['tumours']),
    dbs2=expand("out/{tumour}.mutational_signature_v3_dbs.exposures", tumour=samples['tumours']),
    id=expand("out/{tumour}.mutational_signature_v3_id.exposures", tumour=samples['tumours']),
    id2=expand("out/{tumour}.mutational_signature_v3_id_strelka.filter.exposures", tumour=samples['tumours'])
  output:
    sbs="out/aggregate/mutational_signatures_v3_sbs.filter.combined.tsv",
    sbs2="out/aggregate/mutational_signatures_v3_sbs_strelka.filter.combined.tsv",
    sbs_strand="out/aggregate/mutational_signatures_v3_sbs_strand.filter.combined.tsv",
    sbs2_strand="out/aggregate/mutational_signatures_v3_sbs_strelka_strand.filter.combined.tsv",
    dbs="out/aggregate/mutational_signatures_v3_dbs.filter.combined.tsv",
    dbs2="out/aggregate/mutational_signatures_v3_dbs.combined.tsv",
    id="out/aggregate/mutational_signatures_v3_id.combined.tsv",
    id2="out/aggregate/mutational_signatures_v3_id_strelka.filter.combined.tsv"
  shell:
    "src/combine_tsv.py {input.sbs} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_sbs\\.filter\\.exposures/\\1/' >{output.sbs} && "
    "src/combine_tsv.py {input.sbs2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_sbs_strelka\\.filter\\.exposures/\\1/' >{output.sbs2} && "
    "src/combine_tsv.py {input.sbs_strand} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_sbs_strand\\.filter\\.exposures/\\1/' >{output.sbs_strand} && "
    "src/combine_tsv.py {input.sbs2_strand} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_sbs_strelka_strand\\.filter\\.exposures/\\1/' >{output.sbs2_strand} && "
    "src/combine_tsv.py {input.dbs} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_dbs\\.filter\\.exposures/\\1/' >{output.dbs} && "
    "src/combine_tsv.py {input.dbs2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_dbs\\.exposures/\\1/' >{output.dbs2} && "
    "src/combine_tsv.py {input.id} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_id\\.exposures/\\1/' >{output.id} && "
    "src/combine_tsv.py {input.id2} | sed 's/^out\\/\\(.*\)\\.mutational_signature_v3_id_strelka\\.filter\\.exposures/\\1/' >{output.id2}"


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

# burden exonic snvs
rule mutation_burden_low_af:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.vcf.gz", tumour=samples['tumours']),
    regions=config["regions"]
  output:
    "out/aggregate/mutation_rate.low_af.tsv"
  log:
    stderr="log/mutation_rate.stderr"
  shell:
    "src/mutation_rate.py --verbose --vcfs {input.vcfs} --bed {input.regions} --min_af 0.01 >{output} 2>{log.stderr}"

# burden exonic snvs
rule mutation_burden_high_af:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']),
    regions=config["regions"],
  output:
    "out/aggregate/mutation_rate.high_af.tsv"
  log:
    stderr="log/mutation_rate.stderr"
  shell:
    "src/mutation_rate.py --verbose --vcfs {input.vcfs} --bed {input.regions} --min_af 0.2 >{output} 2>{log.stderr}"

# burden exonic snvs
rule mutation_burden:
  input:
    vcfs=expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']),
    regions=config["regions"],
  output:
    "out/aggregate/mutation_rate.tsv"
  log:
    stderr="log/mutation_rate.stderr"
  shell:
    "src/mutation_rate.py --verbose --vcfs {input.vcfs} --bed {input.regions} >{output} 2>{log.stderr}"

# burden exonic indels 
rule msi_burden:
  input:
    vcfs=expand("out/{tumour}.strelka.somatic.indels.norm.vep.pass.vcf.gz", tumour=samples['tumours']),
    regions="reference/msi.regions.bed"
  output:
    "out/aggregate/msi_burden.tsv"
  log:
    stderr="log/msi_burden.stderr"
  shell:
    "src/mutation_rate.py --verbose --vcfs {input.vcfs} --bed {input.regions} --indels_only --sample_name TUMOR >{output} 2>{log.stderr}"

# af distribution
rule plot_af_strelka:
  input:
    "out/{tumour}.strelka.somatic.snvs.af.vcf.gz"
  output:
    "out/{tumour}.strelka.somatic.af.png"
  shell:
    "src/plot_af.py --log --sample TUMOR --target {output} --dp {config[dp_threshold]} --info_af --title 'Variant count as a function of somatic allele fraction for {wildcards.tumour}' < {input}"

# af distribution

rule annotate_strelka_signatures:
  input:
    reference=config["genome"],
    vcf="out/{tumour}.strelka.somatic.snvs.af.vcf.gz",
    exposures="out/{tumour}.mutational_signature_v3_sbs.filter.exposures"
  output:
    vcf="out/{tumour}.strelka.somatic.snvs.af.signatures.vcf.gz",
    plot="out/{tumour}.strelka.somatic.mutational_signature_v3_sbs.png"
  shell:
    "{config[module_htslib]} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/annotate_context.py --genome {input.reference} --vcf {input.vcf} > out/{wildcards.tumour}.somatic.snvs.af.contexts.vcf.gz && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/assign_signatures.py --definition tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_sbs.txt --signatures {input.exposures} --vcf out/{wildcards.tumour}.somatic.snvs.af.contexts.vcf.gz --threshold 0.01 --plot {output.plot} > {output.vcf}"
  
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
    "out/{tumour}.mutect2.filter.norm.annot.vcf.gz"
  output:
    "out/{tumour}.mutect2.somatic.af.png"
  shell:
    "src/plot_af.py --log --sample {wildcards.tumour} --target {output} --dp {config[dp_threshold]} < {input}"

#####################
# loh
rule loh:
  input:
    snvs="out/{tumour}.strelka.somatic.snvs.af.norm.annot.vcf.gz", # loh requires strelka for now
    indels="out/{tumour}.strelka.somatic.indels.norm.annot.vcf.gz"
  output:
    tsv="out/{tumour}.loh.tsv",
    bed="out/{tumour}.loh.bed"
  log:
    stderr="log/{tumour}.loh.stderr"
  params:
    tumour="{tumour}",
    regions=' '.join(config["regions_of_interest"]),
    region_names=' '.join(config["genes_of_interest"]),
    region_padding=' '.join([str(x) for x in config["loh_region_padding"]])
  shell:
    "(tools/loh_caller-{config[loh_version]}/loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 20 --neutral --min_af 0.1 < {input.snvs} > tmp/{params.tumour}.loh.snvs.tsv && "
    "tools/loh_caller-{config[loh_version]}/loh.py --germline NORMAL --tumour TUMOR --filtered_variants --min_dp_germline 10 --min_dp_tumour 20 --neutral --min_af 0.1 < {input.indels} > tmp/{params.tumour}.loh.indels.tsv && "
    "sort -k1,1 -k2,2n tmp/{params.tumour}.loh.snvs.tsv tmp/{params.tumour}.loh.indels.tsv > {output.tsv} && "
    "tools/loh_caller-{config[loh_version]}/loh_merge.py --verbose --noheader --min_len 1000 --min_prop 0.1 --plot out/{params.tumour}.loh --regions {params.regions} --region_names {params.region_names} --region_padding {params.region_padding} --plot_chromosomes <{output.tsv} >{output.bed}) 2>{log.stderr}"

#####################
# cnv
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
    stderr="log/{tumour}.msisenser.stderr"
  params:
    tumour="{tumour}",
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
    exposures="out/{tumour}.mutational_signature_v3_sbs.filter.exposures"
  output:
    vcf="out/{tumour}.intersect.pass.filter.signatures.vcf.gz",
    png="out/{tumour}.signatures_cosmic_v3_sbs.png"
  shell:
    "{config[module_htslib]} && "
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/annotate_context.py --genome {input.reference} --vcf {input.vcf} > tmp/{wildcards.tumour}.intersect.pass.filter.context.vcf.gz && " # annotates
    "python tools/mutational_signature-{config[signature_version]}/mutational_signature/assign_signatures.py --definition tools/mutational_signature-{config[signature_version]}/data/signatures_cosmic_v3_sbs.txt --signatures {input.exposures} --vcf tmp/{wildcards.tumour}.intersect.pass.filter.context.vcf.gz --threshold 0.01 --plot {output.png} > {output.vcf}"

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
    "src/mutation_rate.py --verbose --vcfs {input.vcfs} --bed {input.regions} --signature_artefacts tools/mutational_signature-{config[signature_version]}/data/signature_summary.tsv --signature_artefact_penalty 0.5 >{output} 2>{log.stderr}"

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

### not working
# tumour purity - currently not working
rule purity:
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
    "samtools view -h -b {input.bam} {config[regions_of_interest]} > {output} && "
    "samtools index {output}"

# somatic
rule somatic_variant_summary:
  input:
    vcfs1=expand("out/{tumour}.intersect.pass.filter.vcf.gz", tumour=samples['tumours']),
    vcfs2=expand("out/{tumour}.strelka.somatic.indels.norm.vep.pass.vcf.gz", tumour=samples['tumours']),
    lohs=expand("out/{tumour}.loh.bed", tumour=samples['tumours'])

  output:
    tsv="out/aggregate/targetted_gene_summary.somatic.tsv",
    png="out/aggregate/targetted_gene_summary.somatic.png"

  shell:
    "src/targetted_gene_summary.py --vcfs {input.vcfs1} {input.vcfs2} --lohs {input.lohs} --genes {config[genes_of_interest]} --loci {config[regions_of_interest]} > {output.tsv} && "
    "src/targetted_gene_plot.py --target {output.png} < {output.tsv}"


# somatic
rule germline_variant_summary:
  input:
    vcf="out/aggregate/germline_joint.hc.normalized.annot.vcf.gz"

  output:
    tsv="out/aggregate/targetted_gene_summary.germline.tsv",
    png="out/aggregate/targetted_gene_summary.germline.png"

  shell:
    "src/targetted_gene_summary.py --vcfs {input.vcf} --genes {config[genes_of_interest]} --loci {config[regions_of_interest]} --multisample > {output.tsv} && "
    "src/targetted_gene_plot.py --target {output.png} < {output.tsv}"

