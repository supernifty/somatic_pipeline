# Somatic pipeline

## Installation
* Python 3 is required.

Additional libraries required:
* curl

```
python -m venv somatic_venv
. ./somatic_venv/bin/activate
pip install -r requirements.txt
```

## Spartan
Before starting:
```
module load Python/3.6.1-intel-2017.u2
module load cURL/7.58.0-intel-2017.u2
module load icc
. ../../../software/venv_somatic_2/bin/activate
```

In config.yaml
```
module_bedtools: 'module load BEDTools/2.26.0-intel-2016.u3'
module_bwa: 'module load BWA/0.7.17-intel-2016.u3'
module_htslib: 'module load HTSlib/1.8-intel-2017.u2'
module_java: 'module load Java/1.8.0_152'
module_python2: 'module load Python/2.7.13-intel-2017.u2'
module_R: 'module load R/3.5.0-GCC-6.2.0'
module_samtools: 'module load SAMtools/1.8-intel-2016.u3-HTSlib-1.8'
module_network: 'module load web_proxy'
```

## Dependencies
* reference/genome.fa: this file needs to be bwa indexed.
* reference/msi.regions.bed: TODO
* reference/regions.bed: TODO
* reference/mantis.bed: ./tools/mantis/tools/RepeatFinder -i genome.fa -o mantis.bed
* msiseq
```
  # using tandem_repeats_nim
  ./main --repeat=1 --min=5 < ../somatic_pipeline/reference/genome.fa > hg19.repeats.1.bed
  sed 's/ [^\t]*//' < /data/cephfs/punim0567/peter/src/tandem_repeats_nim/hg19.repeats.1.bed > ./reference/hg19repeats.msiseq.tsv 
```

* msiseq - old instructions
```
  wget http://steverozen.net/data/Hg19repeats.rda
  Rscript
  > load('Hg19repeats.rda')
  > write.table(Hg19repeats, file='Hg19repeats.tsv', quote=FALSE, sep='\t', row.names=F)
  > quit()
```
  * had to update exp formatted numbers

```
  sed '1d' < reference/Hg19repeats.tsv | bedtools sort | awk '{ $3 += 1; print }' | sed 's/ /\t/g' > reference/Hg19repeats.sorted.tsv
```

* revel
```
cd references
wget https://rothsj06.u.hpc.mssm.edu/revel_grch38_all_chromosomes.csv.zip
unzip < revel_grch38_all_chromosomes.csv.zip > revel_grch38_all_chromosomes.csv.gz
```

Modules
* bwa
* java
* samtools

Tools directory
* picard

## Install
run install.sh to install:
* strelka
* cnv_caller
* fastqc
* gatk 3.8.1
* gatk 4


### GATK 4 Bundle
```
wget "gsapubftp-anonymous@ftp.broadinstitute.org:/bundle/b37/1000G_omni2.5.b37.sites.vcf.*"
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/beta/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz
```

### DKFZBiasFilter
git clone https://github.com/supernifty/DKFZBiasFilter.git
pip install -r requirements.txt

### Conpair
```
git clone https://github.com/supernifty/Conpair.git
```

### pindel
wget https://github.com/genome/pindel/archive/v0.2.5b8.tar.gz
tar xvfj v0.2.5b8.tar.gz
cd pindel-0.2.5b8
module load gcc/6.4.0
module load htslib-intel/1.8
./INSTALL /usr/local/htslib/1.8-intel/

### Trimmomatic
```
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip
```

### verifyBamID
```
wget https://github.com/statgen/verifyBamID/releases/download/v1.1.3/verifyBamIDLibStatGen.1.1.3.tgz
cd tools/verifyBamID_1.1.3
make
```

### platypus
* request from http://www.well.ox.ac.uk/platypus
```
cd tools
git clone https://github.com/andyrimmer/Platypus
cd -
```

### varscan
cd tools
wget https://downloads.sourceforge.net/project/varscan/VarScan.v2.3.9.jar
cd -

The following R packages need to be installed:
```
#!/usr/bin/env Rscript

# installation
install.packages('optparse', repos = "http://cran.us.r-project.org")
source("http://bioconductor.org/biocLite.R");
biocLite("DNAcopy");
```

## Configuration

* cfg/config.yaml: set sample details
* cfg/cluster.json: set cluster resources

## Usage

```
./run.sh
```

## Outputs

### Variant Calls

#### Mutect2
* out/{tumour}.mutect2.vcf.gz: output from mutect2
* out/{tumour}.mutect2.filter.vcf.gz: application of mutect2 filter and BAM_DEPTH annotation
* out/{tumour}.mutect2.filter.norm.vep.vcf.gz: annotated normalised version of mutect2.filter.vcf.gz
* out/{tumour}.mutect2.filter.norm.vep.pass.vcf.gz: annotated normalised version of mutect2.filter.vcf.gz

#### Strelka
* out/{tumour}.strelka.somatic.snvs.vcf.gz: output from strelka
* out/{tumour}.strelka.somatic.snvs.af.vcf.gz: add AF and BAM_DEPTH annotation
* out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz: normalised
* out/{tumour}.strelka.somatic.snvs.af.norm.vep.vcf.gz: normalised and annotated
* out/{tumour}.strelka.somatic.snvs.af.norm.vep.pass.vcf.gz: normalised, annotated, pass

Indels
* out/{tumour}.strelka.somatic.indels.vcf.gz: output from strelka
* out/{tumour}.strelka.somatic.indels.norm.vcf.gz: normalised
* out/{tumour}.strelka.somatic.indels.norm.vep.vcf.gz: normalised, annotated
* out/{tumour}.strelka.somatic.indels.norm.vep.pass.vcf.gz: filtered on pass

#### Intersection
* out/{tumour}.intersect.vcf.gz: intersection of out/{tumour}.mutect2.filter.norm.vep.vcf.gz and out/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz
* out/{tumour}.intersect.pass.vcf.gz: intersected pass calls
* out/{tumour}.intersect.pass.filter.pass.vcf.gz: intersected filter on config AF and DP settings
* out/{tumour}.pass_one.vcf.gz: one of the two callers passes

### Directories
* cfg: configuration files
* in: input files (fastq)
* log: command logs
* out: generated files
* reference: reference files
* tools: 3rd party tools

### Informative
* out/{tumour}.strelka.somatic.af.png: AF distribution
* out/{tumour}.mutect2.somatic.af.png: AF distribution

## Add a new batch
* mkdir batchname
* want to link tmp to /data/scratch/projects/punim0567/tmp-batchname:
  * mkdir -p /data/scratch/projects/punim0567/tmp-batchname && ln -s /data/scratch/projects/punim0567/tmp-batchname tmp
* ./deploy.sh full_path_to_batchname
* symlink fastqs into in directory: for merging can use ./util/merge_fastq.py
  * fastq format is expected to be: in/S1_RG_2_R1.fastq.gz
  * if sample looks like sample_BC then sample components is 2. If it looks like sample-BC then sample components is 1.
  * e.g.
```
python /data/gpfs/projects/punim0567/peter/src/somatic_pipeline/util/merge_fastq.py \
  --files /data/projects/punim0567/data/AGRF_CAGRF20073103_HKVGTDSXY/AGRF_CAGRF20073103_HMCG2DSXY_additional_data/*_R*.fastq.gz \
    /data/projects/punim0567/data/AGRF_CAGRF20073103_HKVGTDSXY/AGRF_CAGRF20073103_HKVGTDSXY/*_R*.fastq.gz \
  --sample_components 2 \
  --outdir in
```
  * e.g. for f in /data/gpfs/projects/punim0567/data/AGRF_CAGRF20073103_HJV72DSXY/merged/*.fastq.gz; do ln -s $f .; done
* update samples.yaml: can use ./util/samples.py > ./cfg/samples.yaml
  * e.g. /data/gpfs/projects/punim0567/peter/src/somatic_pipeline/util/samples.py --tumours in/*_T_*_R* in/*_T[0-9]_*_R* --normals in/*_BC_*_R* --sample_components 2 --verbose > cfg/samples.yaml
  * e.g. /data/gpfs/projects/punim0567/peter/src/somatic_pipeline/util/samples.py --tumours in/*_T*_R* --normals in/*_BC*_R* > cfg/samples.yaml
* ./run.sh

## Add new samples to a batch
```
rm out/aggregate/* out/germline_joint_*.vcf out/tumour_joint_*.vcf out/germline_joint.hc.normalized.vcf out/tumour_joint.hc.normalized.vcf out/mutect2.pon.vcf.gz
```

## TODO
* extract all high impact to tsv
* include commonality in extracted tsv
