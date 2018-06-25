# Somatic pipeline

## Installation
* Python 3 is required.

```
python -m venv somatic_venv
. ./somatic_venv/bin/activate
pip install -r requirements.txt
```

## Dependencies
* reference/genome.fa: this file needs to be bwa indexed.

Modules
* bwa
* java
* samtools

Tools directory
* picard

### Strelka ###
```
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
bunzip2 < strelka-2.9.2.centos6_x86_64.tar.bz2 | tar xvf -
```

### fastqc ###
```
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod +x tools/FastQC/fastqc
```

### GATK 4 ###
```
wget https://github.com/broadinstitute/gatk/releases/download/4.0.0.0/gatk-4.0.0.0.zip
```

### GATK 3.8.1 ###
```
wget 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef'
tar xvfj GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 
```

### GATK 4 Bundle
```
wget "gsapubftp-anonymous@ftp.broadinstitute.org:/bundle/b37/1000G_omni2.5.b37.sites.vcf.*"
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/beta/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz
```

### DKFZBiasFilter
git clone https://github.com/supernifty/DKFZBiasFilter.git
pip install -r requirements.txt

### Conpair ###
```
git clone https://github.com/supernifty/Conpair.git
```

### verifyBamID ###
```
wget https://github.com/statgen/verifyBamID/releases/download/v1.1.3/verifyBamIDLibStatGen.1.1.3.tgz
cd tools/verifyBamID_1.1.3
make
```

## Configuration

* cfg/config.yaml: set sample details
* cfg/cluster.json: set cluster resources

## Usage

```
./run.sh
```

## Directories
* cfg: configuration files
* in: input files (fastq)
* log: command logs
* out: generated files
* reference: reference files
* tools: 3rd party tools
