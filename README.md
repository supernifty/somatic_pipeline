
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
