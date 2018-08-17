#!/usr/bin/env bash
#
# usage:
# $0 input_vcf output_vcf reference threads

set -e

echo "starting vep at $(date)"

module load perl/5.18.0

INPUT=$1
OUTPUT=$2
REFERENCE=$3
THREADS=$4

THREADS=1 # ignore threads parameter due to vep errors

VEPPATH=/vlsci/UOM0040/shared/km/programs/ensembl-vep
CACHE=/vlsci/UOM0040/shared/km/programs/ensembl-vep/data/
export PERL5LIB=$PERL5LIB:/vlsci/VR0002/kmahmood/Programs/vep/vep_87/ensembl-tools-release-87/scripts/variant_effect_predictor:tools/vep

$VEPPATH/vep \
    --cache \
    --refseq \
    --offline \
    --dir_cache $CACHE \
    --fasta $REFERENCE \
    -i $INPUT \
    -o ${OUTPUT}.tmp.vcf \
    --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs \
    --exclude_predicted \
    --af_gnomad \
    --format vcf \
    --force_overwrite --vcf \
    --fields Consequence,IMPACT,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,HGVSc,HGVSp,cDNA_position,CDS_position,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MaxEntScan_alt,MaxEntScan_diff,MaxEntScan_ref,PICK \
    --fork $THREADS \
    --flag_pick \
    --plugin MaxEntScan,/vlsci/UOM0040/shared/km/programs/ensembl-vep/data/MaxEntScan/

bgzip < ${OUTPUT}.tmp.vcf > $OUTPUT
rm ${OUTPUT}.tmp.vcf

echo "finishing vep at $(date)"
