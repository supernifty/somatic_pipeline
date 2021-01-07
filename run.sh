#!/bin/bash

set -o errexit

# run on the cluster
MAXJOBS=32
PREFIX="tl-$(basename `pwd`)"

export R_LIBS_USER="$(pwd)/tools/R_libs"

echo "cleaning up..."

[ -e log.bck ] && rm -r log.bck
[ -e log ] && mv log log.bck

mkdir -p log
touch log/dummy
rm log/*

mkdir -p tmp
touch tmp/dummy
rm -r tmp/*

mkdir -p out/aggregate
chmod -R +w out

# dry
echo "dry run..."
#snakemake --verbose -n -j $MAXJOBS --cluster-config cfg/cluster.yaml --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"
snakemake --verbose -n -j $MAXJOBS --cluster-config cfg/cluster.yaml --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"
echo "return to start pipeline, ctrl-c to quit"
read -n 1 -p "Continue?"

# real
echo "starting live run at $(date)..."
#snakemake -p -j $MAXJOBS --cluster-config cfg/cluster.yaml --stats log/snakemake_stats.json --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"
snakemake -p -j $MAXJOBS --cluster-config cfg/cluster.yaml --stats log/snakemake_stats.json --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n} --nodes={cluster.nodes} -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"

echo "finished at $(date)"

# remove intermediate files
# rm out/*.fq.gz
# rm out/*.bqsr.bam

# mark as read only
#echo "marking read only at $(date)"
#chmod -R -w out
