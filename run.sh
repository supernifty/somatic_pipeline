#!/bin/bash

# run on the cluster
MAXJOBS=512
PREFIX="pg"

echo "cleaning up..."

[ -e log.bck ] && rm -r log.bck
[ -e log ] && mv log log.bck

mkdir -p log
rm log/*

mkdir -p tmp
rm -r tmp/*

mkdir -p out/aggregate

# dry
echo "dry run..."
#snakemake --verbose -n -j $MAXJOBS --cluster-config cfg/cluster.json --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"
snakemake --verbose -n -j $MAXJOBS --cluster-config cfg/cluster.json --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"
echo "return to start pipeline, ctrl-c to quit"
read -n 1 -p "Continue?"

# real
echo "starting live run at $(date)..."
#snakemake -p -j $MAXJOBS --cluster-config cfg/cluster.json --stats log/snakemake_stats.json --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n}  -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"
snakemake -p -j $MAXJOBS --cluster-config cfg/cluster.json --stats log/snakemake_stats.json --rerun-incomplete --jobname "${PREFIX}-{rulename}-{jobid}" --cluster "sbatch -A {cluster.account} -p {cluster.partition} --ntasks={cluster.n} --nodes={cluster.nodes} -t {cluster.time} --mem={cluster.memory} --output=log/slurm-%j.out --error=log/slurm-%j.out"

echo "finished at $(date)"
