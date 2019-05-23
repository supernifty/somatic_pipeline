#!/usr/bin/env bash

set -o errexit

# remove files so that a new sample can run
BCK=out.bck.$(date +%Y-%m-%d)
echo "backing up to $BCK..."

mkdir -p $BCK

# pon means all mutect is invalid
mv out/*mutect* $BCK
mv out/*intersect* $BCK
mv out/*pass_one* $BCK
mv out/aggregate $BCK
mv out/germline_joint* $BCK

echo "backing up to $BCK: done"
