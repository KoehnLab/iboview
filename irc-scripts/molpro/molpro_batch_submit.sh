#!/bin/bash -f

#  Sun Grid Engine directives:
#$ -j y
#$ -l h_cpu=12,vf=16000M,rslots=8

export LD_LIBRARY_PATH=$HOME/bin/lib:$LD_LIBRARY_PATH
export PATH=$HOME/dev/molpro/bin:$PATH

echo running %(FileName)s.inp on `hostname`.
%(MolproCmd)s

