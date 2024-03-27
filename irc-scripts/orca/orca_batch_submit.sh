#!/bin/bash -f

#  Sun Grid Engine directives:
#$ -j y
#$ -l h_cpu=12,vf=16000M,rslots=8

export ORCA={OrcaPath}
export PATH=$HOME/bin:$ORCA:$PATH
export LD_LIBRARY_PATH=$HOME/lib:$ORCA:$LD_LIBRARY_PATH

echo running {FileName}.inp on `hostname`.
cd {FilePath}
{OrcaPath}/orca {FileName1}.inp
{OrcaPath}/orca_2mkl {FileName1} -molden
