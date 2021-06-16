#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 168:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 15   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Wagad"   # job name
#SBATCH --mem=30G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-WAGAD.oneCode.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

cd B.dir/

perl ../build_dictionary.pl --unknown --rm anom_juiced.V1.1.fasta.out > B.dictionary
perl ../one_code_to_find_them_all.pl --rm anom_juiced.V1.1.fasta.out --ltr B.dictionary --unknown
grep "###" *.copynumber.csv > B.summary




