#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH --cpus-per-task=150
#SBATCH --mem=300G 
#SBATCH --time=12-02:30:00
#SBATCH --output=job.BERM.%J.out 
#SBATCH --error=job.BERM.%J.err
#SBATCH -J "BERM" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 


module load repeatmasker/4.0.7-openmpi3-vruhfro

RepeatMasker -pa 75 -s -nolow -xsmall -lib cottonRB23.04.fa -no_is -frag 500000 -nocut -noisy -dir B.dir -html -gff anom_juiced.V1.1.fasta
RepeatMasker -pa 75 -s -nolow -xsmall -lib cottonRB23.04.fa -no_is -frag 500000 -nocut -noisy -dir E.dir -html -gff stocksii.V1.0.hic_juiced.fasta

ml >> job.modules.RM.list


