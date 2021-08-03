# Scripts relating to the genome assembly of _Gossypium anomalum_

**Overview**

This repository contains scripts and/or bioinformatic details regarding the analysis of the _G. anomalum_ genome. The directories include:

_RepeatExplorerAnalysis_: R script (Repeat.B.R) for analyzing the RepeatExplorer output (B.clusters.counts.txt) and high-resolution versions of the output images it generates

_RepeatMaskerAnalysis_: slurm submission scripts for running RepeatMasker/One Code and to make the output table. Submission scripts load any necessary programs (as modules) and include program versions. `build_dictionary.pl` and `one_code_to_find_them_all.pl` contain the versions of the code from http://doua.prabi.fr/software/one-code-to-find-them-all which were used here

_annotations_: run.sh contains a record of what was run to generate annotations. Also links out to a GitLab pipeline (https://gitlab.com/IGBB/pipelines/annotation)

_dotplots_: contains high-resolution versions of the dotplots of _G. anomalum_ versus other genomes, the code used to generate the dotplots, and the list of genomes/accessions included. Also links out to a GitLab pipeline (https://gitlab.com/IGBB/pipelines/dotplots)

_introgression_: contains a list of genomes/reads and how these were used to infer introgression. High-resolution images are included.

_repeat-clustering_: contains information regarding the species/samples clustered via RepeatExplorer, including SRA numbers, genome size, number of reads clustered, and the code used. Also links out to a GitLab pipeline (https://gitlab.com/IGBB/pipelines/repeat-clustering)
\
\
\
\
**About the assembly**

We report a high-quality de novo genome assembly for _G. anomalum_ (B1). This genome was initially assembled using 55x coverage of PacBio reads, yielding a draft assembly of 229 contigs with an N50 of 11 Mb. HiC reads (140.5 million) were integrated to the final, contiguous assembly, consisting of 13 chromosomes with an average length of 92 Mb and containing only 20.8 kb (0.002%) gap sequence within the chromosomal scaffolds. The total assembly is composed of 68 contigs (including the 13 chromosomal scaffolds) with a total length of 1193 Mb, ~88% of the estimated 1359 Mb genome (Hendrix and Stewart 2005). The final assembly N50 is 97.68 Mb and N90 is 73.95 Mb

The genome and related annotation files are available through NCBI (PRJNA421337) and CottonGen. This genome report is published here:

Grover, et al., (in revision) The Gossypium anomalum genome as a resource for cotton improvement and evolutionary analysis of hybrid incompatibility. G3: Genes, Genomes, Genetics 
