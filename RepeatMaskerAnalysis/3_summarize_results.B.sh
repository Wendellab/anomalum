cd B.dir

cat <(head -n1 A1.wagad.PGA.fasta.out_Chr1.copynumber.csv) *.summary > genomes.summary
sed -i 's/[.]fasta[.]out_/\tChr/g' genomes.summary
sed -i 's/Family/Genome\tFamily/' genomes.summary
sed -i 's/[.]copynumber[.]csv:\#\#*/\t/g' genomes.summary
sed -i 's/Family/Chr\tFamily/g' genomes.summary
sed -i 's/ChrChr/Chr/g' genomes.summary
sed -i 's/ChrD5/Chr/g' genomes.summary
sed -i 's/ChrD10/Chr/g' genomes.summary
sed -i 's/Chr_/Chr/g' genomes.summary
sed -i 's/[.]wagad[.]PGA//g' genomes.summary

