make
./cgal genome.scf.fasta
sort -n -k1,1 -k2,2 info.txt > infoOutput.txt
wc -l infoOutput.txt
#commited
