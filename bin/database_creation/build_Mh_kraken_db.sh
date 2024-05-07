

# Use NCBI-download
#Downloading Pasteurellaceae Refseq
ncbi-genome-download --taxids pasteurellaceae_taxids.txt -s refseq -p 2 --formats fasta -o pasteurellaceae_refseq bacteria

#Downloading Moraxellaceae Refseq
ncbi-genome-download --taxids moraxellaceae_taxids.txt -s refseq -p 2 --formats fasta -o moraxellaceae_refseq bacteria

#Downloading Mycoplasmataceae Refseq
ncbi-genome-download --taxids mycoplasmataceae_taxids.txt -s refseq -p 2 --formats fasta -o mycoplasmataceae_refseq bacteria


#Downloading Taxonomy Names-this line of command downloads all the taxonomy in the NCBI file
kraken2-build --download-taxonomy --db BRD_kraken_db

#Downloads Bacteria Refseq genomes
# https://github.com/DerrickWood/kraken2/issues/518
# Change line 46 in rsync_from_ncbi.pl from "ftp" to "https"
kraken2-build --download-library bacteria --db BRD_kraken_db/

#Adding BRD family genomes into kraken database-this line is moving files from individual directory to the database on our choice
find /s/angus/f/nobackup/Doster_tools/family_genome_database/*/ -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db BRD_kraken_db

#  family_genome_database/GCF_014947145.1/GCF_014947145.1_ASM1494714v1_genomic.fna.gz: unexpected end of file



#Building Kraken Database-This line of command works on building actual database by breaking our files into Kmers
#Threads can be changed similar to --ntasks up in the first line
kraken2-build --build --db BRD_kraken_db --threads 4
