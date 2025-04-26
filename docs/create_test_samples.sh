iss generate --genomes temp_conf0_by_sample/iter_1_k_2_PSVs_combined.fna --n_reads 1000000 --model NovaSeq --cpus 48 --output temp_conf0_by_sample/iter_1XX1_PSVsXXk_2_iss --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes GCF_002263795.3_ARS-UCD2.0_genomic.fna --n_reads 10000 --model NovaSeq --cpus 48 --output S1_test_bos --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes GCF_002263795.3_ARS-UCD2.0_genomic.fna --n_reads 5000 --model NovaSeq --cpus 48 --output S2_test_bos --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes GCF_002263795.3_ARS-UCD2.0_genomic.fna --n_reads 1000 --model NovaSeq --cpus 48 --output S3_test_bos --compress --quiet --fragment-length 311 --fragment-length-sd 50 


iss generate --genomes megares_database_v3.00.fasta --n_reads 1000 --model NovaSeq --cpus 24 --output S1_test_AMR --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes megares_database_v3.00.fasta --n_reads 5000 --model NovaSeq --cpus 24 --output S2_test_AMR --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes megares_database_v3.00.fasta --n_reads 10000 --model NovaSeq --cpus 24 --output S3_test_AMR --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes GCF_007965905.1_ASM796590v1_genomic.fna --n_reads 5000 --model NovaSeq --cpus 24 --output S1_test_Mh --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes GCF_007965905.1_ASM796590v1_genomic.fna --n_reads 10000 --model NovaSeq --cpus 24 --output S2_test_Mh --compress --quiet --fragment-length 311 --fragment-length-sd 50 

iss generate --genomes GCF_007965905.1_ASM796590v1_genomic.fna --n_reads 1000 --model NovaSeq --cpus 24 --output S3_test_Mh --compress --quiet --fragment-length 311 --fragment-length-sd 50 


GCF_007965905.1_ASM796590v1_genomic.fna

cat S1_test_AMR_R1.fastq.gz S1_test_bos_R1.fastq.gz S1_test_Mh_R1.fastq.gz > ../test_pipeline/AMRplusplus/data/raw/S1_test_R1.fastq.gz
cat S1_test_AMR_R2.fastq.gz S1_test_bos_R2.fastq.gz S1_test_Mh_R2.fastq.gz > ../test_pipeline/AMRplusplus/data/raw/S1_test_R2.fastq.gz

cat S2_test_AMR_R1.fastq.gz S2_test_bos_R1.fastq.gz S2_test_Mh_R1.fastq.gz > ../test_pipeline/AMRplusplus/data/raw/S2_test_R1.fastq.gz
cat S2_test_AMR_R2.fastq.gz S2_test_bos_R2.fastq.gz S2_test_Mh_R2.fastq.gz > ../test_pipeline/AMRplusplus/data/raw/S2_test_R2.fastq.gz

cat S3_test_AMR_R1.fastq.gz S3_test_bos_R1.fastq.gz S3_test_Mh_R1.fastq.gz> ../test_pipeline/AMRplusplus/data/raw/S3_test_R1.fastq.gz
cat S3_test_AMR_R2.fastq.gz S3_test_bos_R2.fastq.gz S3_test_Mh_R2.fastq.gz > ../test_pipeline/AMRplusplus/data/raw/S3_test_R2.fastq.gz

# Variant++

cat S1_test_AMR_R1.fastq.gz S1_test_bos_R1.fastq.gz S1_test_Mh_R1.fastq.gz > ../test_pipeline/VARIANTplusplus/data/raw/S1_test_R1.fastq.gz
cat S1_test_AMR_R2.fastq.gz S1_test_bos_R2.fastq.gz S1_test_Mh_R2.fastq.gz > ../test_pipeline/VARIANTplusplus/data/raw/S1_test_R2.fastq.gz

cat S2_test_AMR_R1.fastq.gz S2_test_bos_R1.fastq.gz S2_test_Mh_R1.fastq.gz > ../test_pipeline/VARIANTplusplus/data/raw/S2_test_R1.fastq.gz
cat S2_test_AMR_R2.fastq.gz S2_test_bos_R2.fastq.gz S2_test_Mh_R2.fastq.gz > ../test_pipeline/VARIANTplusplus/data/raw/S2_test_R2.fastq.gz

cat S3_test_AMR_R1.fastq.gz S3_test_bos_R1.fastq.gz S3_test_Mh_R1.fastq.gz> ../test_pipeline/VARIANTplusplus/data/raw/S3_test_R1.fastq.gz
cat S3_test_AMR_R2.fastq.gz S3_test_bos_R2.fastq.gz S3_test_Mh_R2.fastq.gz > ../test_pipeline/VARIANTplusplus/data/raw/S3_test_R2.fastq.gz