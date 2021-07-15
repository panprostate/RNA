#!/bin/bash
echo "Downloading resources..."
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf .
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf.idx .
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf .
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Mills_and_1000G_gold_standard.indels.b37.sites.vcf .
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Mills_and_1000G_gold_standard.indels.b37.sites.vcf.idx .
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta .
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta.fai .
gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dict .
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.transcripts.fa.gz
echo "Unzipping some files..."
gunzip gencode.v38lift37.annotation.gtf.gz gencode.v38lift37.transcripts.fa.gz
sed -i -e 's/^chrM/MT/' -e 's/^chr//g' gencode.v38lift37.annotation.gtf
echo "Done!"
