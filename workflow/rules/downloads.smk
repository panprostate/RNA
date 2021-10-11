rule download_resources:
    output:
        reference="resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta",
        dictn="resources/Homo_sapiens_assembly19_1000genomes_decoy.dict",
        faidx="resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta.fai",
        transcripts="resources/gencode.v38lift37.transcripts.fa",
        tx2gene="resources/tx2gene.tsv.gz",
        gtf="resources/gencode.v38lift37.annotation.gtf",
        paSites="resources/pasites_hg19.bed",
        fc="resources/FANTOM_CAT.lv3_robust.gtf",
        blacklist="resources/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz",
        knowFusions="resources/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz",
        proteinDomains="resources/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3",
        dbSNP="resources/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf",
        dbSNPidx="resources/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf.idx",
        knowIndels="resources/Mills_and_1000G_gold_standard.indels.b37.sites.vcf",
        knowIndelsidx="resources/Mills_and_1000G_gold_standard.indels.b37.sites.vcf.idx"
    priority: 3
    threads: 1
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/downloads/resources.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        set -e
        echo "Downloading resources..."
        [[ -f {output.dbSNP} ]] || gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf {output.dbSNP}
        [[ -f {output.dbSNPidx} ]] || gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf.idx {output.dbSNPidx}
        [[ -f {output.knowIndels} ]] || gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Mills_and_1000G_gold_standard.indels.b37.sites.vcf {output.knowIndels}
        [[ -f {output.knowIndelsidx} ]] || gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Mills_and_1000G_gold_standard.indels.b37.sites.vcf.idx {output.knowIndelsidx}
        [[ -f {output.reference} ]] || gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta {output.reference}
        [[ -f {output.faidx} ]] || gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta.fai {output.faidx}
        [[ -f {output.dictn} ]] || gsutil cp gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dict {output.dictn}
        [[ -f {output.paSites} ]] || wget -O resources/pasites_hg19.bed.gz https://www.dropbox.com/s/i55pzcuh6113vvq/pasites_hg19.bed.gz
        [[ -f {output.tx2gene} ]] || wget -O resources/tx2gene.tsv.gz https://www.dropbox.com/s/ays94bytt7r02ex/tx2gene.tsv.gz
        [[ -f {output.fc} ]] || wget -O resources/FANTOM_CAT.lv3_robust.gtf.gz https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.gtf.gz
        [[ -f {output.gtf} ]] || wget -O resources/gencode.v38lift37.annotation.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz
        [[ -f {output.transcripts} ]] || wget -O resources/gencode.v38lift37.transcripts.fa.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.transcripts.fa.gz
        [[ -f {output.blacklist} ]] || wget -O resources/arriba_v2.1.0.tar.gz https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz
        tar -xvf resources/arriba_v2.1.0.tar.gz --directory resources/
        mv resources/arriba_v2.1.0/database/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz resources/arriba_v2.1.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz resources/arriba_v2.1.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3 resources/
        rm -rf resources/arriba_v2.1.0/ resources/arriba_v2.1.0.tar.gz
        echo "Unzipping some files..."
        gunzip resources/gencode.v38lift37.annotation.gtf.gz resources/gencode.v38lift37.transcripts.fa.gz resources/FANTOM_CAT.lv3_robust.gtf.gz resources/pasites_hg19.bed.gz
        sed -i 's/^chr//g' resources/gencode.v38lift37.annotation.gtf
        sed -i 's/^chr//g' resources/FANTOM_CAT.lv3_robust.gtf
        echo "Done!"
        """

rule star_index_download:
    output:
        index = "resources/STARindex_hg19/SA"
    priority: 3
    threads: 1
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/STARindex/STARindex.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        set -e
        wget -O resources/STARindex_hg19.tar https://www.dropbox.com/s/7zf1ad0hokpizfi/STARindex_hg19.tar
        tar -xvf resources/STARindex_hg19.tar --directory resources/
        """

rule salmon_index_download:
    output:
        "resources/salmon_hg19/ctable.bin"
    threads: 1
    priority: 3
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/salmonIndex/salmonIndex.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        set -e
        wget -O resources/salmon_hg19.tar https://www.dropbox.com/s/mybii6z71v9pt57/salmon_hg19.tar
        tar -xvf resources/salmon_hg19.tar --directory resources/
        """

rule download_ctatLib:
    output:
        ctat = "resources/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
    threads: 1
    priority: 2
    resources:
        mem_mb = 1000,
        runtime_min = "24:00:00"
    log:
        "logs/STARfusion/download_CTAT.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        wget -O resources/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
        tar -xzvf resources/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz --directory resources/
        """
