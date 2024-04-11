#!/bin/bash

wget http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#wget http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
#gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz

#Reference files needed for CICERO fusion caller tool
curl -o reference.tar.gz 'https://zenodo.org/record/5088371/files/reference.tar.gz?download=1' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7' \
  -H 'Accept-Language: sr-RS,sr;q=0.9,en-US;q=0.8,en;q=0.7,hr;q=0.6,bs;q=0.5' \
  -H 'Connection: keep-alive' \
  -H 'Cookie: session=2980316b97ce4005_63ff4ad2.mAoUXmba28yXDd_NThYLeN2Bydw; __atuvc=1%7C12; _pk_id.57.a333=0d4f36c1567b011d.1679494278.3.1679504177.1679496774.; _pk_ses.57.a333=*' \
  -H 'Referer: https://zenodo.org/record/5088371' \
  -H 'Sec-Fetch-Dest: document' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Sec-Fetch-User: ?1' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/111.0.0.0 Safari/537.36' \
  -H 'sec-ch-ua: "Google Chrome";v="111", "Not(A:Brand";v="8", "Chromium";v="111"' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'sec-ch-ua-platform: "Windows"' \
  -H 'sec-gpc: 1' \
  --compressed
tar -xzf reference.tar.gz
rm reference.tar.gz


#Index genome in GEM format
#Transcriptome index in GEM format
#Keys convert transcriptome-genome coordinates
~/projects/alkFus_bench/downloads/fusion_callers/ChimPipe/bin/gemtools-1.7.1-i3/gemtools index -i Homo_sapiens.GRCh38.dna.primary_assembly.fa -t 8
~/projects/alkFus_bench/downloads/fusion_callers/ChimPipe/bin/gemtools-1.7.1-i3/gemtools t-index -i Homo_sapiens.GRCh38.dna.primary_assembly.gem -a Homo_sapiens.GR>

#Transcriptome GRCh38 v103
wget http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz


#STAR index for GRCh38
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STARindex_GRCh38 --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.v103.gtf
