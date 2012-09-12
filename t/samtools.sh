#!/bin/sh
file=$1; samtools faidx reference.fasta; samtools view -bt reference.fasta ${file}.sam -o ${file}.bam.un; samtools sort ${file}.bam.un $file; samtools index ${file}.bam.un; samtools mpileup -uf reference.fasta ${file}.bam | bcftools view -cg - > pileup/${file}.vcf; bgzip pileup/${file}.vcf; tabix -p vcf pileup/${file}.vcf.gz; rm ${file}.bam.un*;
