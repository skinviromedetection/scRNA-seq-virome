# scRNA-seq-virome

## Getting started
This workflow is for scRNA-seq experiments run on the 10X Genomics Platform. The first step is the align output .fastq files using CellRanger (https://www.10xgenomics.com/support/software/cell-ranger/latest). Reference genomes can be downloaded from 10X genomics (https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes). CellRanger The version must be greater than CellRanger 7.0 such that the alignment output of CellRanger contains unmapped reads.

The other file needed for this workflow is a file of viral genomes of interest. For this work we have used a previously published viral genome database (Selitsky et al, 2020) containing 1893 vertebrate viral genomes, masked at repetitive locations as well as locations that align to the human genome, available at [https://github.com/dmarron/virdetect/tree/master/reference](https://github.com/dmarron/virdetect/tree/master/reference)

## Align data to human genome
> cellranger count --transcriptome=/path-to-transcriptome/ --fastqs=/path-to-fastq-folder/ --id=SAMPLENAMES --nopreflight --jobmode=local --localcores=32 --localmem=128 --nosecondary --sample LANE1.fastq,LANE2.fastq,LANE3,fastq,LANE4.fastq
>
Typically after running CellRanger the output includes a file with path SAMPLE/outs/possorted_genome_bam.bam. This is the input file for downstream steps

## Extract reads that did not map to human genome
> samtools fastq -f 4 possorted_genome_bam.bam -T CR,CY,CB,UR,UY,UB > ${output}_unmapped.fastq
>
> 
> ## Extract reads that did not map to human genome
