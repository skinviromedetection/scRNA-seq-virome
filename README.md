# scRNA-seq-virome

## Getting started
This workflow is for scRNA-seq experiments run on the 10X Genomics Platform. The first step is the align output .fastq files using CellRanger (https://www.10xgenomics.com/support/software/cell-ranger/latest). Reference transcriptomes can be downloaded from 10X genomics (https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes). CellRanger The version must be greater than CellRanger 7.0 such that the alignment output of CellRanger contains unmapped reads.

The other file needed for this workflow is a file of viral genomes of interest. For this work we have used a previously published viral genome database (Selitsky et al, 2020) containing 1893 vertebrate viral genomes, masked at repetitive locations as well as locations that align to the human genome, available at [https://github.com/dmarron/virdetect/tree/master/reference](https://github.com/dmarron/virdetect/tree/master/reference). 

Finally this workflow requires blastn from the NCBI BLAST+ library of tools: [https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
## Align data to human genome
``` cellranger count --transcriptome=</path-to-transcriptome/> --<fastqs=/path-to-fastqs-folder/> --id=<SAMPLENAME> --nopreflight --jobmode=local --localcores=32 --localmem=128 --nosecondary --sample <LANE1,LANE2,LANE3,LANE4> ```

Typically after running CellRanger the output includes a file with path SAMPLE/outs/possorted_genome_bam.bam. This is the input file for downstream steps

## Extract reads that did not map to human genome
The next step filters the output of CellRanger count (possorted_genome_bam.bam) to extract the reads that did not align to the human genome. The -T flag is used to specifiy a list of tags to be included in the output. These tags are necessary in order to match reads to the cells from which they originated. 

``` samtools fastq -f 4 <possorted_genome_bam.bam> -T CR,CY,CB,UR,UY,UB > <output_unmapped.fastq> ```

 ## Prepare files for BLAST
 Next we prepare the unmapped reads file to be aligned to the database of viral genomes as BLAST only accepts fasta files:
 
``` sed -n '1~4s/^@/>/p;2~4p' <output_unmapped.fastq> | sed 's/\t//g' > <output.fasta>```

We also need to create a BLAST database from the viral genome database. To do this we use the makeblastdb tool from the NCBI BLAST+ library of tools: [https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)

``` makeblastdb -in <your_genome_file.fa> -dbtype nucl -out <db_name>```

## Align unmapped reads to viral genomes
Next we run blastn to align unmapped reads to viral genomes: 

``` blastn -query <output.fasta> -db <path-to-makeblastdb-output/db_name.db> -out <output_BLAST.txt> -outfmt 6 ```

The output of this is a text file with columns as below:

```
Column headers:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
 1.  qseqid      query or source (gene) sequence id
 2.  sseqid      subject or target (reference genome) sequence id
 3.  pident      percentage of identical positions
 4.  length      alignment length (sequence overlap)
 5.  mismatch    number of mismatches
 6.  gapopen     number of gap openings
 7.  qstart      start of alignment in query
 8.  qend        end of alignment in query
 9.  sstart      start of alignment in subject
 10.  send        end of alignment in subject
 11.  evalue      expect value
 12.  bitscore    bit score
```

From this output file we can see what viral genome each read mapped to and various scoring metrics to assess the quality of the alignment. This file can be filtered to ensure high quality alignments by putting requirements on metrics including evalue, pident, length, mismatch,gapopen.

## Integrating with scRNA-seq data analysis

Next we are interested in integrating this data with single cell data analysis. To do this step, it requires that the user has already analyzed their scRNA-seq data using the R Seurat package [https://satijalab.org/seurat/](https://satijalab.org/seurat/), including filtering steps, and have created a seurat object, which we will refer to as ```filtered_seurat_obj``` for downstream work.

In this workflow we first filter the aligments in the <output_BLAST.txt> to only include reads that originate from cells that pass QC filterers after analysis with Seurat in R. Commonly used QC metrics include removing cells with few genes and removing cells with a high percentage of reads mapping to the mitochondrial genome. From the ```filtered_seurat_obj``` of filtered cells, we extract the barcodes of all cells that passed the user-degined filters and then filter the <output_BLAST.txt> file to include only reads that originate from the filtered cells. An example of R code for extracting the barcodes into a file filtered_barcodes.txt is below:

```
seurat_data <- Read10X(data.dir = args)
seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                 min.features = 10, 
                                 project = sample_name)

# Add number of genes per UMI for each cell to metadata
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

# Compute percent mito ratio
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- seurat_obj@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
seurat_obj@meta.data <- metadata

filtered_seurat <- subset(x = seurat_obj, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.2))

write.table(filtered_seurat$cells, file = paste0("filtered_barcodes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")
```

We can then filter <output_BLAST.txt> using the command line tool grep:

```grep -f filtered_barcodes.txt <output_BLAST.txt> > <output_BLAST_filtered.txt>```

This file can then be used for downstream analysis including counting number of reads that align to viral genomes.

Finally, we can add a metadata column to ```filtered_seurat_obj``` in R that indicates whether each cell contained reads aligning to the viral genome, which allows us to filter cells by those that do or do not have viral transcripts. To do this we extract the cell barcodes from reads in  <output_BLAST_filtered.txt>. The user can do this a number of ways; below is code for obtaining the cell barcodes of reads that map to Rubella virus using in python3 into output file output_Rubellapos_barcodes.txt. 

```
import pandas as pd
import glob

fil = <output_BLAST_filtered.txt>
barcodes=pd.read_csv(fil,sep='\t',header=None)
barcodes=list(barcodes[barcodes[1].str.contains('Rubella')][0])
barcodes=[n for n in barcodes if 'CB:Z:' in n]
barcodes = [n.split('CB:Z:')[1].split('-1UR')[0] for n in barcodes]
with open('BLAST_MCC/'output_Rubellapos_barcodes.txt', 'w') as file:
    if len(barcodes)>0:
        for word in barcodes:
            file.write(word + '\n')
```

Then in R, the metadata column can be added to ```filtered_seurat_obj``` as such:

```
Rubella_codes <- readLines('<output_Rubellapos_barcodes.txt>')

contains_rubella <- function(value, codes){
  any(sapply(codes, function(pattern) grepl(pattern, value)))
}

filtered_seurat_obj@meta.data$Rubella  <- ifelse(sapply(filtered_seurat_obj@meta.data$cells, contains_Rubella, Rubella_codes), "yes", "no")
```

