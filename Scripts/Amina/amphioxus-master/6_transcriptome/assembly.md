# Data
We got data from Ferdi. Paired-end RNAseq (Illumina), three embryo stages.

# Genome
Got their genome assembly from: http://amphiencode.github.io/Data/. Can't do better for now. 

# Transcriptome assembly
Cannot trust or use their data (need transcripts, genes and no CDS for the Bgee pipeline). Asked Ferdi, he said they run default settings for Cufflinks (maybe true) and extracted the exons (so, not enough for me). I just rerun the assemblies: 

## Cufflinks 

First, align reads to the reference genome:
```
hisat2 -p30 -q --rna-strandness RF --dta-cufflinks -x ./hisat_idx -1 ../15hpf_A_1.fq.gz,../15hpf_B_1.fq.gz,../36hf_A_1.fq.gz,../36hf_B_1.fq.gz,../7hpf_A_1.fq.gz,../7hpf_B_1.fq.gz -2 ../15hpf_A_2.fq.gz,../15hpf_B_2.fq.gz,../36hf_A_2.fq.gz,../36hf_B_2.fq.gz,../7hpf_A_2.fq.gz,../7hpf_B_2.fq.gz -S hisat_amphio.sam
```

Run was successful, but curious log file:
```
163239054 reads; of these:
  163239054 (100.00%) were paired; of these:
    88312480 (54.10%) aligned concordantly 0 times
    50436953 (30.90%) aligned concordantly exactly 1 time
    24489621 (15.00%) aligned concordantly >1 times
    ----
    88312480 pairs aligned concordantly 0 times; of these:
      835129 (0.95%) aligned discordantly 1 time
    ----
    87477351 pairs aligned 0 times concordantly or discordantly; of these:
      174954702 mates make up the pairs; of these:
        143781171 (82.18%) aligned 0 times
        22561434 (12.90%) aligned exactly 1 time
        8612097 (4.92%) aligned >1 times
55.96% overall alignment rate
```

So basically half of the reads could not be aligned? Anyways we run assembly.

### Cufflinks with gtf

Output:

```
# transcripts
cat transcripts.gtf | cut -d ';' -f2 | sort | uniq | wc -l
231139
# genes
cat transcripts.gtf | cut -d ';' -f1 | cut -f9 | sort | uniq | wc -l
88632
```

Remove transcripts with fpkm=0:

```
cat transcripts.gtf | grep -v 'FPKM "0.0000000000"' > transcripts_nozero.gtf
# transcripts
cat transcripts_nozero.gtf | cut -d ';' -f2 | sort | uniq | wc -l
71171
# genes
cat transcripts_nozero.gtf | cut -d ';' -f1 | cut -f9 | sort | uniq | wc -l
51763
```

### Cufflinks without gtf

```
# transcripts
cat transcripts.gtf | cut -d ';' -f2 | sort | uniq | wc -l
48149
# genes
cat transcripts.gtf | cut -d ';' -f1 | cut -f9 | sort | uniq | wc -l
36171
```

Remove transcripts with fpkm=0:

```
cat transcripts.gtf | grep -v 'FPKM "0.0000000000"' > transcripts_nozero.gtf
# transcripts
cat transcripts_nozero.gtf | cut -d ';' -f2 | sort | uniq | wc -l
48149
# genes
cat transcripts_nozero.gtf | cut -d ';' -f1 | cut -f9 | sort | uniq | wc -l
36171
```

### Todo, maybe 

To rempap to reference genome? Or can we just use these dataset?

## Trinity

```
# transcripts 
cat Trinity.fasta | grep '^>' | wc -l 
576481
# genes
cat Trinity.fasta | grep '^>' | cut -f1 -d' ' | cut -f1-4 -d'_' | sort | uniq | wc -l
197755
```

## Compare to their 'final' assembly (Cufflinks+EVM)

```
# total transcripts
cat Bla_annot_final.gtf | cut -d ';' -f2 | sort | uniq | wc -l
218070
# transcripts only from cufflinks
cat Bla_annot_final.gtf | grep cuf | cut -d ';' -f2 | sort | uniq | wc -l
139904
# total genes
cat Bla_annot_final.gtf | cut -d ';' -f1 | cut -f9 | sort | uniq | wc -l
90927
# genes only from cufflinks
cat Bla_annot_final.gtf | grep cuf | cut -d ';' -f1 | cut -f9 | sort | uniq | wc -l
80382
```

remove fpkm=0: can't really do that because they don't report it in their assembly file!

## Some comments?

- why are there more genes in our cufflinks assembly (with gtf) than what's reported in the official annotation? check IDS in our assembly which are not present in the annot by cufflinks: 

```
# print list of gene IDs from cufflinks in the official annotation
cat Bla_annot_final.gtf | grep cuf | cut -d ';' -f2 | sort | uniq > official_genes_cuf.txt
# print list of gene IDs from gtf-guided cufflinks assembly
cat transcripts.gtf | cut -d ';' -f1 | cut -f9 | sort | uniq > cufflinks-gtf_genes.txt
# check what's in the cufflinks-gtf and not in the official annotation (cufflinks only)
comm -13 ../../../../../official_genes_cuf.txt cufflinks-gtf_genes.txt > genes_cufflinks-gtf_non-original-cufflinks.txt
# check where they come from 
grep -f genes_cufflinks-gtf_non-original-cufflinks.txt ../../../../../Bla_annot_final.gtf > genes_cufflinks-gtf_non-original-cufflinks_origin.txt
# check for evidence
cat genes_cufflinks-gtf_non-original-cufflinks_origin.txt | grep '_cuf' # no output!
```

all these come from evm evidence! -> so no problem. 
we should then either (1) use the complete annotation as reference, or (2) re-run assembly with just annotation from cufflinks as guide.

- then, are all transcripts assembled with cufflinks-gtf with fpkm=0 based on the annotation with evm?

```
cat transcripts.gtf | awk '$3=="transcript"{print $0}' | grep 'FPKM "0.0000000000"' | wc -l
159968
cat transcripts.gtf | awk '$3=="transcript"{print $0}' | grep '_evm' | wc -l
78162
```

the answer is no. 

- then, are there any transcripts assembled with cufflinks-gtf based on evm with fpkm>0? 

```
cat transcripts_nozero.gtf | awk '$3=="transcript"{print $0}' | grep '_evm' | wc -l
9710
```

the answer is yes. 

so I really think they put in their annotation some transcripts with fpkm=0. or, since our RNA dataset only includes material from embryo stages, then some of the genes annotated in the adult stage are not here. 
