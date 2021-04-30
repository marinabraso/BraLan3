# Merging with Illumina assembly

...


### input data

- reference haplotype of PacBio (PB) assembly by Canu: amphio_PB_ref_notpolished.fa
- reference haplotype of Illumina (IL) assembly: Bl71nemr.fa

### merging

the easiest would be to treat those assemblies as independently assembled haplotypes. Therefore i will just cat those assemblies togheder and do a second round of HaploMerger2.

To keep a track of the origin of scaffolds in the concatinated file i will add `IL_` and `PB_` before headers. I will also erase all spaces, : and ;.

```
cat Bl71nemr.fa | tr -s ' ' | tr ' ' '_' | sed 's/>/>IL_/' | tr ':' '_' | tr ';' '_' > amphio_PBIL.fa
cat amphio_PB_ref_notpolished.fa | tr -s ' ' | tr ' ' '_' | sed 's/>/>PB/' | tr ':' '_' | tr ';' '_' >> amphio_PBIL.fa
gzip amphio_PBIL.fa
```

then all HaploMerger2 script are executed

### Using only mp to scaffold

#### mapping

using BWA and script provided together with [BESST] and script I wrote for assembly of timema (LINK TO TIMEMA REPO).

```
GENOME=/scratch/beegfs/monthly/kjaron/amphioxus/data/asm_pb/amphio_polished.fa
MAP=~/timema_assembly/G_scaffolding/1_mapp_mp_lsf.sh
READS=/scratch/beegfs/monthly/kjaron/amphioxus/data/raw_reads/reads_il/

for i in 3 5 8; do
        $MAP  BOSC_mp_is$i""k $GENOME $READS""ATQ_BOSC_IND$i""_1.fq.gz $READS""ATQ_BOSC_IND$i""_2.fq.gz
done;

for i in 2 5 6; do
        $MAP BOSN_mp_is$i""k $GENOME $READS""ATQ_BOSN_IND$i""_1.fq.gz $READS""ATQ_BOSN_IND$i""_2.fq.gz
done;
```

#### scaffolding

using BESST (and interactive job on Vital-it)

```bash
python ~/src/BESST/runBESST -c amphio_polished.fa -f BOSN_mp_is2k.bam BOSC_mp_is3k.bam BOSC_mp_is5k.bam BOSN_mp_is5k.bam BOSN_mp_is6k.bam BOSC_mp_is8k.bam -o amphio_Canu_hm_BESST --no_score
```

#### Gapfilling

using PBjelly

load on Vital-it:

```bash
source ~/scripts/load_PBsuite.sh
```

create a setting file `amphio_PBjelly_protocol.xml`. 

(TO BE ADDED)


### visualisation

Using minimap the relation of original PB assembly was plotteg against merged assembly

[[amphio_PB2merged.pdf]]

### alternative approach

I would like to keep the information in PacBio assembly, therefore I might consider to use mate pairs to scaffold.


### In both cases there are gaps

Gap Filling will be performed using PBjelly. Load PB jelly

```
source ~/scripts/load_PBsuite.sh
```
