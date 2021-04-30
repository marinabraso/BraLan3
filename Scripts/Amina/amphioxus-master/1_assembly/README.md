
# Assembly

Several tests were performent, only relevant assemblies are reported here.


### Canu with in-build error correction

Imput is set of all raw subreads. The datasets contain seqeunces from 83 SMRT cells, 163x in total [link to SRA].
Error rates 0.035 and the default were tested. (add version)


```sh
canu -assemble \
  -p amphio -d amphio-assembly_er025 \
  genomeSize=520m \
  errorRate=0.025 \
  -pacbio-corrected amphio.trimmedReads.fastq \
  -useGrid=false \
  -maxMemory=256 \
  -maxThreads=32
```

I ve created assembly of total 951448017 bases in 5505 contigs (about 10k contigs were discarded by Canu, then the sum would be about twice the genome length). The longest contig has 6.4Mpb, N50 is 706164, average length is 172833.


### Falcon with in-built error correction

```
module add UHTS/PacBio/FALCON-integrate/1.8.6
```

Falcon takes particulary named seqeunces in fasta format which can be prepared using [FALCON-formater](https://github.com/zyndagj/FALCON-formatter) :

```sh
FALCON-formatter amphix_complete.fastq
```

Every SMRT cell was saved in separated fasta.
Falcon request a list of those files.

```sh
mkdir -p data/falcon_asm
ls <full_path>/cells/ | awk '{print "<full_path>/cells/" $0}' > data/falcon_asm/input.fofn
# ls data/raw_reads/reads_pb/*fasta | awk '{print "/scratch/local/kjaron/amphioxus/" $0}' > input.fofn
```

finaly using `input.fofn` and `fc_run_init.cfg` I run falcon (add version, still runing)

```sh
cd data/falcon_asm/
fc_run fc_run_init.cfg
```

- I used a cutoff for filtering reads calculated to keep 50x of reads

#### Update

It crushed, question why. Log referees to three failed jobs :

```
2017-12-13 04:41:13,642 - pypeflow.controller - DEBUG - task status: task://localhost/d_000b_raw_reads, 'fail'
2017-12-13 04:41:13,683 - pypeflow.controller - DEBUG - task status: task://localhost/d_01fd_raw_reads, 'fail'
2017-12-13 04:41:13,798 - pypeflow.controller - DEBUG - task status: task://localhost/d_0755_raw_reads, 'fail'
```

The first attempt for solution is simply erase failed jobs, erase watcher dir and restart the program (solution inspired by this [thread](https://github.com/PacificBiosciences/FALCON/issues/392)).

```
rm -r 0-rawreads/job_0{00b,1fd,755}
rm -rf mypwatcher
```

I add `~/bin` as primary location of binaries (at beginning of PATH), so my updated version of DEALIGNER I just installed will be used

```
export PATH=~/bin:$PATH
```

and restart...

```
fc_run /scratch/beegfs/monthly/kjaron/amphioxus/1_assembly/fc_run_init.cfg
```

#### harakiri of daligner versions

I had to rename files I generated with newer version of daligner :

```
kjaron@dee-serv04:/scratch/local/kjaron/amphioxus/data/falcon_asm/0-rawreads/job_0755$ for i in 24 25 26; do mv raw_reads."$i".raw_reads.77.las L1."$i".77.las; done
kjaron@dee-serv04:/scratch/local/kjaron/amphioxus/data/falcon_asm/0-rawreads/job_0755$ for i in 24 25 26; do mv raw_reads.77.raw_reads."$i".las L1.77."$i".las; done
```

I done it for all three failed jobs and touched "job done".

I also put back the old daligner, since apparently this version of falcon is used to it.
