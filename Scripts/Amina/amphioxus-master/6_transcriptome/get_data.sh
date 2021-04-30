# get reference genome & annotation
wget https://www.dropbox.com/s/r9s9wd3osyf6zud/Bl71nemr.fa.gz; gzip -d Bl71nemr.fa.gz 
wget https://www.dropbox.com/s/d4fqnoa8gdix3pa/Bla_annot_final.gtf.gz; gzip -d Bla_annot_final.gtf.gz

# prepare gtf for bgee pipeline: 
# remove CDS
cat Bla_annot_final.gtf | awk '$3!="CDS"{print $0}' > Bla_annot_final_noCDS.gtf
# compute transcript features
module add UHTS/Assembler/cufflinks/2.2.1
gffread -E Bla_annot_final_noCDS.gtf -o- > Bla_annot_final_noCDS_features.gff3
python /home/aechchik/software/gff2gtf.py Bla_annot_final_noCDS_features.gff3 > Bla_annot_final_noCDS_features.gtf
# compute gene features 
cat Bla_annot_final_noCDS_features.gtf | awk '$3=="transcript"{print $0}' > transcripts.gtf # lucky enough, just one transcript per gene in the annotation
cat transcripts.gtf | sed 's/; transcript_id.*$/;/' | sed 's/transcript/gene/' > genes.gtf 
# merge all & sort by coordinate
cat Bla_annot_final_noCDS_features.gtf genes.gtf | sort -k4n > Bla_annot_final_allfeatures.gtf

# test before validation
