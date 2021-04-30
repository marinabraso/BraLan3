Description of files in this section. 

- `get_data.sh`: download reference genome & 'final' annotation (FR), modify annotation to include transcript and gene features (note: 1 transcript reported per gene locus - realistic?)
- `get_intergenic_regions.R`: inspired from bgee pipeline to prepare the annotation for the computation of intergenic regions
- `get_abundance.sh`: compute tpm of genic & intergenic regions using kallisto
- `amphio_genic-intergenic_plot.sh`: generate density plots from genic & intergenic regions from abundance
