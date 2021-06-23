
### dNdS between paralogs

##### Scripts:
> Extract orthologous group sequences (AA and DNA) for Amphioxus and vertebrae species
- Scripts/dNdS/ExtractOrthologousGroupSequences.sh
> MSA of the AA sequences with mafft + backtranslation to DNA (from the DNA sequences)
- Scripts/dNdS/MSA_AA_backtranslation_DNA.sh
- Scripts/dNdS/Run_MSA_AA_backtranslation_DNA.sh
> evaluate the AA MSA done with mafft t_coffe, extract score and filter mafft alignments with t_coffee score
- Scripts/dNdS/MSA_cleaning_tcoffe.sh
- Scripts/dNdS/Run_MSA_cleaning_tcoffe_PerOG.sh
> Tree from alignment with RAxML and dNdS analysis with Godon (model M8)
- Scripts/dNdS/Calculate_dDdSBetweenParalogs_Godon.sh
- Scripts/dNdS/Run_Calculate_dDdSBetweenParalogs_Godon.sh

##### Usage:
```
sbatch -t 05:00:00 --mem=8000 -J ExtrAmpVert -o tmp/ExtractSeq_AmphVerteb.out -e tmp/ExtractSeq_AmphVerteb.err Scripts/dNdS/ExtractOrthologousGroupSequences.sh AmphVerteb

./Scripts/dNdS/Run_MSA_AA_backtranslation_DNA.sh AmphVerteb

./Scripts/dNdS/Run_MSA_cleaning_tcoffe.sh AmphVerteb

./Scripts/dNdS/Run_Calculate_dDdSBetweenParalogs_Godon.sh AmphVerteb
```
