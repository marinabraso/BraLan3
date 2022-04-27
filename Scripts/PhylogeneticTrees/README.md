

##### Scripts:
> Extract orthologous group sequences (AA and DNA) for Amphioxus and vertebrae species
- Scripts/PhylogeneticTrees/ExtractOrthologousGroupSequences.sh
> MSA of the AA sequences with mafft + Tree with RAxML
- Scripts/PhylogeneticTrees/MSA_TreeReconstruction.sh
- Scripts/PhylogeneticTrees/Run_MSA_TreeReconstruction.sh

##### Usage:
```
sbatch -t 05:00:00 --mem=8000 -J ExtrAmpVert -o tmp/ExtractSeq_AmphVerteb.out -e tmp/ExtractSeq_AmphVerteb.err Scripts/dNdS/ExtractOrthologousGroupSequences.sh AmphVerteb

./Scripts/PhylogeneticTrees/Run_MSA_TreeReconstruction.sh AmphVerteb

```
