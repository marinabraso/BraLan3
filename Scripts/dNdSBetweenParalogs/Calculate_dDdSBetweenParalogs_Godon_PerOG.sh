#!/bin/bash

# Scripts
pwd=$(pwd)
godon=${pwd}/../Software/godon
RAxML=${pwd}/../Software/standard-RAxML-master/raxmlHPC


group=$1
ResultsFolder=$2

cd ${ResultsFolder}/TreeReconstruction
${RAxML} -s ../MSA_cleaning_TCoffee/${group}_DNA_clean.fa -n ${group}_GTRGAMMA -m GTRGAMMA -p 12345
cd ${pwd}
if [[ -s ${ResultsFolder}/TreeReconstruction/RAxML_result.${group}_GTRGAMMA ]]
then
	${godon} test M8 --m0-tree ${ResultsFolder}/MSA_cleaning_TCoffee/${group}_DNA_clean.fa ${ResultsFolder}/TreeReconstruction/RAxML_result.${group}_GTRGAMMA > ${ResultsFolder}/Godon_M8/${group}_M8_likelihood_ratio.txt
else
	echo ${group} >> ${ResultsFolder}/MissingOG_TooFewGenes.txt
fi


