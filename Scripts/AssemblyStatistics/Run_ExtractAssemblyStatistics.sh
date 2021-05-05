#!/bin/bash




./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Branchiostoma_lanceolatum.BraLan3
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Asterias_rubens.eAstRub1.3
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Branchiostoma_belcheri.Haploidv18h27
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Branchiostoma_floridae.Bfl_VNyyK
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Danio_rerio.GRCz11
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Gallus_gallus.GRCg6a
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Homo_sapiens.GRCh38
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Mus_musculus.GRCm39
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Saccoglossus_kowalevskii.Skow1.1
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Strongylocentrotus_purpuratus.Spur5.0
./Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh Bl71nemr_Bralan2_fromAmphiencode


# Merge results
line="Statistic"
for f in Branchiostoma_lanceolatum.BraLan3 Bl71nemr_Bralan2_fromAmphiencode Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1 
do 
	line=$line"	"$f
done

echo $line > Results/AssemblyStatistics/Joined_Statistics.txt
for l in 1 2 3 4 5 6 7 8 9 10 11 12 # For each statistic
do
	line=$(awk -v l=${l} '{if(NR==l){print $1}}' Results/AssemblyStatistics/Branchiostoma_lanceolatum.BraLan3/Branchiostoma_lanceolatum.BraLan3_statistics.txt)
	for f in Branchiostoma_lanceolatum.BraLan3 Bl71nemr_Bralan2_fromAmphiencode Branchiostoma_belcheri.Haploidv18h27 Branchiostoma_floridae.Bfl_VNyyK Danio_rerio.GRCz11 Gallus_gallus.GRCg6a Homo_sapiens.GRCh38 Mus_musculus.GRCm39 Strongylocentrotus_purpuratus.Spur5.0 Asterias_rubens.eAstRub1.3 Saccoglossus_kowalevskii.Skow1.1 
	do
		line=$line"	"$(awk -v l=${l} '{if(NR==l){print $2}}' Results/AssemblyStatistics/${f}/${f}_statistics.txt)
	done
	echo $line >> Results/AssemblyStatistics/Joined_Statistics.txt
done

# Transpose the information
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' Results/AssemblyStatistics/Joined_Statistics.txt | sed 's/C:\([0-9\.]\+\)%\[S:\([0-9\.]\+\)%,D:\([0-9\.]\+\)%\].*$/\1\t\2\t\3/g' | sed 's/ /\t/g' > Results/AssemblyStatistics/Joined_Statistics_t.tmp	
awk '{if(NR==1){print $0"_C\tBUSCO_S\tBUSCO_D"}else{print $0}}' Results/AssemblyStatistics/Joined_Statistics_t.tmp | sed 's/proteinmetazoa//g' > Results/AssemblyStatistics/Joined_Statistics_t.txt 
rm Results/AssemblyStatistics/Joined_Statistics_t.tmp

