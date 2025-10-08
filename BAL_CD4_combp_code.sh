# PvC
awk -F"," '{print $13 "\t" $14 "\t" $15 "\t" $8}' PvC_DifferentialMethylation_FullResults_peer_correctedsex.csv > test.txt
gsed -i 's/\"//g' test.txt
sort -k 1,1 -k2,2n test.txt > srt.bed
gsed -i '$d' srt.bed
echo "#chrom\tstart\tend\tpvalue" | cat - srt.bed > ./comb-p/PvC_input.bed


# NPvC
awk -F"," '{print $13 "\t" $14 "\t" $15 "\t" $8}' NPvC_DifferentialMethylation_FullResults_peer_correctedsex.csv > test.txt
gsed -i 's/\"//g' test.txt
sort -k 1,1 -k2,2n test.txt > srt.bed
gsed -i '$d' srt.bed
echo "#chrom\tstart\tend\tpvalue" | cat - srt.bed > ./comb-p/NPvC_input.bed

# PvNP
awk -F"," '{print $13 "\t" $14 "\t" $15 "\t" $8}' PvNP_DifferentialMethylation_FullResults_peer_correctedsex.csv > test.txt
gsed -i 's/\"//g' test.txt
sort -k 1,1 -k2,2n test.txt > srt.bed
gsed -i '$d' srt.bed
echo "#chrom\tstart\tend\tpvalue" | cat - srt.bed > ./comb-p/PvNP_input.bed

# SvC
awk -F"," '{print $13 "\t" $14 "\t" $15 "\t" $8 "\t" $17}' redo_SvC_DifferentialMethylation_FullResults_peer_correctedsex.csv > test.txt
gsed -i 's/\"//g' test.txt
sort -k 1,1 -k2,2n test.txt > srt.bed
gsed -i '$d' srt.bed
echo "#chrom\tstart\tend\tpvalue" | cat - srt.bed > ./comb-p/redo_SvC_input.bed

# run comb-p
comb-p pipeline --dist 300 -p PvC PvC_input.bed
comb-p pipeline --dist 300 -p NPvC NPvC_input.bed
comb-p pipeline --dist 300 -p SvC SvC_input.bed
comb-p pipeline --dist 300 -p PvNP PvNP_input.bed