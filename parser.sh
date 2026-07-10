list=$(ls chunk_*)
for i in $list
do
	#sort $i > tmp
	join -t $'\t' -1 1 -2 1 $i ../NMPFam2_table.tsv > $i"_subset.tsv"
done
