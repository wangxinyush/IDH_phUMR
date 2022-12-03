cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/oncogene_histone/bed
awk '{OFS="\t"; print $1,$2,$3,"bed"NR}' hyper.bed >hyper_4col.bed
awk '{OFS="\t"; print $1,$2,$3,"bed"NR}' umr.bed >umr_4col.bed
awk '{OFS="\t"; print $1,$2,$3,"bed"NR}' downstream.bed >downstream_4col.bed
sed -i s/\\r//g hyper_4col.bed
sed -i s/\\r//g umr_4col.bed
sed -i s/\\r//g downstream_4col.bed

#定量 - b149 vs. g076
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/oncogene_histone/histone

bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
bed=(hyper umr downstream)
for((i=0;i<${#histone[@]};i++))
do
	for((j=0;j<${#bed[@]};j++))
	do
	b149_file=$bw_path/norm_b149_${histone[i]}_adjusted.wig.bw
	g076_file=$bw_path/norm_g076_${histone[i]}_adjusted.wig.bw

	bigWigAverageOverBed -bedOut=b149_${histone[i]}_${bed[j]}.bed $b149_file ../bed/${bed[j]}_4col.bed b149_${histone[i]}_${bed[j]}.tab
	bigWigAverageOverBed -bedOut=g076_${histone[i]}_${bed[j]}.bed $g076_file ../bed/${bed[j]}_4col.bed g076_${histone[i]}_${bed[j]}.tab
	done
done

#定量 - b150 vs. g213
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/oncogene_histone/histone

bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
bed=(hyper umr downstream)
for((i=0;i<${#histone[@]};i++))
do
	for((j=0;j<${#bed[@]};j++))
	do
	b150_file=$bw_path/norm_b150_${histone[i]}_adjusted.wig.bw
	g213_file=$bw_path/norm_g213_${histone[i]}_adjusted.wig.bw

	bigWigAverageOverBed -bedOut=b150_${histone[i]}_${bed[j]}.bed $b150_file ../bed/${bed[j]}_4col.bed b150_${histone[i]}_${bed[j]}.tab
	bigWigAverageOverBed -bedOut=g213_${histone[i]}_${bed[j]}.bed $g213_file ../bed/${bed[j]}_4col.bed g213_${histone[i]}_${bed[j]}.tab
	done
done
