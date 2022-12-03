#calucaltion signals for partially Hyper和flanking UMRs

#1) 准备输入数据 - partially/fully refumrs
cut -f 1-3 partial_IDH_Hyper_refumr.bed >partial_Hyper_refumr.bed
cut -f 1-3 partial_Common_Hyper_refumr.bed >>partial_Hyper_refumr.bed
sort -k 1,1 -k 2,2n partial_Hyper_refumr.bed | uniq >sort_partial_Hyper_refumr.bed
rm partial_Hyper_refumr.bed

cut -f 1-3 full_IDH_Hyper_refumr.bed >full_Hyper_refumr.bed
cut -f 1-3 full_Common_Hyper_refumr.bed >>full_Hyper_refumr.bed
sort -k 1,1 -k 2,2n full_Hyper_refumr.bed | uniq >sort_full_Hyper_refumr.bed
rm full_Hyper_refumr.bed

awk '{OFS="\t"; print $1,$2,$3,NR}' sort_partial_Hyper_refumr.bed >sort_partial_Hyper_refumr_bed4.bed
awk '{OFS="\t"; print $1,$2,$3,NR}' sort_full_Hyper_refumr.bed >sort_full_Hyper_refumr_bed4.bed

#3276   9828  78361 sort_full_Hyper_refumr.bed
#2373   7119  56725 sort_partial_Hyper_refumr.bed
 
#2) 计算信号均值
#Normal
cd /share/pub/wangxy/IDH_Glioma/2_Annotation_IDH_Hyper/202207_partial_full_refumrs/signal
Roadmap_path=/share/pub/wangxy/Project_Data/Roadmap_EGG/ChIP-seq/Brain_Hippocampus_Middle
samples=(149_H3K27ac 149_H3K27me3 149_H3K36me3 149_H3K4me1 149_H3K4me3 149_H3K9me3)
for((i=0;i<${#samples[@]};i++))
do
bigWigAverageOverBed -bedOut=partial_${samples[i]}.bed $Roadmap_path/${samples[i]}_sf_v1.bw ../sort_partial_Hyper_refumr_bed4.bed partial_hyper_${samples[i]}.tab
bigWigAverageOverBed -bedOut=full_${samples[i]}.bed $Roadmap_path/${samples[i]}_sf_v1.bw ../sort_full_Hyper_refumr_bed4.bed full_hyper_${samples[i]}.tab
done

samples=(150_H3K27ac 150_H3K27me3 150_H3K36me3 150_H3K4me1 150_H3K4me3 150_H3K9me3)
for((i=0;i<${#samples[@]};i++))
do
bigWigAverageOverBed -bedOut=partial_${samples[i]}.bed $Roadmap_path/${samples[i]}_sf_v1.bw ../sort_partial_Hyper_refumr_bed4.bed partial_hyper_${samples[i]}.tab
bigWigAverageOverBed -bedOut=full_${samples[i]}.bed $Roadmap_path/${samples[i]}_sf_v1.bw ../sort_full_Hyper_refumr_bed4.bed full_hyper_${samples[i]}.tab
done

#IDH
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
samples=(AK066 AK076 AK124 AK199 AK213 AK231)
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq
for((i=0;i<${#samples[@]};i++))
do
	for((j=0;j<${#histone[@]};j++))
	do
	file=$(ls $GBM_path/*${samples[i]}_${histone[j]}*)
	bigWigAverageOverBed -bedOut=full_hyper_${samples[i]}_${histone[j]}.bed $file ../sort_partial_Hyper_refumr_bed4.bed partial_hyper_${samples[i]}_${histone[j]}.tab
	bigWigAverageOverBed -bedOut=partial_hyper_${samples[i]}_${histone[j]}.bed $file ../sort_full_Hyper_refumr_bed4.bed full_hyper_${samples[i]}_${histone[j]}.tab
	done
done

#IDH_v2
cd /share/pub/wangxy/IDH_Glioma/2_Annotation_IDH_Hyper/202207_partial_full_refumrs/signal_IDH_v2
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
samples=(AK076 AK124)
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq
for((i=0;i<${#samples[@]};i++))
do
	for((j=0;j<${#histone[@]};j++))
	do
	file=$(ls $GBM_path/*${samples[i]}_${histone[j]}*)
	echo $file
	bigWigAverageOverBed -bedOut=full_hyper_${samples[i]}_${histone[j]}.bed $file ../sort_partial_Hyper_refumr_bed4.bed partial_hyper_${samples[i]}_${histone[j]}.tab
	bigWigAverageOverBed -bedOut=partial_hyper_${samples[i]}_${histone[j]}.bed $file ../sort_full_Hyper_refumr_bed4.bed full_hyper_${samples[i]}_${histone[j]}.tab
	done
done

