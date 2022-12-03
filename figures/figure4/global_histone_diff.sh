#1) 对partially/fully up/down-regulated refumrs的甲基化水平可视化 (只选择了TSS近端2kb的)
#2) 对partially/fully up/down-regulated genes的TSS上下游3kb进行组蛋白修饰的可视化

#一、refumrs甲基化水平可视化
#1. refumrs区域准备
cd E:\IDH_Hyper\Results\4_Gene_Expression\202208_hyper_DEG_TSS\hyper_bed
#*_hyper.txt 
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/global_hist_diff/regions_TSS
#带正负链的
#把TSS近端2000bp的筛选出来:
awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){ print $1,$2,$3, "partial_down_umr"NR, 0, $23 } }' partial_down_hyper.txt >partial_down_refumr.bed
sed -i s/hyper_in_right/-/g partial_down_refumr.bed
sed -i s/hyper_in_left/+/g partial_down_refumr.bed

awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){ print $1,$2,$3, "partial_up_umr"NR, 0, $23 } }' partial_up_hyper.txt >partial_up_refumr.bed
sed -i s/hyper_in_right/-/g partial_up_refumr.bed
sed -i s/hyper_in_left/+/g partial_up_refumr.bed

awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){ print $1,$2,$3 } }' full_down_hyper.txt >full_down_refumr.bed

#   68   204  1620 full_down_refumr.bed
#  233  1398 11061 partial_down_refumr.bed
#  134   804  6046 partial_up_refumr.bed

#2. 可视化
#每组各三个区域
#dp_wig准备
#/share/pub/wangxy/IDH_Glioma/1_Identificaiton_IDH_Hyper/202110_PCA/Data/Brain_Normal_wig_list.txt
cd /share/pub/wangxy/IDH_Glioma/WGBS_analysis/deeptools_wig
wigToDeeptools -i /share/pub/wangxy/dlj/IDH/ref_umr2/brain_1_methy.wig -o brain_1_methy.dp.wig
wigToDeeptools -i /share/pub/wangxy/dlj/IDH/ref_umr2/brain_2_methy.wig -o brain_2_methy.dp.wig
wigToDeeptools -i /share/pub/wangxy/dlj/IDH/ref_umr2/brain_3_methy.wig -o brain_3_methy.dp.wig
wigToDeeptools -i /share/pub/wangxy/dlj/IDH/IDH-GBM1.wig -o IDH-GBM1.dp.wig
wigToDeeptools -i /share/pub/wangxy/dlj/IDH/IDH-GBM2.wig -o IDH-GBM2.dp.wig
wigToDeeptools -i /share/pub/wangxy/dlj/IDH/IDH-GBM3.wig -o IDH-GBM3.dp.wig
wigToBigWig brain_1_methy.dp.wig /share/pub/wangxy/software/ucsc_tools/hg19.chrom.sizes brain_1_methy.dp.bw
wigToBigWig brain_2_methy.dp.wig /share/pub/wangxy/software/ucsc_tools/hg19.chrom.sizes brain_2_methy.dp.bw
wigToBigWig brain_3_methy.dp.wig /share/pub/wangxy/software/ucsc_tools/hg19.chrom.sizes brain_3_methy.dp.bw
wigToBigWig IDH-GBM1.dp.wig /share/pub/wangxy/software/ucsc_tools/hg19.chrom.sizes IDH_GBM1_methy.dp.bw
wigToBigWig IDH-GBM2.dp.wig /share/pub/wangxy/software/ucsc_tools/hg19.chrom.sizes IDH_GBM2_methy.dp.bw
wigToBigWig IDH-GBM3.dp.wig /share/pub/wangxy/software/ucsc_tools/hg19.chrom.sizes IDH_GBM3_methy.dp.bw

#deeptools可视化
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/global_hist_diff
wig_dp_path=/share/pub/wangxy/IDH_Glioma/WGBS_analysis/deeptools_wig
computeMatrix scale-regions -p 20 -a 0 -b 0 --regionBodyLength 2000 -R ./regions_TSS/partial_up_refumr.bed ./regions_TSS/partial_down_refumr.bed ./regions_TSS/full_down_refumr.bed -S $wig_dp_path/brain_1_methy.dp.bw $wig_dp_path/brain_2_methy.dp.bw $wig_dp_path/brain_3_methy.dp.bw $wig_dp_path/IDH_GBM1_methy.dp.bw $wig_dp_path/IDH_GBM2_methy.dp.bw $wig_dp_path/IDH_GBM3_methy.dp.bw -o three_refumr_matrix.gz --outFileSortedRegions three_refumr_matrix.bed
plotHeatmap -m three_refumr_matrix.gz --startLabel start --endLabel end -o three_refumr_matrix_v1.pdf --plotFileFormat pdf --colorMap OrRd --yMin 0 --yMax 1 --outFileSortedRegions three_refumr_matrix_v1.bed
#plotHeatmap -m three_refumr_matrix.gz --startLabel start --endLabel end -o three_refumr_matrix_v2.pdf --plotFileFormat pdf --colorMap OrRd --outFileSortedRegions three_refumr_matrix_v2.bed

#对应genes
#deeptools排完序后再取TSS
#perl matrix.bed genes.tab >TSS.bed

#二、组蛋白修饰可视化
#v1
#准备TSS bed? 是否真的应该准备TSS bed? - 暂时先不准备
#根据排序后 - 直接取refumr bed 上下游1kb
#样本: 149, 150; AK076, AK124, AK213
Roadmap_path=/share/pub/wangxy/Project_Data/Roadmap_EGG/ChIP-seq/Brain_Hippocampus_Middle
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
for((i=0;i<${#histone[@]};i++))
do
brain1_file=$(ls $Roadmap_path/149_${histone[i]}*.bw)
brain2_file=$(ls $Roadmap_path/150_${histone[i]}*.bw)
IDH1_file=$(ls $GBM_path/*AK124_${histone[i]}*)
IDH2_file=$(ls $GBM_path/*AK213_${histone[i]}*)

computeMatrix scale-regions -p 20 -a 2000 -b 2000 --regionBodyLength 2000 -R ./regions_TSS/partial_up_refumr.bed ./regions_TSS/partial_down_refumr.bed ./regions_TSS/full_down_refumr.bed -S $brain1_file $brain2_file $IDH1_file $IDH2_file -o three_refumr_matrix_${histone[i]}.gz --outFileSortedRegions three_refumr_matrix_${histone[i]}.bed
plotHeatmap -m three_refumr_matrix_${histone[i]}.gz --startLabel start --endLabel end -o three_refumr_matrix_${histone[i]}_v1.pdf --plotFileFormat pdf --colorMap OrRd --outFileSortedRegions three_refumr_matrix_${histone[i]}_v1.bed
done

#v2
#标准化
#以CGIs计算信号值, 然后根据CGIs信号和 对其进行标准化
Roadmap_path=/share/pub/wangxy/Project_Data/Roadmap_EGG/ChIP-seq/Brain_Hippocampus_Middle
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq

cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/global_hist_diff/CGI_signal
awk '{OFS="\t"; print $1,$2,$3,"CGI"NR}' /share/pub/wangxy/Annotation/CGI/hg19_CGI.bed >hg19_CGI_bed4.bed

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
for((i=0;i<${#histone[@]};i++))
do
brain1_file=$(ls $Roadmap_path/149_${histone[i]}*.bw)
brain2_file=$(ls $Roadmap_path/150_${histone[i]}*.bw)
IDH1_file=$(ls $GBM_path/*AK124_${histone[i]}*)
IDH2_file=$(ls $GBM_path/*AK213_${histone[i]}*)
IDH3_file=$(ls $GBM_path/*AK076_${histone[i]}*)

bigWigAverageOverBed -bedOut=brain1_${histone[i]}.bed $brain1_file ./hg19_CGI_bed4.bed brain1_${histone[i]}.tab
bigWigAverageOverBed -bedOut=brain2_${histone[i]}.bed $brain2_file ./hg19_CGI_bed4.bed brain2_${histone[i]}.tab
bigWigAverageOverBed -bedOut=IDH1_${histone[i]}.bed $IDH1_file ./hg19_CGI_bed4.bed IDH1_${histone[i]}.tab
bigWigAverageOverBed -bedOut=IDH2_${histone[i]}.bed $IDH2_file ./hg19_CGI_bed4.bed IDH2_${histone[i]}.tab
bigWigAverageOverBed -bedOut=IDH3_${histone[i]}.bed $IDH3_file ./hg19_CGI_bed4.bed IDH3_${histone[i]}.tab
done

#三、组蛋白修饰定量 - boxplot
Roadmap_path=/share/pub/wangxy/Project_Data/Roadmap_EGG/ChIP-seq/Brain_Hippocampus_Middle
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq

#准备区域:
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/global_hist_diff/regions_TSS
awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){print $1,$10,$11,"p_up_hyper"NR} }' partial_up_hyper.txt >partial_up_hyper_bed4.bed
awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){print $1,$10,$11,"p_down_hyper"NR} }' partial_down_hyper.txt >partial_down_hyper_bed4.bed
awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){print $1,$10,$11,"f_down_hyper"NR} }' full_down_hyper.txt >full_down_hyper_bed4.bed

#这里应该用hyper还是用refumr好? 感觉用hyper合适点: 因为这部分区域发生了组蛋白修饰的变化
#这里的更名：(暂时不用这个:)
mv region_signal bak_hyper_region_signal
#
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/global_hist_diff/region_signal
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
hyper=(partial_up_hyper partial_down_hyper full_down_hyper)
for((i=0;i<${#histone[@]};i++))
do
	for((j=0;j<${#hyper[@]};j++))
	do
	brain1_file=$(ls $Roadmap_path/149_${histone[i]}*.bw)
	brain2_file=$(ls $Roadmap_path/150_${histone[i]}*.bw)
	IDH1_file=$(ls $GBM_path/*AK124_${histone[i]}*)
	IDH2_file=$(ls $GBM_path/*AK213_${histone[i]}*)
	IDH3_file=$(ls $GBM_path/*AK076_${histone[i]}*)
	bigWigAverageOverBed -bedOut=brain1_${histone[i]}_${hyper[j]}.bed $brain1_file ../regions_TSS/${hyper[j]}_bed4.bed brain1_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=brain2_${histone[i]}_${hyper[j]}.bed $brain2_file ../regions_TSS/${hyper[j]}_bed4.bed brain2_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=IDH1_${histone[i]}_${hyper[j]}.bed $IDH1_file ../regions_TSS/${hyper[j]}_bed4.bed IDH1_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=IDH2_${histone[i]}_${hyper[j]}.bed $IDH2_file ../regions_TSS/${hyper[j]}_bed4.bed IDH2_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=IDH3_${histone[i]}_${hyper[j]}.bed $IDH3_file ../regions_TSS/${hyper[j]}_bed4.bed IDH3_${histone[i]}_${hyper[j]}.tab
	done
done

#四、组蛋白修饰定量 - TSS boxplot
Roadmap_path=/share/pub/wangxy/Project_Data/Roadmap_EGG/ChIP-seq/Brain_Hippocampus_Middle
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq

#准备区域:
#genes 
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/global_hist_diff/regions_TSS
awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){print $6} }' partial_up_hyper.txt >partial_up_genes.txt
awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){print $6} }' partial_down_hyper.txt >partial_down_genes.txt
awk '{OFS="\t"; if( ($5 <= 2000 && $5 >= -2000) || ($5 == "NA") ){print $6} }' full_down_hyper.txt >full_down_genes.txt
#给partial_up_genes.txt添加部分我check过的K27me3 loss genes
CCND1
RUNX1
PDGFRA
SIX1
HOXA13
ZEB1
MYCL
SYK
CXCR4
LIMD1
AJUBA
ACKR3
TET1
sort partial_up_genes.txt | uniq >new_partial_up_genes.txt
perl get_TSS_promoter.pl new_partial_up_genes.txt /share/pub/wangxy/Project_Data/GTEx/GTEx_v19_gene_tab.txt >partial_up_promoter.bed
perl get_TSS_promoter.pl partial_down_genes.txt /share/pub/wangxy/Project_Data/GTEx/GTEx_v19_gene_tab.txt >partial_down_promoter.bed
perl get_TSS_promoter.pl full_down_genes.txt /share/pub/wangxy/Project_Data/GTEx/GTEx_v19_gene_tab.txt >full_down_promoter.bed

sort -k 1,1 -k 2,2n partial_up_promoter.bed >sort_partial_up_promoter.bed
sort -k 1,1 -k 2,2n partial_down_promoter.bed >sort_partial_down_promoter.bed
sort -k 1,1 -k 2,2n full_down_promoter.bed >sort_full_down_promoter.bed

awk '{OFS="\t"; print $1,$2,$3,"p_up"NR}' sort_partial_up_promoter.bed >partial_up_promoter_bed4.bed
awk '{OFS="\t"; print $1,$2,$3,"p_down"NR}' sort_partial_down_promoter.bed >partial_down_promoter_bed4.bed
awk '{OFS="\t"; print $1,$2,$3,"f_down"NR}' sort_full_down_promoter.bed >full_down_promoter_bed4.bed

rm sort_*
rm *_promoter.bed

#组蛋白修饰定量
mkdir promoter_signal
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/global_hist_diff/promoter_signal
Roadmap_path=/share/pub/wangxy/Project_Data/Roadmap_EGG/ChIP-seq/Brain_Hippocampus_Middle
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
hyper=(partial_up_promoter partial_down_promoter full_down_promoter)
for((i=0;i<${#histone[@]};i++))
do
	for((j=0;j<${#hyper[@]};j++))
	do
	brain1_file=$(ls $Roadmap_path/149_${histone[i]}*.bw)
	brain2_file=$(ls $Roadmap_path/150_${histone[i]}*.bw)
	IDH1_file=$(ls $GBM_path/*AK124_${histone[i]}*)
	IDH2_file=$(ls $GBM_path/*AK213_${histone[i]}*)
	IDH3_file=$(ls $GBM_path/*AK076_${histone[i]}*)
	bigWigAverageOverBed -bedOut=brain1_${histone[i]}_${hyper[j]}.bed $brain1_file ../regions_TSS/${hyper[j]}_bed4.bed brain1_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=brain2_${histone[i]}_${hyper[j]}.bed $brain2_file ../regions_TSS/${hyper[j]}_bed4.bed brain2_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=IDH1_${histone[i]}_${hyper[j]}.bed $IDH1_file ../regions_TSS/${hyper[j]}_bed4.bed IDH1_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=IDH2_${histone[i]}_${hyper[j]}.bed $IDH2_file ../regions_TSS/${hyper[j]}_bed4.bed IDH2_${histone[i]}_${hyper[j]}.tab
	bigWigAverageOverBed -bedOut=IDH3_${histone[i]}_${hyper[j]}.bed $IDH3_file ../regions_TSS/${hyper[j]}_bed4.bed IDH3_${histone[i]}_${hyper[j]}.tab
	done
done

