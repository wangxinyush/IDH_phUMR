#149, 150, AK076, AK213 甲基化水平(deeptools)
#wigToDeeptools -i ../GSM3444668_AK213.wig -o GSM3444668_AK213_dp.wig
#wigToBigWig GSM3444668_AK213_dp.wig /share/pub/wangxy/software/genome/ucsc/hg19/hg19.chrom.sizes GSM3444668_AK213_dp.bw

methy149_bw=/share/pub/wangxy/IDH_Glioma/WGBS_analysis/deeptools_wig/brain_3_methy.dp.bw
methy150_bw=/share/pub/wangxy/IDH_Glioma/WGBS_analysis/deeptools_wig/brain_2_methy.dp.bw
methy076_bw=/share/pub/wangxy/Project_Data/60_GBM_Data/WGBS/deeptools_bw/GSM3444634_AK076_dp.bw
methy213_bw=/share/pub/wangxy/Project_Data/60_GBM_Data/WGBS/deeptools_bw/GSM3444668_AK213_dp.bw

#基因bed
gene_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/gene_bed

#2) 可视化
#TSS - TES
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/methy_compare/paired_TSS_methy
outname=b149_g076_TSS_TES
computeMatrix scale-regions -p 20 \
 -R $gene_path/b149_g076_partial_up_gene.bed $gene_path/b149_g076_partial_down_gene.bed $gene_path/b149_g076_full_up_gene.bed $gene_path/b149_g076_full_down_gene.bed \
 -S $methy149_bw $methy076_bw \
 -a 3000 -b 3000 --regionBodyLength 9000 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax 1 \
 --missingDataColor \#DFDBD9 \
 --colorMap coolwarm --plotFileFormat pdf --perGroup

outname=b150_g213_TSS_TES
computeMatrix scale-regions -p 20 \
 -R $gene_path/b150_g213_partial_up_gene.bed $gene_path/b150_g213_partial_down_gene.bed $gene_path/b150_g213_full_up_gene.bed $gene_path/b150_g213_full_down_gene.bed \
 -S $methy150_bw $methy213_bw \
 -a 3000 -b 3000 --regionBodyLength 9000 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax 1 \
 --missingDataColor \#DFDBD9 \
 --colorMap coolwarm --plotFileFormat pdf --perGroup

#修改颜色
outname=b149_g076_TSS_TES
plotHeatmap -m $outname.gz -o $outname\_bwrr.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax 1 \
 --missingDataColor \#F2CBB7 \
 --colorList '#3C4DC3,#7A9DF8,#BFD3F6,#F2CBB7,#EE8468,#B80427' \
 --plotFileFormat pdf --perGroup

outname=b150_g213_TSS_TES
plotHeatmap -m $outname.gz -o $outname\_bwrr.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax 1 \
 --missingDataColor \#F2CBB7 \
 --colorList '#3C4DC3,#7A9DF8,#BFD3F6,#F2CBB7,#EE8468,#B80427' \
 --plotFileFormat pdf --perGroup


#TSS up, down 3kb
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/methy_compare/paired_TSS_methy
outname=b149_g076_TSS_3kb
computeMatrix reference-point -p 20 \
 -R $gene_path/b149_g076_partial_up_gene.bed $gene_path/b149_g076_partial_down_gene.bed $gene_path/b149_g076_full_up_gene.bed $gene_path/b149_g076_full_down_gene.bed \
 -S $methy149_bw $methy076_bw \
 -a 3000 -b 3000 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax 1 \
 --colorMap coolwarm --plotFileFormat pdf --perGroup

outname=b150_g213_TSS_3kb
computeMatrix reference-point -p 20 \
 -R $gene_path/b150_g213_partial_up_gene.bed $gene_path/b150_g213_partial_down_gene.bed $gene_path/b150_g213_full_up_gene.bed $gene_path/b150_g213_full_down_gene.bed \
 -S $methy150_bw $methy213_bw \
 -a 3000 -b 3000 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax 1 \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
