#149 vs. AK076; 150 vs. AK213
##1. 识别partially/fully hyper
#1) 识别UMR
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/WGBS

#1.1) umr more than 200bp
awk '{OFS="\t"; if($5 >= 200){print $0} }' /share/pub/wangxy/IDH_Glioma/1_Identificaiton_IDH_Hyper/202111_new_refumrs/brain_refumrs/samples/brain_3/brain_3_UM.bed | cut -f 1-3 >brain149_mr200_UMR.bed
awk '{OFS="\t"; if($5 >= 200){print $0} }' /share/pub/wangxy/IDH_Glioma/1_Identificaiton_IDH_Hyper/202111_new_refumrs/brain_refumrs/samples/brain_1/brain_1_UM.bed | cut -f 1-3 >brain150_mr200_UMR.bed
sed -i '1d' brain149_mr200_UMR.bed #26665
sed -i '1d' brain150_mr200_UMR.bed #25509

#2) 识别Hyper
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/WGBS
ls /share/pub/wangxy/dlj/IDH/ref_umr2/brain_3_methy.wig >brain149_AK076_wl.txt
ls /share/pub/wangxy/dlj/IDH/ref_umr2/brain_1_methy.wig >brain150_AK213_wl.txt
ls /share/pub/wangxy/Project_Data/60_GBM_Data/WGBS/GSM3444634_AK076.wig >>brain149_AK076_wl.txt
ls /share/pub/wangxy/Project_Data/60_GBM_Data/WGBS/GSM3444668_AK213.wig >>brain150_AK213_wl.txt

perl get_umr_hyper.pl brain149_mr200_UMR.bed /share/pub/wangxy/Annotation/CG/hg19_CpG.txt -1 brain149_AK076_wl.txt 1,1 brain149_AK076_umr_hyper.bed
perl get_umr_hyper.pl brain150_mr200_UMR.bed /share/pub/wangxy/Annotation/CG/hg19_CpG.txt -1 brain150_AK213_wl.txt 1,1 brain150_AK213_umr_hyper.bed

#1.2) umr more than 200bp, mean methylation level < 0.1
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/WGBS/umr_200_mean0.1
awk '{OFS="\t"; if($5 >= 200 && $7 < 0.1){print $0} }' /share/pub/wangxy/IDH_Glioma/1_Identificaiton_IDH_Hyper/202111_new_refumrs/brain_refumrs/samples/brain_3/brain_3_UM.bed | cut -f 1-3 >brain149_mr200_mean0.1_UMR.bed #26578
awk '{OFS="\t"; if($5 >= 200 && $7 < 0.1){print $0} }' /share/pub/wangxy/IDH_Glioma/1_Identificaiton_IDH_Hyper/202111_new_refumrs/brain_refumrs/samples/brain_1/brain_1_UM.bed | cut -f 1-3 >brain150_mr200_mean0.1_UMR.bed #25453
cp ../umr_200/*_wl.txt ./
perl ../get_umr_hyper.pl brain149_mr200_mean0.1_UMR.bed /share/pub/wangxy/Annotation/CG/hg19_CpG.txt -1 brain149_AK076_wl.txt 1,1 brain149_AK076_umr_hyper.bed #10083
perl ../get_umr_hyper.pl brain150_mr200_mean0.1_UMR.bed /share/pub/wangxy/Annotation/CG/hg19_CpG.txt -1 brain150_AK213_wl.txt 1,1 brain150_AK213_umr_hyper.bed #9811

#3) Annotation
E:\IDH_Hyper\Results\4_Gene_Expression\202208_oncogenes_Track\pairs_hist_diff\get_umr_hyper

##2. 差异表达基因
E:\IDH_Hyper\Results\4_Gene_Expression\202208_oncogenes_Track\pairs_hist_diff\get_umr_hyper\region_anno.R

##3. 组蛋白修饰
#bw to wig
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/ChIP
Roadmap_path=/share/pub/wangxy/Project_Data/Roadmap_EGG/ChIP-seq/Brain_Hippocampus_Middle
GBM_path=/share/pub/wangxy/Project_Data/60_GBM_Data/ChIP-seq
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)

for((i=0;i<${#histone[@]};i++))
do
b149_bw=$(ls $Roadmap_path/149_${histone[i]}*.bw)
b150_bw=$(ls $Roadmap_path/150_${histone[i]}*.bw)
g076_bw=$(ls $GBM_path/*AK076_${histone[i]}*)
g213_bw=$(ls $GBM_path/*AK213_${histone[i]}*)

bigWigToWig $b149_bw b149_${histone[i]}.wig
bigWigToWig $b150_bw b150_${histone[i]}.wig
bigWigToWig $g076_bw g076_${histone[i]}.wig
bigWigToWig $g213_bw g213_${histone[i]}.wig
done

#4) 可视化
#plot by group, gene_bed_v2
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/deeptools_v3
bed_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/gene_bed_v2
bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
max=(15 10 10 10 30 10)
for((i=0;i<${#histone[@]};i++))
do
#partial - up/down
yMax=${max[i]}
outname=b149_g076_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b149_g076_partial_up_gene.bed $bed_path/b149_g076_partial_down_gene.bed $bed_path/b149_g076_full_up_gene.bed $bed_path/b149_g076_full_down_gene.bed \
 -S $bw_path/norm_b149_${histone[i]}_adjusted.wig.bw $bw_path/norm_g076_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 9000 \
 --scale 2 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
outname=b150_g213_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b150_g213_partial_up_gene.bed $bed_path/b150_g213_partial_down_gene.bed $bed_path/b150_g213_full_up_gene.bed $bed_path/b150_g213_full_down_gene.bed \
 -S $bw_path/norm_b150_${histone[i]}_adjusted.wig.bw $bw_path/norm_g213_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 9000 \
 --scale 2 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
done

#print H3K4me3
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
max=(15 10 10 10 35 10)
i=4
yMax=${max[i]}
outname=b149_g076_${histone[i]}
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
outname=b150_g213_${histone[i]}
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup

<<deeptools_V2
#plot by group
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/deeptools_v2
bed_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/gene_bed
bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
max=(15 10 10 10 25 10)
for((i=0;i<${#histone[@]};i++))
do
#partial - up/down
yMax=${max[i]}
outname=b149_g076_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b149_g076_partial_up_gene.bed $bed_path/b149_g076_partial_down_gene.bed $bed_path/b149_g076_full_up_gene.bed $bed_path/b149_g076_full_down_gene.bed \
 -S $bw_path/norm_b149_${histone[i]}_adjusted.wig.bw $bw_path/norm_g076_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 9000 \
 --scale 2 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
outname=b150_g213_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b150_g213_partial_up_gene.bed $bed_path/b150_g213_partial_down_gene.bed $bed_path/b150_g213_full_up_gene.bed $bed_path/b150_g213_full_down_gene.bed \
 -S $bw_path/norm_b150_${histone[i]}_adjusted.wig.bw $bw_path/norm_g213_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 9000 \
 --scale 2 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
done
deeptools_V2

<<deeeptools_v1
#deeptools - 可视化 - (bak__V1)
#plot by sample
#将基因列表生成deeptools需要的bed (chr start end strand)
#在R里就可以做了:
write_gene_tab( intersect(pair1_up_DEG, pair1_gene_list$partial_promoter), "gene_bed\\b149_g076_partial_up_gene.bed" )
write_gene_tab( intersect(pair1_down_DEG, pair1_gene_list$partial_promoter), "gene_bed\\b149_g076_partial_down_gene.bed" )
write_gene_tab( intersect(pair1_up_DEG, pair1_gene_list$full_promoter), "gene_bed\\b149_g076_full_up_gene.bed" )
write_gene_tab( intersect(pair1_down_DEG, pair1_gene_list$full_promoter), "gene_bed\\b149_g076_full_down_gene.bed" )

write_gene_tab <- function(symbols, outfile){
  #symbols <- intersect(pair1_up_DEG, pair1_gene_list$partial_promoter)
  #outfile <- "gene_bed\\b149_g076_partial_up_gene.bed"
  GTEx_gene_tab <- read.table("E:\\IDH_Hyper\\Data\\annotation\\GTEx_v19_gene_tab_ucsc.txt",
                              header=F, sep="\t", fill=T, check.names=F, quote="")
  sub_tab <- GTEx_gene_tab[match(symbols, GTEx_gene_tab$V7), ]
  sort_sub_tab <- sub_tab[order(sub_tab$V1, sub_tab$V2), ]
  write.table(sort_sub_tab, outfile, sep="\t", row.names=F, col.names=F, quote=F)
}

cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/deeptools
bed_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/gene_bed
bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP

histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
max=(15 10 10 10 25 10)
for((i=0;i<${#histone[@]};i++))
do
#partial - up/down
computeMatrix scale-regions -p 20 -R $bed_path/b149_g076_partial_up_gene.bed $bed_path/b149_g076_partial_down_gene.bed -S $bw_path/norm_b149_${histone[i]}_adjusted.wig.bw $bw_path/norm_g076_${histone[i]}_adjusted.wig.bw -a 3000 -b 3000 --regionBodyLength 9000 --scale 2 -o b149_g076_partial_${histone[i]}.gz --outFileNameMatrix b149_g076_partial_${histone[i]}.tab --outFileSortedRegions b149_g076_partial_${histone[i]}.bed
plotHeatmap -m b149_g076_partial_${histone[i]}.gz --colorMap coolwarm --startLabel start --endLabel end -o b149_g076_partial_${histone[i]}.pdf --plotFileFormat pdf --yMin 0 --yMax ${max[i]}
#full - up/down
computeMatrix scale-regions -p 20 -R $bed_path/b149_g076_full_up_gene.bed $bed_path/b149_g076_full_down_gene.bed -S $bw_path/norm_b149_${histone[i]}_adjusted.wig.bw $bw_path/norm_g076_${histone[i]}_adjusted.wig.bw -a 3000 -b 3000 --regionBodyLength 9000 --scale 2 -o b149_g076_full_${histone[i]}.gz --outFileNameMatrix b149_g076_full_${histone[i]}.tab --outFileSortedRegions b149_g076_full_${histone[i]}.bed
plotHeatmap -m b149_g076_full_${histone[i]}.gz --colorMap coolwarm --startLabel start --endLabel end -o b149_g076_full_${histone[i]}.pdf --plotFileFormat pdf --yMin 0 --yMax ${max[i]}

computeMatrix scale-regions -p 20 -R $bed_path/b150_g213_partial_up_gene.bed $bed_path/b150_g213_partial_down_gene.bed -S $bw_path/norm_b150_${histone[i]}_adjusted.wig.bw $bw_path/norm_g213_${histone[i]}_adjusted.wig.bw -a 3000 -b 3000 --regionBodyLength 9000 --scale 2 -o b150_g213_partial_${histone[i]}.gz --outFileNameMatrix b150_g213_partial_${histone[i]}.tab --outFileSortedRegions b150_g213_partial_${histone[i]}.bed
plotHeatmap -m b150_g213_partial_${histone[i]}.gz --colorMap coolwarm --startLabel start --endLabel end -o b150_g213_partial_${histone[i]}.pdf --plotFileFormat pdf --yMin 0 --yMax ${max[i]}

computeMatrix scale-regions -p 20 -R $bed_path/b150_g213_full_up_gene.bed $bed_path/b150_g213_full_down_gene.bed -S $bw_path/norm_b150_${histone[i]}_adjusted.wig.bw $bw_path/norm_g213_${histone[i]}_adjusted.wig.bw -a 3000 -b 3000 --regionBodyLength 9000 --scale 2 -o b150_g213_full_${histone[i]}.gz --outFileNameMatrix b150_g213_full_${histone[i]}.tab --outFileSortedRegions b150_g213_full_${histone[i]}.bed
plotHeatmap -m b150_g213_full_${histone[i]}.gz --colorMap coolwarm --startLabel start --endLabel end -o b150_g213_full_${histone[i]}.pdf --plotFileFormat pdf --yMin 0 --yMax ${max[i]}
done
deeeptools_v1

#5) 生成joint hyper需要的bed
#从hyper anno, hyper bed, DEGs中生成joint hyper需要的bed:
#cd E:\IDH_Hyper\Results\4_Gene_Expression\202208_oncogenes_Track\pairs_hist_diff\get_umr_hyper\extract_hyper_umr_bed
#perl extract_hyper_umr_bed.pl ../umr_200_mean0.1/brain149_AK076_umr_hyper.bed ../region_anno/brain149_AK076_partial_hyper_gene_anno.txt ../gene_bed/b149_g076_partial_up_gene.bed brain149_AK076_partial_up partial
#perl extract_hyper_umr_bed.pl ../umr_200_mean0.1/brain149_AK076_umr_hyper.bed ../region_anno/brain149_AK076_full_hyper_gene_anno.txt ../gene_bed/b149_g076_full_up_gene.bed brain149_AK076_full_up full

cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/get_umr_hyper/extract_hyper_umr_bed
type=(partial full)
for((i=0;i<${#type[@]};i++))
do
hyper_bed=../umr_200_mean0.1/brain149_AK076_umr_hyper.bed
hyper_anno=../region_anno/brain149_AK076_${type[i]}_hyper_gene_anno.txt
up_deg=../gene_bed/b149_g076_${type[i]}_up_gene.bed
down_deg=../gene_bed/b149_g076_${type[i]}_down_gene.bed
perl extract_hyper_umr_bed.pl $hyper_bed $hyper_anno $up_deg b149_g076_${type[i]}_up ${type[i]}
perl extract_hyper_umr_bed.pl $hyper_bed $hyper_anno $down_deg b149_g076_${type[i]}_down ${type[i]}

hyper_bed=../umr_200_mean0.1/brain150_AK213_umr_hyper.bed
hyper_anno=../region_anno/brain150_AK213_${type[i]}_hyper_gene_anno.txt
up_deg=../gene_bed/b150_g213_${type[i]}_up_gene.bed
down_deg=../gene_bed/b150_g213_${type[i]}_down_gene.bed
perl extract_hyper_umr_bed.pl $hyper_bed $hyper_anno $up_deg b150_g213_${type[i]}_up ${type[i]}
perl extract_hyper_umr_bed.pl $hyper_bed $hyper_anno $down_deg b150_g213_${type[i]}_down ${type[i]}
done

#partially - joint
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/get_umr_hyper/deeptools_joint
bed_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/get_umr_hyper/extract_hyper_umr_bed
bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
max=(30 10 10 10 40 10)
for((i=0;i<${#histone[@]};i++))
do
yMax=${max[i]}
#b149_g076
outname1=b149_g076_partial_hyper_${histone[i]}
outname2=b149_g076_partial_umr_${histone[i]}
outname=b149_g076_partial_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b149_g076_partial_up_hyper.bed $bed_path/b149_g076_partial_down_hyper.bed\
 -S $bw_path/norm_b149_${histone[i]}_adjusted.wig.bw $bw_path/norm_g076_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 2000 \
 --scale 2 \
 -o $outname1.gz --outFileNameMatrix $outname1.tab --outFileSortedRegions $outname1.bed

computeMatrix scale-regions -p 20 \
 -R $bed_path/b149_g076_partial_up_umr.bed $bed_path/b149_g076_partial_down_umr.bed\
 -S $bw_path/norm_b149_${histone[i]}_adjusted.wig.bw $bw_path/norm_g076_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 2000 \
 --scale 2 \
 -o $outname2.gz --outFileNameMatrix $outname2.tab --outFileSortedRegions $outname2.bed
perl new_joint_hyper_umr.pl $outname1.gz $outname2.gz $outname.gz
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
#b150_g213
outname1=b150_g213_partial_hyper_${histone[i]}
outname2=b150_g213_partial_umr_${histone[i]}
outname=b150_g213_partial_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b150_g213_partial_up_hyper.bed $bed_path/b150_g213_partial_down_hyper.bed\
 -S $bw_path/norm_b150_${histone[i]}_adjusted.wig.bw $bw_path/norm_g213_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 2000 \
 --scale 2 \
 -o $outname1.gz --outFileNameMatrix $outname1.tab --outFileSortedRegions $outname1.bed

computeMatrix scale-regions -p 20 \
 -R $bed_path/b150_g213_partial_up_umr.bed $bed_path/b150_g213_partial_down_umr.bed\
 -S $bw_path/norm_b150_${histone[i]}_adjusted.wig.bw $bw_path/norm_g213_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 2000 \
 --scale 2 \
 -o $outname2.gz --outFileNameMatrix $outname2.tab --outFileSortedRegions $outname2.bed
perl new_joint_hyper_umr.pl $outname1.gz $outname2.gz $outname.gz
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
done

#modify K4me3 - 150
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/get_umr_hyper/deeptools_joint
bed_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/get_umr_hyper/extract_hyper_umr_bed
bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
i=4
outname=b150_g213_partial_${histone[i]}
echo $outname
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax 50 \
 --colorMap coolwarm --plotFileFormat pdf --perGroup

#fully
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/get_umr_hyper/deeptools_fully
bed_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/get_umr_hyper/extract_hyper_umr_bed
bw_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/pairs_hist_diff/Data/adjust_ChIP
histone=(H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9me3)
max=(15 10 10 10 25 10)
for((i=0;i<${#histone[@]};i++))
do
#partial - up/down
yMax=${max[i]}
outname=b149_g076_fully_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b149_g076_full_up_hyper.bed $bed_path/b149_g076_full_down_hyper.bed \
 -S $bw_path/norm_b149_${histone[i]}_adjusted.wig.bw $bw_path/norm_g076_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 4000 \
 --scale 2 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup

outname=b150_g213_fully_${histone[i]}
computeMatrix scale-regions -p 20 \
 -R $bed_path/b150_g213_full_up_hyper.bed $bed_path/b150_g213_full_down_hyper.bed \
 -S $bw_path/norm_b150_${histone[i]}_adjusted.wig.bw $bw_path/norm_g213_${histone[i]}_adjusted.wig.bw \
 -a 3000 -b 3000 --regionBodyLength 4000 \
 --scale 2 \
 -o $outname.gz --outFileNameMatrix $outname.tab --outFileSortedRegions $outname.bed
plotHeatmap -m $outname.gz -o $outname.pdf \
 --startLabel start --endLabel end --yMin 0 --yMax $yMax \
 --colorMap coolwarm --plotFileFormat pdf --perGroup
done

