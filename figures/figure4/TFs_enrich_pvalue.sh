#记录使用Homer对甲基化文件进行known motifs富集分析

#v1.0
#1) p_up; 背景: p_down
#参数说明:
#findMotifsGenome.pl 输入bed区域 基因组 输出文件夹
#-size given 给定输入区域
#-nomotif 不de novo发现motifs (速度很快, 几百个区域几分钟即可完成注释)
#-mknown 使用脊椎动物的已知motifs进行注释
#-bg 背景control区域
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/homer
/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl ../regions/sort_partial_up_hyper.bed hg19 ./p_up_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg ../regions/sort_partial_down_hyper.bed -p 10

#2) p_down; 背景: p_up
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/homer
/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl ../regions/sort_partial_down_hyper.bed hg19 ./p_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg ../regions/sort_partial_up_hyper.bed -p 10

#3) f_down; 背景: p_up
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/homer
/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl ../regions/sort_full_down_hyper_hyper.bed hg19 ./f_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg ../regions/sort_partial_up_hyper.bed -p 10

#v1.1
#在基础上设置为超几何分布
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/homer/homer1.1
region_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/regions
/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_partial_up_hyper.bed hg19 ./p_up_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_down_hyper.bed -p 10 -h

/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_partial_down_hyper.bed hg19 ./p_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_up_hyper.bed -p 10 -h

/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_full_down_hyper_hyper.bed hg19 ./f_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_up_hyper.bed -p 10 -h

#v1.2
#设置float
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/homer/homer1.2
region_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/regions
/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_partial_up_hyper.bed hg19 ./p_up_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_down_hyper.bed -p 10 -float 

/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_partial_down_hyper.bed hg19 ./p_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_up_hyper.bed -p 10 -float

/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_full_down_hyper_hyper.bed hg19 ./f_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_up_hyper.bed -p 10 -float

#v1.3
#关闭GC bias
#几乎没有差别
cd /share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/homer/homer1.3
region_path=/share/pub/wangxy/IDH_Glioma/4_Gene_Expression/202208_oncogenes_Track/motif_analysis/regions
/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_partial_up_hyper.bed hg19 ./p_up_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_down_hyper.bed -p 10 -noweight

/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_partial_down_hyper.bed hg19 ./p_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_up_hyper.bed -p 10 -noweight

/share/pub/wangxy/software/miniconda3/share/homer/bin/findMotifsGenome.pl $region_path/sort_full_down_hyper_hyper.bed hg19 ./f_down_homer -size given -nomotif -mknown /share/pub/wangxy/software/homer/data/knownTFs/vertebrates/known.motifs -bg $region_path/sort_partial_up_hyper.bed -p 10 -noweight
