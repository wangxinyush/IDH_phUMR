* Here is the computational framework to identify partially and fully hypermethylated UMRs.
* **average_methy:** Calculated average methylation level of WGBS samples for normal and cancer samples.
* Example of Usage: 
```
average_methy -n normal_brain_wig_list.txt -c IDH_glioma_wig_list.txt -o mean_methy_output
```
* **hyper_state:** identify three methylation states sites by using HMM.
* Example of Usage: 
```
RScript hyper_state.r mean_methy.wig HMM_state.txt
```
* **merge_hmm_DMC:** merge hmm and DMC(from statistical test) into a matrix.
* Example of Usage: 
```
merge_hmm_DMC HMM_states.txt DMC_state.txt refumr_CG_DMC_stat_HMM.txt
```
* **refumr_CG_methy:** get CpGs in refumr and extract methylation matrix with methylation status.
* Example of Usage:
```
refumr_CG_methy reference_UMR.bed wig_list.txt refumr_CG_methy_matrix.txt
```
* **CG_to_region:** merge single hyper CpG into a region.
* Example of Usage:
```
CG_to_region hyper_refumr_state_stat.txt refumr_CG_DMC_stat_HMM.txt hyper_region_in_refumr.txt
```
* The above 5 commands were used to identifiy phUMRs and fhUMRs for multiple noraml and cancer samples
* **phumr_twosamples.pl:** get phUMRs and fhUMRs for only two WGBS samples.
* Example of Usage:
```
perl phumr_twosamples.pl Roadmap150_umr.bed hg38_CpG.txt wig_list.txt >Roadmap150_umr_hyper.bed
```
