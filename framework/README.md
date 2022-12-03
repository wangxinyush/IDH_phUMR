* Here is the computational framework to identify partially and fully hypermethylated UMRs.
* 1) average_methy.py Calculated average methylation level of WGBS samples for normal and cancer samples.
* 2) hyper_state.r. identify three methylation states sites by using HMM.
* 3) merge_hmm_DMC.pl. merge hmm and DMC(from statistical test) into a matrix.
* 4) refumr_CG_methy.pl. get CpGs in refumr and extract methylation matrix with methylation status.
* 5) CG_to_region.pl. merge single hyper CpG into a region.
* phumr_twosamples.pl: get phUMRs and fhUMRs for only two WGBS samples.
