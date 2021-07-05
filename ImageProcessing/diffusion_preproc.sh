#!/bin/bash

# arguments not required!

# dwi       - multishell diffusion dataset (folder containing DICOM)
# dwi_negpe - negative phase-encoded b0 image (folder containing DICOM)

# specify file path here
for i in /rawdata/sub-*/ses-preop/dwi; do 
    if [ -d $i ]; then
	echo $i
	cd $i
    
	# extract b0 image and negpe b0 image
	mrconvert negPE dwi_negpe.nii.gz -stride 1,2,3
	mrconvert dwi dwi.nii.gz -stride 1,2,3,4
   
	mv dwi_negpe.nii.gz b0_flip.nii.gz
	fslroi dwi.nii.gz b0.nii.gz 0 1
   
	fslmerge -t b0s_all.nii.gz b0.nii.gz b0_flip.nii.gz
	mkdir mrtrix
    
	mrconvert b0s_all.nii.gz mrtrix/dwi.mif

	#denoising

	dwidenoise mrtrix/dwi.mif mrtrix/out.mif -noise noise.mif
	
	#eddy correction
	
	dwipreproc -rpe_pair -se_epi mrtrix/out.mif -pe_dir AP -eddy_options "--repol " dwi mrtrix/dndwi.mif
	mrconvert mrtrix/dndwi.mif  mrtrix/data.nii.gz -stride 1,2,3,4
	fslroi  mrtrix/data.nii.gz  mrtrix/nodif.nii.gz 0 1
	bet  mrtrix/nodif.nii.gz  mrtrix/brain -m -n -f 0.3
	mrconvert mrtrix/brain_mask.nii.gz mrtrix/mask.mif

	# add bioas correction step
    
	dwibiascorrect -ants -mask mrtrix/mask.mif mrtrix/dndwi.mif mrtrix/dnbcdwi.mif
	dwi2tensor mrtrix/dnbcdwi.mif -mask mrtrix/mask.mif mrtrix/dt.mif

	# get tensors
    
	tensor2metric mrtrix/dt.mif -fa mrtrix/fa.mif -adc mrtrix/md.mif -vector mrtrix/ev.mif
   
	# ODF response using CSD

	dwi2response dhollander mrtrix/dnbcdwi.mif mrtrix/wm_response.txt mrtrix/gm_response.txt mrtrix/csf_response.txt -nocleanup
	dwi2fod msmt_csd -mask mrtrix/mask.mif mrtrix/dnbcdwi.mif mrtrix/wm_response.txt mrtrix/wm.mif mrtrix/gm_response.txt gm.mif mrtrix/csf_response.txt mrtrix/csf.mif
    
	cd -
    fi
done
