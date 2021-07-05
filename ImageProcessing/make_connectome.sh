#!/bin/bash

# cd into data directory
cd /rawdata


for i in sub-*; do
    echo $i

    # set location of structural and diffusion data
    atlas_base=../sourcedata/derivatives/nipype/${i}/ses-preop/anatomical_pipeline/parcellation_stage/
    diffusion_base=$i/ses-preop/dwi/mrtrix
    
    # register structural, atlas and diffusion data
    reg_aladin -flo $atlas_base/Lausanne2018_parcellation/brain.nii.gz -ref $diffusion_base/nodif.nii.gz -aff t12diff.txt -res $diffusion_base/diff_space_brain.nii.gz -rigOnly

    reg_resample -flo $atlas_base/parcCombiner/ROIv_Lausanne2018_scale3_final.nii.gz -ref $diffusion_base/nodif.nii.gz -aff t12diff.txt -res $diffusion_base/diff_space_labels.nii.gz -inter 0
    
    reg_resample -flo $atlas_base/Lausanne2018_parcellation/aparc+aseg.native.nii.gz -ref $diffusion_base/nodif.nii.gz -aff t12diff.txt -res $diffusion_base/diff_space_aparc+aseg.nii.gz -inter 0

    # run 5ttgen
    
    5ttgen freesurfer $diffusion_base/diff_space_aparc+aseg.nii.gz $diffusion_base/5tt.mif

    # generate cortical targets
    
    labelconvert $diffusion_base/diff_space_labels.nii.gz Scale3_OldLabels.txt Scale3_NewLabels.txt $diffusion_base/targets.nii.gz -force

    #create structural connectome
    
    tckgen $diffusion_base/wm.mif -act $diffusion_base/5tt.mif -select 5000000 -seed_dynamic $diffusion_base/wm.mif $diffusion_base/5M.tck

    tcksift2 $diffusion_base/5M.tck $diffusion_base/wm.mif $diffusion_base/5M_sift.txt

    tck2connectome -symmetric -tck_weights_in $diffusion_base/5M_sift.txt $diffusion_base/5M.tck $diffusion_base/targets.nii.gz $diffusion_base/newconnectome.csv -force
    
done
