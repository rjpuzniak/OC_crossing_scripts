 #!/bin/bash

General_postprocessing_mrtrix () {

subj=${1}

analysis=/home/auguser2016/Projects/0001_Chiasm/$subj

mrconvert $analysis/${subj}\_150mm_clean_aligned_trilin_noMEC.nii.gz -fslgrad $analysis/${subj}\_150mm_clean_aligned_trilin_noMEC.bvecs $analysis/${subj}\_150mm_clean_aligned_trilin_noMEC.bvals $analysis/${subj}\_dwi_ACPC.mif

#mrconvert  $analysis/${subj}\_150mm_clean_aligned_trilin_noMEC.nii.gz -grad $analysis/${subj::-2}_gtab.b $analysis/${subj}\_dwi_ACPC_mrdif.mif

dwi2mask $analysis/${subj}\_dwi_ACPC.mif $analysis/${subj}\_acpc_mask.nii.gz -force

# Estimating single fibre response (SFR) for lmax=6
dwi2response tournier $analysis/${subj}\_dwi_ACPC.mif $analysis/${subj}\_acpc_SFR.txt -shell 1600 -lmax 6 -mask $analysis/${subj}\_acpc_mask.nii.gz  -voxels $analysis/${subj}\_acpc_SFR_voxels.nii.gz 

# Estimation of Fiber Orientation Distribution function for 3 lmax
for lmax in 6 8 10; do  
  dwiextract $analysis/${subj}\_dwi_ACPC.mif - | dwi2fod msmt_csd - $analysis/${subj}\_acpc_SFR.txt -lmax $lmax $analysis/${subj}\_acpc_FOD_lmax_$lmax.mif -mask $analysis/${subj}\_acpc_mask.nii.gz 
done
}

General_postprocessing_mrdif () {

subj=${1}

analysis=/home/auguser2016/Projects/Dir_test/$subj

mrconvert  $analysis/${subj::-2}\_clean_150mm_aligned_trilin_noMEC.nii.gz -grad $analysis/${subj::-2}_gtab.b $analysis/${subj}\_dwi_ACPC_mrdif.mif

dwi2mask $analysis/${subj}\_dwi_ACPC_mrdif.mif $analysis/${subj::-2}\_acpc_mask_mrdif.nii.gz -force

# Estimating single fibre response (SFR) for lmax=6
dwi2response tournier $analysis/${subj}\_dwi_ACPC_mrdif.mif $analysis/${subj}\_acpc_SFR_mrdif.txt -shell 1600 -lmax 6 -mask $analysis/${subj::-2}\_acpc_mask_mrdif.nii.gz  -voxels $analysis/${subj}\_acpc_SFR_voxels_mrdif.nii.gz 

# Estimation of Fiber Orientation Distribution function for 3 lmax
for lmax in 6 8 10 12; do  
  dwiextract $analysis/${subj}\_dwi_ACPC_mrdif.mif - | dwi2fod msmt_csd - $analysis/${subj}\_acpc_SFR_mrdif.txt -lmax $lmax $analysis/${subj}\_acpc_FOD_lmax_$lmax\_mrdif.mif -mask $analysis/${subj::-2}\_acpc_mask_mrdif.nii.gz 
done
}


current=$(pwd)

cd /home/auguser2016/Projects/0001_Chiasm

for i in */; do 
General_postprocessing_mrtrix ${i::4};
done

cd $current
