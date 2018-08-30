COMMIT_filtering () {

# subject
subj=$1

# model type
model_type=$2

# type of input tracks 
tracks_type=$3

# input and output folders
current_folder=$(pwd)
input_folder=/home/auguser2016/Projects/0001_Chiasm/$subj
output_folder=/home/auguser2016/Projects/Leicester_vs_LiFE

data_folder=/home/auguser2016/dMRI_DATA/PREPROCESSED_DATA

# initialization and preparation of input

# DWI
gunzip $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC.nii.gz -k #$input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC.nii
mrtransform $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC.nii -template $data_folder/Anatomies_ACPC_Aligned_Anonymized_mri_deface/${subj}_t1_acpc_anonim.nii.gz $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC_resized.nii -force
dwi=$input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC_resized.nii
#mrinfo $dwi

#dwi_scheme
fsl2scheme -bvecfile $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC.bvecs -bvalfile $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC.bvals > $input_folder/${subj}\_camino_grad_scheme.scheme
dwi_scheme=$input_folder/${subj}\_camino_grad_scheme.scheme

# peaks
# create peaks image and adjust it's dimensions so it corresponds to anatomy
#run calculations on new resized dwi to extract matching peaks
dwi2response tournier $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC_resized.nii -fslgrad $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvecs $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvals -shell 1600 -lmax 6 -mask $(dwi2mask $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC_resized.nii -fslgrad $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvecs $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvals - ) $input_folder/${subj}_COMMIT_sfr_lmax_6.txt -force
dwi2fod msmt_csd $(dwiextract $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC_resized.nii -fslgrad $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvecs $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvals - ) $input_folder/${subj}_COMMIT_sfr_lmax_6.txt -lmax 8 $input_folder/${subj}\_COMMIT_acpc_FOD_lmax_8.mif -mask $(dwi2mask $input_folder/${subj}\_150mm_clean_aligned_trilin_noMEC_resized.nii -fslgrad $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvecs $input_folder/${subj}_150mm_clean_aligned_trilin_noMEC.bvals - ) -force
sh2peaks $input_folder/${subj}\_COMMIT_acpc_FOD_lmax_8.mif $input_folder/${subj}_camino_sh2peaks.nii.gz -force
peaks=$input_folder/${subj}_camino_sh2peaks.nii.gz
#mrinfo $peaks

# tracks
# convert tracks from mrtrix do camino's .trk
python /home/auguser2016/Software/tck2trk.py $data_folder/Anatomies_ACPC_Aligned_Anonymized_mri_deface/${subj}_t1_acpc_anonim.nii.gz $output_folder/Leicester_$tracks_type/${subj}/${subj}\_ACT.tck -f
tracks=$output_folder/Leicester_$tracks_type/${subj}/${subj}\_ACT.trk

# wm
# create non-cropped 5tt image, transform correcting patch to it's dimensions, adjust wm definition with manually prepared patch and extract wm from given 5tt
5ttgen fsl $data_folder/Anatomies_ACPC_Aligned_Anonymized_mri_deface/${subj}_t1_acpc_anonim.nii.gz $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_t1_acpc_5tt_nocrop.nii.gz -nocrop -force
mrtransform $data_folder/5TTs_ACPC_Aligned_Chiasm_T1_Corrected/${subj}_patch.mif -template $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_t1_acpc_5tt_nocrop.nii.gz $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_patch_fit_to_nocrop.nii.gz -force
5ttedit $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_t1_acpc_5tt_nocrop.nii.gz -wm $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_patch_fit_to_nocrop.nii.gz $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_t1_acpc_5tt_nocrop_corrected.nii.gz -force
mrconvert -coord 3 2 -axes 0,1,2 $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_t1_acpc_5tt_nocrop_corrected.nii.gz $data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_wm_from_t1_corrected_5tt2.nii.gz -force
wm=$data_folder/WM_Mask_from_T1_Corrected_5TT/${subj}_wm_from_t1_corrected_5tt2.nii.gz
#mrinfo $wm

# COMMIT filtering run in python

python /home/auguser2016/Projects/Leicester_vs_LiFE/COMMIT/1_apply_commit.py $dwi $dwi_scheme $peaks $tracks $wm $output_folder/COMMIT/$model_type/$tracks_type/$subj -model $model_type
python /home/auguser2016/Projects/Leicester_vs_LiFE/COMMIT/2_present_results.py $subj $tracks_type $model_type $output_folder/COMMIT/$model_type/$tracks_type/$subj
python /home/auguser2016/Software/trk2tck.py /home/auguser2016/Projects/Leicester_vs_LiFE/COMMIT/$model_type/$tracks_type/${subj}/${subj}_${tracks_type}_${model_type}_streamlines.trk
# extraction of fibers (output organized exactly as within LIFE workframe)

# run script based on LiFE-based filtering script (no-filtering variant)
    
    }

models=" StickZeppelinBall_Model LiFE_Model"
tracking_type='DT_seed DT_selected CSD_seed CSD_selected'

for a in $models; do

	mkdir /home/auguser2016/Projects/Leicester_vs_LiFE/COMMIT/$a

	for b in $tracking_type; do
		mkdir /home/auguser2016/Projects/Leicester_vs_LiFE/COMMIT/$a/$b
	

		for probanden in /home/auguser2016/Projects/0001_Chiasm/*; do
			COMMIT_filtering ${probanden: -4} $a $b 
		done

	done
	
done





