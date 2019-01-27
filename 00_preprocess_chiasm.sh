# Initializing

# Folder with raw and preprocessed data
data=/home/auguser2016/dMRI_DATA
#scripts=/home/auguser2016/Scripts
mkdir $data/RDY_DATA/1.5mm_iso

General_preprocessing () {

	# Preprocessing image in original resolution (1.5 mm iso), inside each subjects' folder all preprocessed data sets with masks, map of distortions etc.
	
	# Get code for given subject
	subj=${1::-6} 
	
	# create folders for preprocessing
	low_res=RDY_DATA/1.5mm_iso/$subj

	mkdir $data/$low_res
	mkdir $data/$low_res/Temporal
	mkdir $data/$low_res/Temporal/Distortions

	#mkdir $data/$high_res
	#mkdir $data/$high_res/Temporal
	
	cd $data/RAW_DATA/${1::-1}/s*/

	# Extract raw data using mrtrix designated command - mrconvert
	mrconvert *ep2d*AP $data/$low_res/Temporal/$subj\_raw_AP.mif
	mrconvert *ep2d*PA $data/$low_res/Temporal/$subj\_raw_PA.mif

	cd $data/$low_res/Temporal
		
	# Denoise data
	dwidenoise $data/$low_res/Temporal/$subj\_raw_AP.mif $data/$low_res/Temporal/$subj\_denoised_AP.mif -noise $data/$low_res/Temporal/Distortions/$subj\_noise_AP.mif
	dwidenoise $data/$low_res/Temporal/$subj\_raw_PA.mif $data/$low_res/Temporal/$subj\_denoised_PA.mif -noise $data/$low_res/Temporal/Distortions/$subj\_noise_PA.mif

	# Remove Gibbs ringing
	mrdegibbs $data/$low_res/Temporal/$subj\_denoised_AP.mif $data/$low_res/Temporal/$subj\_denoised_deringed_AP.mif
	mrdegibbs $data/$low_res/Temporal/$subj\_denoised_PA.mif $data/$low_res/Temporal/$subj\_denoised_deringed_PA.mif

	# Combine data with both encoding direction into one (as requested by dwipreproc)
	mrcat $data/$low_res/Temporal/$subj\_denoised_deringed_AP.mif $data/$low_res/Temporal/$subj\_denoised_deringed_PA.mif $data/$low_res/Temporal/$subj\_AP_PA_preeddy.mif

	# Run dwipreproc on prepared data
	dwipreproc $data/$low_res/Temporal/$subj\_AP_PA_preeddy.mif $data/$low_res/Temporal/$subj\_postproc.mif -rpe_all -pe_dir AP -cuda -debug 

	# Obtain mask
	dwi2mask $data/$low_res/Temporal/$subj\_postproc.mif $data/$low_res/Temporal/$subj\_150mm_mask.mif

	# Masked is used for bias field correction using ANTS software
	dwibiascorrect -ants $data/$low_res/Temporal/$subj\_postproc.mif $data/$low_res/$subj\_150mm_clean.mif -mask $data/$low_res/Temporal/$subj\_150mm_mask.mif -bias $data/$low_res/Temporal/Distortions/$a\_bias_field.mif
	
}

cd $data/RAW_DATA

# List of all data_sets used in studies
#list=(ce04_2054/ dn20_0713/ fe21_1325/ hw91_0844/ ib57_0731/ kw99_0633/ la21_1353/ lw37_0977/ nb30_1185/ ow93_0974/ ps94_1516/ rx88_1234/ sj22_1218/ ta14_1065/ tq63_1214/ uf97_1072/ uh47_1309/ xn78_1085/ xs62_1217/ )

list=(nh50_2693/)

for i in $list; do General_preprocessing $i;done	

for i in ${list[@]:0:2}; do General_preprocessing $i & done
wait
for i in ${list[@]:2:4}; do General_preprocessing $i & done
wait
for i in ${list[@]:4:6}; do General_preprocessing $i & done
wait
for i in ${list[@]:6:8}; do General_preprocessing $i & done
wait
for i in ${list[@]:8:10}; do General_preprocessing $i & done
wait
for i in ${list[@]:10:12}; do General_preprocessing $i & done
wait
for i in ${list[@]:12:14}; do General_preprocessing $i & done
wait
for i in ${list[@]:14:16}; do General_preprocessing $i & done
wait
for i in ${list[@]:16:18}; do General_preprocessing $i & done
wait
for i in ${list[@]:18:20}; do General_preprocessing $i & done
wait
for i in ${list[@]:20:22}; do General_preprocessing $i & done
wait


cd $data/RDY_DATA


#cd $data/RAW_DATA




# This way 4 sets are processed at once - we can't process all at once due to RAM usage, however preprocessing one at a time is ineffective because e.g. eddy from dwipreproc is using only one core. There is a version of eddy operating on multicore processors, but it's not running on this machine



# Returning to original folder
#cd /media/auguser2016/Volume/Test/chiasm

END





: <<'END'



	# Preprocessing with FSL, direct approach (!)
	
	for b in AP PA;
	do
		#convert to FSL format
		#mrconvert $data/$low_res/Temporal/$subj\_denoised_deringed_$b.mif $data/$low_res/Temporal/$subj\_denoised_deringed_$b.nii -export_grad_fsl $data/$low_res/Temporal/bvecs_$b $data/$low_res/Temporal/bvals_$b

		# create temporary folder for splitting
		#mkdir $data/$low_res/Temporal/$b
		#fslsplit $data/$low_res/Temporal/$subj\_denoised_deringed_$b.nii $data/$low_res/Temporal/$b/$b\_
		
		#cd $data/$low_res/Temporal/$b
		#cp *0000.nii.gz nodif00.nii.gz
		#cp *0001.nii.gz nodif01.nii.gz
		#cp *0018.nii.gz nodif02.nii.gz
		#cp *0035.nii.gz nodif03.nii.gz
		#cp *0052.nii.gz nodif04.nii.gz
		#cp *0069.nii.gz nodif05.nii.gz
		#cp *0086.nii.gz nodif06.nii.gz
		#cp *0103.nii.gz nodif07.nii.gz
		#cp *0120.nii.gz nodif08.nii.gz
		#cp *0137.nii.gz nodif09.nii.gz
		#fslmerge -t b0_$b nodif*
		#rm $b\_*
		#rm nodif*
	echo ' ' 
	done

	# Combine all extracted b0 for topup (we already have DWI for AP and PA, bvecs and bvals for AP and PA)
	cd $data/$low_res/Temporal 
	#fslmerge -t $subj\_all_b0 $data/$low_res/Temporal/AP/b0_AP $data/$low_res/Temporal/PA/b0_PA

	# Cleaning unneeded files

	# Preparing text files for topup (index file describing aquisition parameters - line per image)

	# Prepare merged bvecs and bvals files
	touch bvals_up_t;
	touch bvecs_up_t;
	touch bvals_dn_t;
	touch bvecs_dn_t;
	touch bvals_double_t;
	touch bvecs_double_t;

	awk -f $scripts/Additional_Required_Files/transpose.txt bvals_AP > bvals_up_t;
	awk -f $scripts/Additional_Required_Files/transpose.txt bvecs_AP > bvecs_up_t;
	awk -f $scripts/Additional_Required_Files/transpose.txt bvals_PA > bvals_dn_t;
	awk -f $scripts/Additional_Required_Files/transpose.txt bvecs_PA > bvecs_dn_t;

	cat bvals_up_t bvals_dn_t > bvals_double_t;
	cat bvecs_up_t bvecs_dn_t > bvecs_double_t;

	awk -f $scripts/Additional_Required_Files/transpose.txt bvals_double_t > bvals_ready;
	awk -f $scripts/Additional_Required_Files/transpose.txt bvecs_double_t > bvecs_ready;

	rm bvals_up_t;
	rm bvecs_up_t;
	rm bvals_dn_t;
	rm bvecs_dn_t;
	rm bvals_double_t;
	rm bvecs_double_t;
	
	#mkdir $data/$low_res/Temporal/Topup

	#Topup
	#topup --imain=$data/$low_res/Temporal/$subj\_all_b0.nii.gz --datain=$scripts/Additional_Required_Files/blip_up_dn.txt --config=$scripts/Additional_Required_Files/b02b0.cnf --out=Topup/topup_results --iout=Topup/hifi_b0 --fout=Topup/field --logout=Topup/topup.log -v
	# out - output spline coefficients and movement parameters (2 files actually)
	# iout - 4D imwage file with unwarped images
	# fout - image file with field

	# Obtain brain mask from Topup output
	#mkdir $data/$low_res/Temporal/Brain_Mask
	#fslmaths $data/$low_res/Temporal/Topup/hifi_b0 -Tmean $data/$low_res/Temporal/Brain_Mask/hifi_b0_mean # averaging all b=0 images
	#bet $data/$low_res/Temporal/Brain_Mask/hifi_b0_mean $data/$low_res/Temporal/Brain_Mask/hifi_b0_mean_brain -m -f 0.3 # extracting brain mask

	# Merge AP and PA data
	#fslmerge -t $data/$low_res/Temporal/$subj\_AP_PA_preeddy $data/$low_res/Temporal/$subj\_denoised_deringed_AP.nii $data/$low_res/Temporal/$subj\_denoised_deringed_PA.nii

	# Eddy
	#mkdir $data/$low_res/Temporal/Eddy
	#eddy_openmp --imain=$data/$low_res/Temporal/$subj\_AP_PA_preeddy.nii --mask=$data/$low_res/Temporal/Brain_Mask/hifi_b0_mean_brain_mask.nii.gz --acqp=$scripts/Additional_Required_Files/blip_up_dn.txt --index=$scripts/Additional_Required_Files/index.txt --bvecs=$data/$low_res/Temporal/bvecs_ready --bvals=$data/$low_res/Temporal/bvals_ready --topup=$data/$low_res/Temporal/Topup/topup_results --out=$data/$low_res/Temporal/Eddy/data_ed --verbose
 
	# Correct bias field 
	#dwibiascorrect -fsl $data/$low_res/Temporal/Eddy/data_ed.nii.gz $data/$low_res/$subj\_clean_150mm.nii.gz -mask $data/$low_res/Temporal/Brain_Mask/hifi_b0_mean_brain_mask.nii.gz -bias $data/$low_res/Temporal/Distortions/bias_field.nii.gz -fslgrad $data/$low_res/Temporal/Eddy/data_ed.eddy_rotated_bvecs $data/$low_res/Temporal/bvals_AP -debug

	# Alternative comparison
	

# Re-extraction of brain mask from eddy-corrected data
#bet2 Eddy/data_eddy.nii.gz Brain_Mask/eddy_corrected_data_brain -m -f 0.3

# Upsampling DWI image for sake of registration and future tracking
#mrresize Eddy/data_eddy.nii.gz -scale 2 Eddy/data_corrected_upsampled_075.nii.gz -force

# Executing next step
#source $scripts/4_registration.txt $1

	#Eddy

		


		

	
	# Obtaining mask from preprocessed data

	# Bias field correction

	# 



   Combine DWI series into one - first half AP, second PA along time axis
	mrcat $data/$low_res/Temporal/$subj_denoised_deringed_AP.mif $data/$low_res/Temporal/$subj_denoised_deringed_PA.mif $data/$low_res/$subj_denoise_dering_combined.mif -axis 3

	# Preprocess data (motion-correction, topup, eddy)




	
	#mrcat $data/$low_res/Temporal/*AP_denoised.mif $data/$low_res/Temporal/*PA_denoised.mif $data/$low_res/Temporal/$a\_combined_DWIs.mif -axis 3 -force

	# Script from MRtrix interacting with fsl. Here version for correcting DWI with 2 phase encoding directions is used
	dwipreproc $data/$low_res/Temporal/$a\_combined_DWIs.mif $data/$low_res/Temporal/$a\_denoised_preproc.mif -pe_dir AP -rpe_all -export_grad_mrtrix $data/$low_res/Temporal/Distortions/$a\_gradients 

 #-pe_dir AP -rpe_all $data/$low_res/Temporal/*PA_denoised.mif  $data/$low_res/Temporal/*AP_denoised.mif  $data/$low_res/Temporal/$a\_denoised_preproc.mif -export_grad_mrtrix $data/$low_res/Temporal/Distortions/$a\_gradients 
	
	# Masking preprocessed data. Results should be checked
	dwi2mask $data/$low_res/Temporal/$a\_denoised_preproc.mif $data/$low_res/$a\_clean_150mm_mask.mif

	# Masked is used for bias field correction using ANTS software
	dwibiascorrect -ants $data/$low_res/Temporal/$a\_denoised_preproc.mif $data/$low_res/$a\_clean_150mm.mif -mask $data/$low_res/$a\_clean_150mm_mask.mif -bias $data/$low_res/Temporal/Distortions/$a\_bias_field.mif
 
	# Preprocessing of 1.5 mm iso data sets done. Moving to upsampled resolution
	
	# Upsampling by factor of 2 to 0.75mm. Justified step, in FBA even factor 3 can be recommended
	mrresize $data/$low_res/$a\_clean_150mm.mif -scale 2.0 $data/$high_res/$a\_clean_075mm.mif
	
	# Calculating mask again
	dwi2mask $data/$high_res/$a\_clean_075mm.mif $data/$high_res/$a\_clean_075mm_mask.mif

	# Estimating single fibre response (SFR) for default lmax
	dwi2response tournier $data/$high_res/$a\_clean_075mm.mif $data/$high_res/$a\_clean_075mm_SFR.txt -shell 1600 -lmax 8 -mask $data/$high_res/$a\_clean_075mm_mask.mif -voxels $data/$high_res/Temporal/$a\_clean_075mm_SFR_voxels.mif

	# Estimating FOD for previously obtained SFR
	dwiextract $data/$high_res/$a\_clean_075mm.mif - | dwi2fod msmt_csd - $data/$high_res/$a\_clean_075mm_SFR.txt $data/$high_res/$a\_clean_075mm_FOD.mif -mask $data/$high_res/$a\_clean_075mm_mask.mif

	# Performing coregistration of anatomical image to DWI images and segmenting it. Its easier approach as gradient table is not rotated this way. However for cross-modal studies it would be better to register DWI to anatomy

	cd *MPRAGE*iso
	
	# Converting DWI image from .mif to .nii.gz as FSL likes. This step should be checked as well.
	mrconvert $data/$high_res/$a\_clean_075mm.mif $data/$high_res/Temporal/$a\_clean_075mm.nii.gz
	
	# Coping anatomical image to corresponding folder
	cp o* $data/$high_res/$a\_anatomical.nii.gz

	# Extracting clean brain from anatomical image. Check this step
	bet2 $data/$high_res/$a\_anatomical.nii.gz $data/$high_res/Temporal/$a\_anatomical_brain.nii.gz -f 0.6  	
	
	# Corregister DWI to anatomy
	epi_reg --epi=$data/$high_res/Temporal/$a\_clean_075mm.nii.gz --t1=$data/$high_res/$a\_anatomical.nii.gz --t1brain=$data/$high_res/Temporal/$a\_anatomical_brain.nii.gz --out=$data/$high_res/Temporal/$a\_nodif2mprage_2 --echospacing=0.00023 --pedir=-y

	# Reverse corregistration matrix
	convert_xfm -inverse $data/$high_res/Temporal/$a\_nodif2mprage_2.mat -omat $data/$high_res/Temporal/$a\_mprage2nodif_2.mat

	# Apply reverted matrix to corregister anatomy to DWI
	flirt -in $data/$high_res/$a\_anatomical.nii.gz -ref $data/$high_res/Temporal/$a\_clean_075mm.nii.gz -out $data/$high_res/$a\_075mm_anatomical_coreg2nodif_undist_2 -applyxfm -init $data/$high_res/Temporal/$a\_mprage2nodif_2.mat -interp sinc -sincwidth 7 -sincwindow hanning

	# Perform segmentation into WM, GM, subcortical GM, CSF and possible pathological tissue using FIRST from FSL
	5ttgen fsl $data/$high_res/$a\_075mm_anatomical_coreg2nodif_undist_2.nii.gz $data/$high_res/$a\_075mm_5tt.nii.gz
'


