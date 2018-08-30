% Add required folders to path
addpath(genpath('/home/auguser2016/Software/spm12'))
addpath(genpath('/home/auguser2016/Software/vistasoft'))
addpath(genpath('/home/auguser2016/Software/mba'))
addpath(genpath('/home/auguser2016/Software/encode'))

cd /home/auguser2016/dMRI_DATA/RDY_DATA/1.5mm_iso

subjects=dir

for i=11:size(subjects,1)
	id = subjects(i).name

	% Load data

	anatomy=strcat('/home/auguser2016/dMRI_DATA/PREPROCESSED_DATA/Anatomies_ACPC_Aligned_MIPAV/',id,'_t1_acpc_anonym.nii.gz')
	dwi=    strcat('/home/auguser2016/dMRI_DATA/RDY_DATA/1.5mm_iso/',id,'/',id,'_150mm_clean.nii.gz')
	bvecs=  strcat('/home/auguser2016/dMRI_DATA/RDY_DATA/1.5mm_iso/',id,'/',id,'_bvecs')
	bvals=  strcat('/home/auguser2016/dMRI_DATA/RDY_DATA/1.5mm_iso/',id,'/',id,'_bvals')
      
    % Create output folder

    mkdir(strcat('/home/auguser2016/Projects/0001_Chiasm/',id))
	output=strcat('/home/auguser2016/Projects/0001_Chiasm/',id)
  
	% No flip
  
    	% Set up dtiInit paramters specific for our study
    	dwParams = dtiInitParams;
    	dwParams.rotateBvecsWithRx = 1;
    	dwParams.rotateBvecsWithCanXform = 1;
    	dwParams.bvecsFile = fullfile(bvecs);
    	dwParams.bvalsFile = fullfile(bvals);
    	dwParams.eddyCorrect = -1;
    	dwParams.dwOutMm = [1.5 1.5 1.5]; % change this to 1.5 and process 0.75 anyway
    	%dwParams.dwOutMm = [1.5 1.5 1.5];
    	dwParams.phaseEncodeDir = 2;
    	dwParams.outDir=output
    	%dwParams.outDir = fullfile(pwd,'dtiInitPreprocessed');

    	%Start Diffusion Imaging preprocessing.
    	[dt6FileName, outBaseDir] = dtiInit(dwi,...
                                        anatomy, ...
                                        dwParams); 
    
    	%mrtrix_bfileFromBvecs(strcat('/home/auguser2016/Projects/Dir_test/',id,'_0/nb30_clean_150mm_aligned_trilin_noMEC.bvecs'),strcat('/home/auguser2016/Projects/Dir_test/',id,'_0/nb30_clean_150mm_aligned_trilin_noMEC.bvals'),strcat('/home/auguser2016/Projects/Dir_test/',id,'_0/',id,'_gtab.b'))

end

