# OC_crossing_scripts

## Remarks


## Explanation of pipeline
How to repeat the analysis:

LiFE Analysis 
-------------
1) Preprocessed data is uploaded to BL and is to be accessed here: https://brainlife.io/project/5aaa347ef0b5260027e24a96/detail (You have permission to access it)
2) Run dtiInit on it - script 01_t1_and_dwi_to_acpc_chiasm.m
3) Estimating SFR and ODFs from output of dtiInit - script 02_postprocess_chiasm.sh
4) Script 03_Leicester_vs_LiFE_tracking.txt performs tracking - different approaches (CSD/DT, merged A->B and B->A, multiple parameters, 'seed' indicates fixed number of seeding attempts, 'select' indicates fixed number of streamlines we select for each group). Output will be saved to folder Leicester_XXX_YYY/SUBJ, where XXX is model used for tractography (CSD or DT), YYY is either SEED or SELECTED, SUBJ is subject's ID. Output consists of 5 files - 4 files SUBJ_? _?_ACT.tck where ? can be either r (r=right) or l (l=left). This indicate side of starting and ending ROI. File SUBJ_ACT.tck consists of 4 merged files SUBJ_?_?_ACT.tck. 
5) Script 04_Leicester_vs_LiFE_filtering.m performs LiFE filtering - script creates a pair of output folders - one with non-filtered tractogram broken done into 4 fiber bundles (SUBJ_noLIFE), one with LiFE filtered tractogram broken down into 4 bundles (SUBJ_postLiFE). In order to select those sub-bundles a script 07_connecting_rois.sh is invoked. Files organized as in step 4) - one pair per one combination of CSD/DT and SEED/SELECTED. Altogether 2x2x2 = 8 folders, each containing 5 .tck files. Additionally, MATLAB data files with tables are created in main folder. Those files contain number of streamlines for every fiber group and every participant for given combination of CSD/DT and SEED/SELECTED. 
6) Run statistics on results and plot figures. Returns dot plots for every combination and summary plot with decussation indices - 05_Plot_all_together.m. This script takes also result from COMMIT analysis, so this one must be done before this step.
--------------



COMMIT Analysis
---------------
1) as for LiFE analysis, they're using same data and same scripts.
2) as for LiFE analysis, they're using same data and same scripts.
3) as for LiFE analysis, they're using same data and same scripts.
4) as for LiFE analysis, they're using same data and same scripts.
5) Prepare data to be run through COMMIT framework - 0_optic_chiasm_commit.sh. Upon completion this scripts automatically runs 1_apply_commit.py, which performs COMMIT filtering for Stick-Zeppelin-Ball models and LiFE (single stick) model, as well as script 2_present_results.py which selects streamlines with non-zero weights and saves them in a format to be fed to next script.
6) COMMIT_filtering.m scripts breaks merged fibre bundle into 4 sub-bundles - 2 contra- and 2 ipsilaterally going bundles and does exactly the same what is done in step 5) or LiFE analysis
7) Plots are done in step 6) of LiFE analysis
--------------
