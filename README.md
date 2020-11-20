# MDD_ConnectomeGradient
Core code and relevant toolboxes used for data analysis in the paper titled "Large-scale Gradient Dysfunction of the Functional Connectome in Major Depression" by Xia et al. 2020

## Overview
The content includes complied standalone software, source code, and demo data. Due to the large size of the analyzed data, we can only provide a small portion of data for validating the code. 
The project is structured into five parts corresponding to the major analyses in the paper, including fMRI data preprocessing, gradient analysis, cognitive terms, gene expression association analysis, and symptom predication analysis. 

## Toolboxes
All custom code and toolboxes were tested on a 64-bit Windows 10 PC (i7-6700k, 64GB RAM) with MATLAB R2019b, including

SPM12, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

SeeCAT, a custom developed toolbox, available in this project

MICA diffusion_map_embedding, ver. 20180921, https://github.com/MICA-MNI/micaopen/tree/master/diffusion_map_embedding

ComBatHarmonization, ver. 20180620, https://github.com/Jfortin1/ComBatHarmonization

AHBAprocessing, ver. 20181025, https://github.com/BMHLab/AHBAprocessing

Scripts from Whitaker et al. PNAS 2016, https://github.com/KirstieJane/NSPN_WhitakerVertes_PNAS2016

libsvm, ver. 3.24, https://www.csie.ntu.edu.tw/~cjlin/libsvm/

Scripts from Cui et al. Neuroimage 2018, https://github.com/ZaixuCui/Pattern_Regression_Matlab

Brainnet Viewer, ver. 1.7, http://www.nitrc.org/projects/bnv/


Other softwares and web-based tools include

IBM SPSS Statistics 21, commercial software, https://www.ibm.com/products/spss-statistics

NeuroSynth, web-based decoding functions were used, https://neurosynth.org/

Gorilla, http://cbl-gorilla.cs.technion.ac.il/

We thank the authors and developers for providing these wonderful tools for data analysis. 

## Installation guide
Please use “add path” in MATLAB to add toolboxes and scripts. This procedure is not time consuming. 

## Demo
We fully suggest to read the manual of each toolbox for detailed instruction of how to use them. Here, we provide a brief description of data analysis in our paper. 

### Data preprocessing
We use SeeCAT and SPM to perform resting-state fMRI data preprocessing. SeeCAT is a GUI-based toolbox for resting-state fMRI data preprocessing, functional connectivity analysis, voxel-based degree calculation, and statistical analysis. 
1. Arrange image files as the same structure in the demo folder (FunImg for initial Nifti files).
2. After installing SeeCAT and SPM in MATLAB, call SeeCAT data preprocessing module by typing SeeCAT_PrepfMRI in the MATLAB command window.
3. Locate the data path and start folder, select preprocess steps, and click run button.
4. The final output preprocessed fMRI data in our analysis is stored in FunImgARWSDCFB (A for slice timing, R for realign, W for normalization, S for smooth, D for detrend, C for covariate regression, F for filtering, and B for scrubbing). 

## Gradient analysis
This part was mainly based on our custom script for MATLAB with some functions fulfilled by MICA diffusion_map_embbeding (https://github.com/MICA-MNI/micaopen/tree/master/diffusion_map_embedding) and ComBatHarmonization (https://github.com/Jfortin1/ComBatHarmonization). Please see comments in Script.m for details. The analysis details included
1. Downsample preprocessed fMRI data to a 4-mm isotropic resolution. The input of this step is the preprocessed fMRI data in folder FunImgARWSDCFB
2. Generate voxel-wise connectivity matrix and calculate connectome gradient.
3. Align gradient maps across individuals by calling the function mica_iterativeAlignment.
4. Calculate the explained variance of each gradient for each individual. 
5. Correct the center effect of gradient score by using ComBatHarmonization.
6. Calculate gradient range and spatial variance.
7. Generate nifti files for each gradient map.
8. Calculate displacement of each voxel in gradient space for each MDD. 

Due to the limited upload data size of Github, we cannot provide all data in this analysis. The preprocessed data of one subject was provided for testing matrix generation and gradient calculation (In the folder of previous step). We also provided the global measures (i.e., explained variance, gradient range, and spatial variance, GlobalMetrics.sav), final gradient maps of all subjects, and the displacement maps of patients with MDD (check DownloadForGradientMaps.txt for download link). The statistical analysis for global metrics was done by using SPSS and that for gradient maps was done by SeeCAT (call SeeCAT_Stat and SeeCAT_Viewer in MATLAB), respectively. Both tools are GUI-based. 

## Cognitive terms
We used Neurosynth (https://neurosynth.org/) to assess the topic terms associated with the alterations in the connectome gradient in MDD. 
1. The thresholded Z-maps derived from the between-group comparisons for each gradient were first divided into MDD-positive and MDD-negative maps (see Z-maps in demodata). 
2. The Z-maps were then uploaded to Neurovault and analyzed using the “decoder” function of Neurosynth website. 
3. For each of the maps, the noncognitive terms (e.g., anatomical and demographic terms) were removed and the top 30 cognitive terms were selected. 
4. The cognitive terms were visualized on a word-cloud plot with the font size scaled according to their correlation with corresponding meta-analytic maps generated by Neurosynth.

## Gene expression association
In this analysis, we mainly use the revised script from Whitaker et al. PNAS 2016, (https://github.com/KirstieJane/NSPN_WhitakerVertes_PNAS2016). The main script is named pls_gradient_geneEx.m and see code in this file for more details. 
1. The Gene expression data from the Allen Institute for Brain Science was first preprocessed by using AHBAprocessing (https://github.com/BMHLab/AHBAprocessing), obtaining the gene expression profile for the Glasser-360 atlas (100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat, 4mm_Glasser360.nii).
2. Perform PLS to examine the association between the between-group Z-maps (Z-maps.zip) of the connectome gradients and gene expression profiles. 
3. The weight of PLS components were determined by bootstrap method (PLS_bootstrap_Jin.m) and their significance was determined by permutation test (PLS_calculate_stats_Jin). 
4. Both the descending order and ascending order of PLS weighted genes were submitted to GOrilla (http://cbl-gorilla.cs.technion.ac.il/) for enrichment analysis. See GOrilla website for instruction. 

## SVR-based symptom prediction
We used support vector regression (SVR) and 10-fold cross-validation (where the predictive model was repeatedly trained on 9 folds of the data and tested on the 10th fold) to examine whether the connectome gradient features were able to predict depressive symptoms (Hamilton Depression Rating Scale, HDRS) in the patients. This was done by calling script from Cui et al. Neuroimage 2018 (https://github.com/ZaixuCui/Pattern_Regression_Matlab) and libsvm 3.24 (https://www.csie.ntu.edu.tw/~cjlin/libsvm/). The displacement in the gradient space of the clusters with MDD-related alterations as features. See Script.m for details. The demo data included the extracted features (SVR_Subjects_Data.mat), the HDRS (SVR_Subjects_Scores.mat), and the covariates (SVR_Covariates.mat).

## Reproducibility
We estimated the reproducibility of the identified MDD-related gradient alterations by considering several potential confounding factors. First, we used a leave-one-site-out cross- validation strategy repeating the between-group comparisons on the data, excluding one site at a time. Second, some of the participants were younger than 18 years old, which may explain the between-group difference in brain development. Thus, we reperformed the data analysis for only adult participants. Finally, we repeated the between-group comparisons with the mean framewise displacement as an additional covariate.
These analyses were mainly done via re-arrange data or add additional covariates in statistical analysis. Thus, no further code is required.
