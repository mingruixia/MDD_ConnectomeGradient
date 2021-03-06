# MDD_ConnectomeGradient
This repository provides core code and relevant toolboxes for data analysis in the article entitled "Connectome Gradient Dysfunction in Major Depression and Its Association with Gene Expression Profiles " by Xia et al. 2020

## Overview
Content includes standalone software, source code, and demo data. Due to the large size of the analyzed data, we can only provide a small portion of the data needed for validating the code. 
The project is structured into four parts corresponding to the major analyses in the article, including fMRI data preprocessing, gradient analysis, cognitive terms, gene expression association analysis. 

## Toolboxes
All custom code and toolboxes were tested on two 64-bit Windows 10 PCs (PC1: Intel Core i7-6700k, 64GB RAM; PC2 AMD Ryzen Threadripper 3970x, 256G RAM) with MATLAB R2019b, which were included below. 

SPM12, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

SeeCAT, a custom developed toolbox. https://github.com/mingruixia/MDD_ConnectomeGradient/tree/main/0.Preprocessing/SeeCAT

MICA diffusion_map_embedding, ver. 20180921, https://github.com/MICA-MNI/micaopen/tree/master/diffusion_map_embedding

ComBatHarmonization, ver. 20180620, https://github.com/Jfortin1/ComBatHarmonization

AHBAprocessing, ver. 20181025, https://github.com/BMHLab/AHBAprocessing

Scripts from Whitaker et al. PNAS 2016, https://github.com/KirstieJane/NSPN_WhitakerVertes_PNAS2016

Brainnet Viewer, ver. 1.7, http://www.nitrc.org/projects/bnv/

Other software and web-based tools include:

IBM SPSS Statistics 21, commercial software, https://www.ibm.com/products/spss-statistics

NeuroSynth, web-based decoding functions were used, https://neurosynth.org/

Gorilla, http://cbl-gorilla.cs.technion.ac.il/

BrainSMASH, https://brainsmash.readthedocs.io/en/latest/index.html

We thank the authors and developers for providing these wonderful tools for data analysis. 

## Installation guide
Please use the “add path” method in MATLAB to add toolboxes and scripts. This procedure is not time-consuming. 

## Demo
We strongly suggest reading the manual for each toolbox for detailed instructions on how to use them. Here, we provide a brief description of the data analysis in our paper. 

### Data preprocessing
We used SeeCAT and SPM to perform resting-state fMRI data preprocessing. SeeCAT is a GUI-based toolbox for resting-state fMRI data preprocessing, functional connectivity analysis, voxel-based degree calculation, and statistical analyses. 
1. Arrange image files to have the same structure in the demo folder (FunImg for initial Nifti files).
2. After installing SeeCAT and SPM in MATLAB, call the SeeCAT data preprocessing module by typing SeeCAT_PrepfMRI in the MATLAB command window.
3. Locate the data path and start folder, select preprocess steps and click the run button in the GUI.
4. The final output of preprocessed fMRI data in our analysis is stored in FunImgARWSDCFB (A for slice timing, R for realigning, W for normalization, S for smooth, D for detrending, C for covariate regression, F for filtering, and B for scrubbing). 

## Gradient analysis
This part was primarily carried out using our custom script for MATLAB with some functions fulfilled by MICA diffusion_map_embedding (https://github.com/MICA-MNI/micaopen/tree/master/diffusion_map_embedding) and ComBatHarmonization (https://github.com/Jfortin1/ComBatHarmonization). Please see the comments in Script.m for details. The analysis details include:
1. Downsample preprocessed fMRI data to a 4-mm isotropic resolution. The input for this step is the preprocessed fMRI data in the folder FunImgARWSDCFB
2. Generate the voxel-wise connectivity matrix and calculate the connectome gradient.
3. Align gradient maps across individuals by calling the function mica_iterativeAlignment.
4. Calculate the explained variance of each gradient for each individual. 
5. Correct the center effect of the gradient score by using ComBatHarmonization.
6. Calculate gradient range and spatial variance.
7. Generate Nifti files for each gradient map.
8. Generate group-averaged gradient maps for HC and MDD
9. Estimate spatial correlation of the group-averaged map between HC and MDD for each gradient.

Due to the limited upload data size available on Github, we cannot provide all the data in this analysis. The preprocessed data of one subject was provided for test matrix generation and gradient calculation (In the folder of the previous step). We also provided the global measures (i.e., explained variance, gradient range, and spatial variance, GlobalMetrics.sav), final gradient maps of all subjects, and the displacement maps of patients with MDD (check DownloadForGradientMaps.txt for download link). The statistical analysis for global metrics was done by using SPSS and that for gradient maps was done by SeeCAT (call SeeCAT_Stat and SeeCAT_Viewer in MATLAB), respectively. Both tools are GUI-based. 

## Cognitive terms
We used Neurosynth (https://neurosynth.org/) to assess the topic terms associated with the alterations in the connectome gradient in MDD. 
1. The thresholded Z-maps derived from the between-group comparisons for each gradient were first divided into MDD-positive and MDD-negative maps (see Z-maps in demo data). 
2. The Z-maps were then uploaded to Neurovault and analyzed using the “decoder” function of the Neurosynth website. 
3. For each of the maps, the noncognitive terms (e.g., anatomical and demographic terms) were removed and the top 30 cognitive terms were selected. 
4. The cognitive terms were visualized on a word-cloud plot with the font size scaled according to their correlation with corresponding meta-analytic maps generated by Neurosynth.

## Gene expression association
In this analysis, we used the revised script from Whitaker et al. PNAS 2016, (https://github.com/KirstieJane/NSPN_WhitakerVertes_PNAS2016). The main script is named gene_PLS.m and see the code in this file for more details. 
1. The Gene expression data from the Allen Institute for Brain Science was first preprocessed by using AHBAprocessing (https://github.com/BMHLab/AHBAprocessing), obtaining the gene expression profile for the Glasser-360 atlas (100DS360scaledRobustSigmoidNSGRNAseqQC1LRcortex_ROI_NOdistCorrEuclidean.mat, 4mm_Glasser360.nii).
2. Perform PLS analysis to examine the association between the between-group Z-map (Z-maps.zip) of the connectome gradients and gene expression profiles. 
3. Significance of the PLS component was determined by permutation test in which the spatial autocorrelations were corrected by generative modeling (gen_surrogate_map_for g1z.py)
4. The weights of PLS components were determined by the bootstrap method . 
5. Both the descending order and ascending order of PLS weighted genes were submitted to GOrilla (http://cbl-gorilla.cs.technion.ac.il/) for enrichment analysis. See the GOrilla website for instructions. 

## Validation
We validated our results by considering several potential confounding factors. First, we used a leave-one-site-out cross-validation strategy to examine whether our findings were influenced by specific sites. This was implemented by repeating the between-group comparisons on the data, excluding one site at a time. Second, some of the participants were younger than 18 years, which might explain the between-group differences in brain development. Thus, we repeated the statistical analysis for only adult participants (1,002 patients with MDD and 1,034 HCs). Third, to further control for the effect of head motion on R-fMRI connectivity measures, we repeated the between-group comparisons with the mean framewise displacement as an additional covariate. Finally, the use of a joint embedding framework may increase the alignment of individual connectome gradient maps compared to the Procrustes rotation. Thus, we reperformed the alignment using joint embedding and then repeated the statistical analysis.
