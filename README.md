# Network Controllability in Paediatric Drug Resistant Epilepsy

This repository contains the instructions & code to construct structural connectivity matrices and calculate network controllability and covers the analyses in the follwoing manuscript:

*Drug-resistant focal epilepsy in children is associated with increased modal controllability of the whole brain and epileptogenic regions* by Chari et al. 

# Image Processing

## Structural Imaging & Parcellation

The structural imaging (3T T1 MPRAGE sequences) were organised in [BIDS](https://bids.neuroimaging.io) format and fed through [Connectome Mapper 3](https://connectome-mapper-3.readthedocs.io/en/latest/). We chose the Lausanne Atlas at Scale 3 for this study. 

## Diffusion Imaging Preprocessing

Multi-shell multi-tissue diffusion imaging was preprocessed using standard preprocessing steps according to the script **diffusion_preproc.sh**. The multi-shell multi-tissue dicom files were stored in a folder called 'dwi' and the negative phase encoded series in a folder called 'negPE'. All tools were [mrtrix](https://mrtrix.readthedocs.io/en/latest/) in line with their recommended steps. 

## Constructing the Connectome

This was done by combining the parcellation scheme with the preprocessed diffusion imaging according to the script **make_connectome.sh**. Firstly the structural data was registered to diffusion space and the cortical target labels were amended according to **Scale3_NewLabels.txt**. For this, we fused the hippocampal subparcels into a single label, resulting in 253 cortical, subcortical and brainstem parcels per subject. The connectome used anatomically constrained tractography and the SIFT2 filtering step. The output is a 253x253 weighted, undirected structural connectivity matrix in .csv format.

# Connectome Analyses

## Organising Data

The extracted connectomes were assembled as structures in Matlab v2020b using the scripts **ExtractConnectomes.m** and **LoadConnectomes.m**. Following this, the weighted degree, average and modal controllabilities for each node in each patient were calculated using the script **CalculateControllabilities.m**. These final data files with the connectivity matrices and controllabilities for each subject is avaiable within each group folder (Controls, Patients, PostopPatients and VNS).

## Main Analyses

The main analyses are encompassed by the file **MainAnalyses.m**. This file calculates the ranks, correlation coefficients (WD-AC, WD-MC and AC-MC), Z-scores etc that encompass all the patient data analysis elements of the manuscript. The **GraphTheory.m** file follows on from this and contains the script that calcualtes the graph theory metrics. It requires the [BrainConnectivityToolbox](https://sites.google.com/site/bctnet/) to be installed. 

GLMs were conducted by manually exporting the data to SPSS v24. 

# Modelling and Simulation

The modelling elements are contained in the folder **Modelling**. The code for first modelling exercise, adding random edges to the controls, is contained in **ModellingEx1.m**. The code for the second modelling exercise of adding low edge weights to thalamocortical edges is contained in **ModellingEx2.m**. 

# Visualisation

Elements required for visualising the atlas are contained in the folder **SurfIce** and requires dowloading [SurfIce](https://www.nitrc.org/projects/surfice/) and the [AtlasStatMap](https://github.com/rordenlab/spmScripts/blob/master/AtlasStatMap.m) function. The folder also contains a version of the atlas in SurfIce mz3 format (merge.mz3). Use the code included to display any assigned values to the parcels. For more information, please see the [SurfIce Wiki on viewing atlases](https://www.nitrc.org/plugins/mwiki/index.php/surfice:MainPage#Atlas-based_Region_of_interest_Analyses). 
