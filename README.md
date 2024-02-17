# SC-FC-prediction
MATLAB implementations of two measures to statistically evaluate the extent to which functional connectivity (FC) predictions derived from an individual's structural connectome (SC) capture individual-specific effects. 

# Context
Several recent studies have optimized deep neural networks to learn high-dimensional relationships linking structural and functional connectivity across the human connectome. However, the extent to which these models recapitulate individual-specific characteristics of resting-state functional brain networks remains unclear. We provide two measures and associated benchmark/null conditions to evaluate the extend to which an individual's predicted FC captures individual effects. 

For further details, see Zalesky et al, Predicting an individual’s functional connectivity from their structural connectome: Evaluation of evidence, recommendations and future prospects. 

# Notation
eFC: empirically measured FC

pFC: predicted FC

SC: structural connectome, derived from diffusion MRI tractography

# Measure 1: Intra- vs inter-individual similarity 
A similiarity matrix is first computed to quantify the correlation between an individual's predicted FC (pFC) and the empirically measured FC (eFC) of all other individuals. This yields an N x N matrix, where N is the number of individuals and element (i,j) quantifies the correlation between individual i's pFC and individual j's eFC. Statistical testing is then performed to determine whether intra-individual similarity is greater than inter-individual similarity. This is similar to the measure proposed by Chen et al (2023). Similarity can be alternatively measured with a Riemannian distance such as the log-Euclidean Riemannian metric (LERM). 

# Measure 2: Individual matching
A binary assignment problem is solved in which individual pFC matrices are matched to individual eFC matrices. The Hungarian algorithm is used to compute an optimal individual matching. The proportion of individuals correctly matched (i.e. pFC-eFC matching) is compared to a null condition. Exceedance of the null condition suggests that individual pFC matrices capture individual effects beyond chance levels.  

# Data
SC, eFC and pFC matrices for 1000 individuals comprising the Human Connectome Project (HCP) are stored in full_data_with_pFCs.mat 

This data file contains derivatives from diffusion and functional MRI data made available by the HCP via ConnectomeDB (www.humanconnectome.org) and it is not made available here. 

![image](https://github.com/AndrewZalesky/SC-FC-prediction/assets/57614210/91de1ce9-5706-4863-acab-e8f2f11979a2)
