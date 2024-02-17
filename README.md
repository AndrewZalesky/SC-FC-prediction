# SC-FC-prediction
Two measures to statistically evaluate the extent to which predictions of a functional connectivity (FC) matrix capture individual-specific effects

# Measure 1: Intra- vs inter-individual similarity 
A similiarity matrix is first computed to quantify the correlation between and individual's predicted FC (pFC) and the empirically measured FC (eFC) of all other individuals. This yields an N x N matrix, where N is the number of individuals and element (i,j) quantifies the correlation between individual i's pFC and individual j's eFC. Statistical testing is then performed to determine whether intra-individual similarity is greater than inter-individual similarity. 

# Measure 2: Individual matching

