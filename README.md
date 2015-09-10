# MTBPROM2.0


This repository contains codes to run transcription factor knockout and overexpression simulations with an updated version of the Probabilistic Regulation of Metabolism (PROM) framework (Chandrasekaran et al., PNAS 2010).

The inputs include:
1. A genome-scale flux balance-based metabolic model compatible with the COBRA toolbox,
2. A two-column list of transcriptional regulatory network interactions, and 
3. A gene expression dataset.

 
Simulations in which the conditional probability is calculated from the entire input expression dataset are provided in:
promv2.m (knockout simulations), and promv2TFOE.m (overexpression simulations).

Code for simulating where the conditional probability is calculated with sampling (promSampling.m). 

Code for simulating TF-metabolic gene double knockouts is also provided (PROMdoubleKO.m). 

Additionally, code describing how to generate the carbon and nitrogen source models is included (Metabolite-condition-specific-ModelCode.m).
