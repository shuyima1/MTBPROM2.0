# MTBPROM2.0


This repository contains codes and input files to run transcription factor knockout and overexpression simulations with an updated version of the Probabilistic Regulation of Metabolism (PROM) framework [Chandrasekaran et al., PNAS 2010].

The inputs include:
1. A genome-scale flux balance-based metabolic model compatible with the COBRA toolbox
2. A two-column list of transcriptional regulatory network interactions
3. A gene expression dataset

*Note:* The input files included in this repository are the components of an updated regulatory-metabolic model for Mycobacterium tuberculosis (MTB), featured in [Ma et al., accepted at PLoS Computational Biology 2015]. The genome scale metabolic model is MTB iSM810, the transcriptional regulatory network interations are TF-metabolic gene binding interactions from ChIP-seq data described in [Minch et al., Nature Communications 2015], and the gene expression dataset is the MTB strain H37Rv TF overexpression dataset described in [Rustad et al., Genome Biology 2014].

## Code Contents:

Simulations in which the conditional probability is calculated from the entire input expression dataset are provided in:
`promv2.m` (knockout simulations), and `promv2TFOE.m` (overexpression simulations).

Code for simulating where the conditional probability is calculated with sampling in `promSampling.m`. 

Code for simulating TF-metabolic gene double knockouts is also provided in `PROMdoubleKO.m`. 

Additionally, code describing how to generate the carbon and nitrogen source models is included in `Metabolite-condition-specific-ModelCode.m`.
