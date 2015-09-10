function [fko] = PROMdoubleKO(model,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,KAPPA,datathresh,DATATHRESHVAL,probtfgene,sizeflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Shuyi Ma 2015
%
% This code predicts the growth phenotype and the flux response
% after a TF-metabolic gene double-knockout. 
% To reduce computational time, the conditional probability of TF influence 
% is estimate from the entire expression dataset rather than by sampling. 
%
% Inputs:
% model - genome-scale metabolic tools from COBRA toolbox (Schellenberger et al. 2011)
% expression - matrix of expression data values. Rows = genes, columns = samples
% expressionid - cell array of gene names corresponding to the rows of the 'expression' matrix
% regulator,regulated - cell arrays listing the TF-target gene interactions. 
%			   Regulator and Regulated have the same number of rows; 
%			   each row describes a separate interaction. For each row, 
%			   regulator(i) lists the TF regulator of that interaction, 
% 			   and regulated(i) lists target gene regulated in that interaction.
% sizeflag  - tells PROM if the regulatory network is large. It is 0 for large networks
%  			   and 1 for small networks ( less than 1000 interactions)
% Optional Inputs:
% litevidence,prob_prior - these should have the same length as the regulator/target arrays;
%			   high confidence interactions (not necessarily based on literature) 
%			   should be flagged as "1" in litevidence array, and the rest should be set as 0. 
%			   User can set the desired conditional probability in the corresponding row of prob_prior
% subsets - cell array denoting subsets of TFs to simulate
% v11,v12 - vmin, vmax derived from Flux Variability Analysis
% KAPPA - determines strength of regulatory constraints - default value 1 works for most systems
% datathresh - expression threshold value for ON vs. OFF gene expression
% DATATHRESHVAL - If user does not want to set an explicit expression threshold value for ON/OFF, can set a quantile value instead. Default = 0.33
% probtfgene - vector of conditional probabilities. Can be input directly if calculated from previous simulation.
% 
%
% Outputs:
% f - matrix of simulated growth rates from the TF perturbations. rows = TFs, columns = sampling iterations;
% probtfgene - matrix of conditional probabilities calculated across sampling iterations.
% metregulator,metregulated - arrays describing the interactions between TFs and target genes 
%			   featured in the metabolic model. These are the only genes that the 
%   		   PROM model will generate a solution. Those TFs with no interactions to any 
%			   metabolic genes would yield a predicted growth rate equal to the wild-type rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numtf = size(unique(regulator),1);
fko = nan(size(model.genes,1),numtf);

[grRatio grKO] = singleGeneDeletion(model);
nonlethal = find(grRatio > 0.01);

fko(grRatio <= 0.01,:) = repmat(grKO(grRatio <= 0.01),1,numtf);
clear grKO

for i = 1:size(nonlethal,1)
    [modelDel] = deleteModelGenes(model,model.genes{nonlethal(i)});
    [f] =  promv2(modelDel,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,[],[],KAPPA,datathresh,DATATHRESHVAL,probtfgene,sizeflag);
    fko(nonlethal(i),:) = f;
end

