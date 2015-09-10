function [f,probtfgene,metregulator,metregulated] =  promSampling(model,expression,expressionid,regulator,regulated,TFsamples,TFstatus,litevidence,prob_prior,subsets,v11,v12,KAPPA,datathresh,DATATHRESHVAL,probtfgene,sizeflag,nboots,nsamples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Shuyi Ma 2015
%
% This code predicts the growth phenotype and the flux response
% after transcriptional perturbation (either Overexpression or Knockout) after 
% calculating conditional probabilities of TF influence based on 
% multiple iterations of random sampling the gene expression data. 
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
% TFsamples- a cell array containing the column numbers of the samples belonging to each group. 
%%             In our work, we had input data comprising multiple replicates of overexpression data from each TF. 
%%             Therefore, each cell of TFsamples contained the column numbers of the matrix "expression" corresponding 
%%             to the overexpression data for a different TF.
%%             If the user has an expression set that should be uniformly sampled without grouping, 
%%             then TFsamples should be a 1 x 1 cell array containing all of the column numbers.
% TFstatus -  a string denoting whether knockout (enter 'TFKO') or an overexpression (enter 'TFOE') simulation is desired. 
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
% nboots - number of sampling iterations
% nsamples - minimum number of samples use from each group in TFsamples.
%
% Outputs:
% f - matrix of simulated growth rates from the TF perturbations. rows = TFs, columns = sampling iterations;
% probtfgene - matrix of conditional probabilities calculated across sampling iterations.
% metregulator,metregulated - arrays describing the interactions between TFs and target genes 
%			   featured in the metabolic model. These are the only genes that the 
%   		   PROM model will generate a solution. Those TFs with no interactions to any 
%			   metabolic genes would yield a predicted growth rate equal to the wild-type rate.
%
%
% example call for simulating TF knockouts with sampled conditional probabilities:
% [f,probtfgene,metregulator,metregulated] =  promSampling(iSM810_7H9,expression,expressionid,MinchChIPseqOperonTFints(:,1),MinchChIPseqOperonTFints(:,2),TFsamples,'TFKO',[],[],[],[],[],[],[],[],[],0,500,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter Initialization
if nargin == 7
	litevidence = [];
	prob_prior = [];
	subsets = [];
	v11 = [];
	v12 = [];
	KAPPA = 1;
	datathresh = [];
	DATATHRESHVAL = 0.33;
	probtfgene = [];
	sizeflag = 0;
	nboots = 500;
	nsamples = 3;
elseif (nargin < 7) | (nargin > 20),
    error('Incorrect number of input arguments to %s',mfilename)
end

% To reduce the calculation burden, will calculate only the conditional probabilities for the TF-metabolic gene interactions.
mettarints=cellfun(@(x) any(strcmp(x,model.genes)),regulated);

metregulated = regulated(mettarints);
metregulator = regulator(mettarints);
% metTFs = unique(metregulator);

if isequal(upper(TFstatus),'TFKO')
	[probtfgene] =  promProbTFsampling(model,expression,expressionid,metregulator,metregulated,TFsamples,litevidence,prob_prior,datathresh,probtfgene,DATATHRESHVAL,nboots,nsamples);

	f = zeros(size(unique(metregulator),1),nboots);

	for i = 1:nboots
    	f(:,i) = promv2(model,expression,expressionid,metregulator,metregulated,[],[],subsets,v11,v12,KAPPA,datathresh,DATATHRESHVAL,probtfgene(:,i),sizeflag);
	end
elseif isequal(upper(TFstatus),'TFOE')
	[probtfgene] =  promProbTFOEsampling(model,expression,expressionid,metregulator,metregulated,TFsamples,litevidence,prob_prior,datathresh,probtfgene,DATATHRESHVAL,DATATHRESHVALOE,nboots,nsamples);

	f = zeros(size(unique(regulator),1),nboots);

	for i = 1:nboots
    	f(:,i) = promv2TFOE(model,expression,expressionid,metregulator,metregulated,[],[],subsets,v11,v12,KAPPA,datathresh,DATATHRESHVAL,DATATHRESHVALOE,probtfgene(:,i),sizeflag);
	end
end