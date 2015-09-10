function [probtfgene, data, data1] =  promProbTFOEsampling(model,expression,expressionid,regulator,regulated,TFsamples,litevidence,prob_prior,datathreshoff,datathreshOE,probtfgene,DATATHRESHVAL,DATATHRESHVALOE,nboots,nsamples)

%% Shuyi Ma 2015
%% This function calculates the conditional probability of Target Gene = ON when TF = OVEREXPRESSED by sampling the expression data.
%% Key Variables:
%% nboots = number of sampling iterations
%% TFsamples = a cell array containing the column numbers of the samples belonging to each group. 
%%             In our work, we had input data comprising multiple replicates of overexpression data from each TF. 
%%             Therefore, each cell of TFsamples contained the column numbers of the matrix "expression" corresponding 
%%             to the overexpression data for a different TF.
%%             If the user has an expression set that should be uniformly sampled without grouping, 
%%             then TFsamples should be a 1 x 1 cell array containing all of the column numbers.
%% nsamples = minimum number of samples to include from each group in TFsamples

%% INPUT HANDLING
%===========================================================
if nargin == 6, litevidence = [];prob_prior = [];
elseif (nargin < 5) || (nargin == 14),
    error('Incorrect number of input arguments to %s',mfilename)
end
%===========================================================
%SOME BASIC INITIALIZATION
%===========================================================
disp('initializing data')

% set default value
if (~exist('DATATHRESHVAL','var')) || (isempty(DATATHRESHVAL))
    DATATHRESHVAL = 0.33; 
end

if (~exist('DATATHRESHVALOE','var')) || (isempty(DATATHRESHVALOE))
    DATATHRESHVAL = 0.9; 
end

if (~exist('nboots','var')) | (isempty(nboots))
    nboots = 500;
end

if (~exist('nsamples','var')) | (isempty(nsamples))
    nsamples = 3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

litevidence = logical(litevidence);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
lost_xn = false(size(regulated));
disp('finding probabilities')

data = expression;
data = knnimpute(data);
data = quantilenorm(data);  %its already normalized..

data1 = data; %data will be binned into 3 categories, data1 is the original values

if isempty(datathreshoff)
    datathreshoff = quantile(data(:),DATATHRESHVAL);
end

if isempty(datathreshOE)
    datathreshOE = quantile(data(:),1-DATATHRESHVAL);
end

if datathreshoff < 0,
    data(data1>=datathreshoff) = 1;
    data(data1 < datathreshoff) = 0;
else
    data(data1 < datathreshoff) = 0;
    data(data1>=datathreshoff) = 1;
end


data(data1 >= datathreshOE) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To reduce the calculation burden, will calculate only the conditional probabilities for the TF-metabolic gene interactions.
% This is done in the wrapper function now.
%mettarints=cellfun(@(x) any(strcmp(x,model.genes)),regulated);
%
%metregulated = regulated(mettarints);
%metregulator = regulator(mettarints);

if (~exist('probtfgene','var')) || (isempty(probtfgene))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     try probtfgene = ones(size(regulated));
    try probtfgene = ones(size(regulated,1),nboots);
    catch M1
        probtfgene = 1; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = cellfun(@(x) find(strcmp(x,expressionid)),regulated,'UniformOutput',false);
    l = cellfun(@(x) find(strcmp(x,expressionid)),regulator,'UniformOutput',false);
    for n = 1:nboots
        bootsamples = [];
        for utf = 1:size(TFsamples,1)
            try
                bootsamples = [bootsamples;randsample(TFsamples{utf},nsamples)];
            catch M1
                bootsamples = [bootsamples;TFsamples{utf}];
            end
        end
        
        datatmp = data1(:,bootsamples);
        databintmp = data(:,bootsamples);
        for  i = 1:length(regulated)
%             k = find(strcmp(regulated(i),expressionid));
%             l = find(strcmp(regulator(i),expressionid));
            if ~isempty(k{i}) & ~isempty(l{i})
                tarexp = datatmp(k{i},:);

                tarbin = databintmp(k{i},:);
                tfbin = databintmp(l{i},:);

                try ksh = kstest2(tarexp(tfbin == 2),tarexp(tfbin== 0));

                    if  ksh

                        prob1 = sum(tarbin(tfbin == 2)>0)/length(tarbin(tfbin==2));
                            probtfgene(i,n) = prob1;
                        %   this formula also gives the same answer  - (sum(~tfbin(tarbin == 1))/length(tfbin(tarbin==1))) * (sum(tarbin)/sum(~tfbin))

                    end

                catch ERRLG    % cant be estimated from microarray ; if it has strong evidence, i might consider setting this to zero later on
    %                 probtfgene(i) = 1;
                    lost_xn(i) = 1;

                end
            else
    %             probtfgene(i) = 1;
                lost_xn(i) = 1;
            end
        end
    end
    
    toc
    
    if ~isempty(litevidence)
        probtfgene(litevidence) = prob_prior(litevidence);  % u could set those interactions that u think have strong literature evidence to have predefined
    end                                                            % probabilities
    
    toc
end