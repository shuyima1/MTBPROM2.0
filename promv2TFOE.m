function [f,f_ko,v,v_ko,status1,lostxns,probtfgene,bnumstobekoed] =  promv2TFOE(model,expression,expressionid,regulator,regulated,litevidence,prob_prior,subsets,v11,v12,KAPPA,datathreshoff,datathreshOE,DATATHRESHVAL,DATATHRESHVALOE,probtfgene,sizeflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Shuyi Ma 2015
% Adapted based on PROM code from Chandrasekaran et al., PNAS (2010) 
%
% This code predicts the growth phenotype and the flux response
% after TF overexpression, given a metabolic and regulatory
% network. Conditional probabilities are calculated based on the entire expression matrix, not sampled.
%
% Inputs:
% model - genome-scale metabolic tools from COBRA toolbox (Schellenberger et al. 2011)
% expression - matrix of expression data values. Rows = genes, columns = samples
% expressionid - cell array of gene names corresponding to the rows of the 'expression' matrix
% regulator,regulated - cell arrays listing the TF-target gene interactions. 
%              Regulator and Regulated have the same number of rows; 
%              each row describes a separate interaction. For each row, 
%              regulator(i) lists the TF regulator of that interaction, 
%              and regulated(i) lists target gene regulated in that interaction.
% sizeflag  - tells PROM if the regulatory network is large. It is 0 for large networks
%              and 1 for small networks ( less than 1000 interactions)
% Optional Inputs:
% litevidence,prob_prior - these should have the same length as the regulator/target arrays;
%              high confidence interactions (not necessarily based on literature) 
%              should be flagged as "1" in litevidence array, and the rest should be set as 0. 
%              User can set the desired conditional probability in the corresponding row of prob_prior
% subsets - cell array denoting subsets of TFs to simulate
% v11,v12 - vmin, vmax derived from Flux Variability Analysis
% KAPPA - determines strength of regulatory constraints - default value 1 works for most systems
% datathreshoff - expression threshold value for ON vs. OFF gene expression
% datathreshOE - expression threshold value for ON vs. OVERexpressed gene expression
% DATATHRESHVAL - If user does not want to set an explicit expression threshold value for ON/OFF, can set a quantile value instead. Default = 0.33
% DATATHRESHVALOE - If user does not want to set an explicit expression threshold value for ON/OVERexpressed, can set a quantile value instead. Default = 0.9
% probtfgene - vector of conditional probabilities. Can be input directly if calculated from previous simulation.
%
%%
% Outputs: - the algorithm gives the growth rate (f) and flux response (v) after knock out of all
% f - the simulated growth rate of each TF overexpression perturbation.
% v - The flux response profile of each TF overexpression perturbation. 
%
% status - glpk solver status(should be 5 for glpk for correct function; if not then check solver error log)
%
% lostxns - the interactions that could not be quantified based on the
%              threshold set for binarization
%
% probtfgene -  the conditional probabilities of TF influence estimated for each interaction
%
%
% EXAMPLE
% load MTBPROM2_0Inputs
%
%[f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2TFOE(iSM810_7H9,expression,expressionid,MinchChIPseqOperonTFints(:,1),MinchChIPseqOperonTFints(:,2),[],[],[],[],[],[],0.2,0.7,[],[],[],0);
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT HANDLING
%===========================================================
if nargin == 5, litevidence = [];prob_prior = [];subsets = [];
elseif (nargin < 5) || (nargin == 6),
    error('Incorrect number of input arguments to %s',mfilename)
end
%===========================================================
%SOME BASIC INITIALIZATION
%===========================================================
disp('initializing data')

[tfnames] = unique(regulator);

thresh = 10^(-6); mthresh = 10^(-3);

if (~exist('KAPPA','var')) || (isempty(KAPPA))
    KAPPA = 1;% KAPPA = 100;
end

%% for doing flux variability
%if (length(regulator) < 1000) % for small networks its faster to compute min/max for the subset thats regulated than the entire metab model
%    sizeflag = 1;
%else
%    sizeflag = 0;
%end
%%

if ~sizeflag
    if (~exist('v11','var')) || (isempty(v11))
        [v11,v12] = fluxVariability(model); 
        
    end
end

%fprintf('params used - KAPPA: %d and DATATHRESH: %1.3f \n', KAPPA,DATATHRESHVAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

litevidence = logical(litevidence);
if ~isempty(subsets)
    bnumstobekoed = subsets;
    regulated = regulated(ismember(regulator,subsets));
    regulator = regulator(ismember(regulator,subsets));
else
    bnumstobekoed= tfnames;   % bnumstobekoed -  gives the geneids of the genes to be knocked out - by default it knocksout all the tfs in the model one by one
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scou = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('probtfgene','var')) || (isempty(probtfgene))
    [probtfgene lost_xn datathreshflag] = probtfgenecalc(expression,expressionid,regulated,regulator,datathreshoff,datathreshOE,DATATHRESHVAL,DATATHRESHVALOE,litevidence,prob_prior);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  run PROM for each overexpresser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('running PROM')

[rxnpos,genelist] = find(model.rxnGeneMat);
% u - reaction; v - genes
% finding rxn position

% count = 1; 
% allgenes = [model.genes;tfnames];
[posgenelist,posgenelist] = ismember(regulated,model.genes);

weights = model.c; S = model.S; ctype = repmat('S',size(model.b));lbf = model.lb; ubf = model.ub;dxdt = model.b; param.msglev = 1;

lbf(lbf==ubf) = ubf(lbf == ubf) - thresh;

% [v,f_wt] = glpk(-weights,stoic,dxdt,lbf,ubf,ctype);
% f_wt = -f_wt;

[v,f0] = glpk(-weights,S,dxdt,lbf,ubf,ctype);
% flooring small values to zero
v(abs(v) < thresh) = 0;
vm = zeros(size(v));
%vm = nan(size(v));

%% new additions

% lbg = model.lb; ubg = model.ub;
% lbg(lbg==ubg) = ubg(lbg == ubg) - thresh;
a1 = [S,zeros(size(S,1),numel(ubf)),zeros(size(S,1),numel(ubf))];
a2 = sparse([eye(numel(ubf)), eye(numel(ubf)),zeros(numel(ubf))]);
a3 = sparse([eye(numel(ubf)), zeros(numel(ubf)),-eye(numel(ubf))]);
A = [a1;a2;a3];
weights11 = [weights;zeros(2*numel(lbf),1)];
weights00 = [weights;zeros(2*numel(lbf),1)];
lb11 = [-1000*ones(numel(lbf),1);zeros(numel(lbf),1);zeros(numel(lbf),1)];
ub11 = [1000*ones(numel(lbf),1);zeros(numel(lbf),1);zeros(numel(lbf),1)];
lb11(lb11==ub11) = ub11(lb11 == ub11) - thresh;
dxdt0 = [zeros(size(S,1),1);lbf;ubf];
ctype1 = [repmat('S',size(S,1),1);repmat('L',size(lbf,1),1);repmat('U',size(lbf,1),1)];

% [f0,f0] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);

%%



% kapavals = [0,0.0001,0.001,0.05,0.1,0.25,0.33,0.5,1,5,10];
%kappa = kapavals(scou);
% kappa = KAPPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hw = waitbar(0);
% knock out the genes


% bnumstobekoed = [tfs_that_are_off;tf_that u want to knockout];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if ~sizeflag
    % disp('running fast fva for the entire network'); tic;
    %         [v11,v12] = fastFVA(model); toc;
    v11(abs(v11) < thresh) = 0;
    v12(abs(v12) < thresh) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v00 = zeros(length(bnumstobekoed),numel(lb11));
f00 = zeros(length(bnumstobekoed),1);
status1 = zeros(length(bnumstobekoed),1);
f = zeros(length(bnumstobekoed),1);
f_ko = zeros(length(bnumstobekoed),1);
v_ko = zeros(length(bnumstobekoed),numel(lbf));

for ci = 1:length(bnumstobekoed)
    
%     disp(ci)
    
    %% resetting the flux bounds
    lbg = lbf; ubg = ubf;
    lb11 = [-1000*ones(numel(lbg),1);zeros(numel(lbg),1);zeros(numel(lbg),1)];
    ub11 = [1000*ones(numel(lbg),1);zeros(numel(lbg),1);zeros(numel(lbg),1)];
    % check if its a metabolic or regulatory gene or both
    
    if any(strcmpi(model.genes,bnumstobekoed(ci)))
        temppos = rxnpos(genelist == find(strcmp(model.genes,bnumstobekoed(ci))));
        for jj = 1:numel(temppos)
            if model.rev(temppos(jj))
                lbg(temppos) = -thresh;
                ubg(temppos) = thresh;
            else
                lbg(temppos) = -thresh;
            end
        end
        
    end
    
%     [v1,fk(ci)]  = glpk(-weights,S,dxdt,lbg,ubg,ctype);
%     % flooring small values to zero
%     v1(abs(v1) < thresh) = 0;
    
    if any(ismember(tfnames,bnumstobekoed(ci))),
        
        
%         tfstate = false(size(tfnames));
%         tfstate(ismember(tfnames,bnumstobekoed(ci))) = true;
%         tfstate = logical(ismember(tfnames,bnumstobekoed(ci)));
%         k = find(ismember(regulator,tfnames(tfstate)));

        %% This section isolates the interactions associated with the TF being knocked out
        k = ismember(regulator,bnumstobekoed(ci));
        %    k(tempgeneprobs == 1) = '';
        tempgene = regulated(k);
        tempgeneprobs = probtfgene(k);
        tempgenepos = posgenelist(k);
        temprxnpos = rxnpos(ismember(genelist,tempgenepos));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% this section is for gene-protein-reaction relationship
        x = true(size(model.genes));
        [geneInd,geneInd] = ismember(tempgene,model.genes);
%         x(geneInd) = false;
        x(geneInd(geneInd > 0)) = false;
        
        constrainRxn = false(numel(temprxnpos),1);
        % Figure out if any of the reaction states is changed
        for j = 1:numel(temprxnpos)
            if (~eval(model.rules{temprxnpos(j)}))
                constrainRxn(j) = true;
            end
        end
        % Constrain flux through the reactions associated with these genes
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tempgeneprobs(tempgenepos == 0)  = '';
        tempgene(tempgenepos == 0) = '';
        tempgenepos(tempgenepos == 0)  = '';
        
        % temprxnpos has the rxns that are going to  be affected by this tf
        % krxnpos are the rxns that will be affected by this target gene alone..
        % we loop around all the genes..
        
        for l = 1:numel(tempgenepos)
            if ~isnan(tempgeneprobs(l))
                krxnpos = ismember(temprxnpos,rxnpos(ismember(genelist,tempgenepos(l))));
                for m = 1:numel(temprxnpos)
                    if krxnpos(m)
                        if constrainRxn(m)
                            if (tempgeneprobs(l) < 1)    % if its 1 no use in changing the bounds - might as well save time
                                
                                if (tempgeneprobs(l) ~= 0)   % if its zero no point in estimating vm again - saves time.. but cant include in the above statement coz u have to change the bounds
                                    %if v(temprxnpos(m))
                                    if ~(vm(temprxnpos(m)))    % done to save time - if estimated already use it
                                        
                                        if sizeflag
                                            %% redoing FVA with the
                                            % temporary reactions minimized
                                            % and maximized
                                            weights1 = weights; lbv = lbf; ubv = ubf;
                                            grwthpos = find(weights == 1);
                                            lbv(grwthpos) = v(grwthpos);
                                            weights1(temprxnpos(m)) = -1;
                                            [v11,fva1]  = glpk(-weights1,S,dxdt,lbv,ubv,ctype);
                                            weights1(temprxnpos(m)) = 1;
                                            [v12,fva2]  = glpk(-weights1,S,dxdt,lbv,ubv,ctype);
                                        end
                                        
                                        %% modifying the fva flux bounds
                                        if  v(temprxnpos(m)) < 0
                                            vm(temprxnpos(m)) = min([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        elseif v(temprxnpos(m)) > 0
                                            vm(temprxnpos(m)) = max([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        else
                                            vm(temprxnpos(m)) = max([abs(v11(temprxnpos(m))),abs(v12(temprxnpos(m))),abs(v(temprxnpos(m)))]);
                                        end
                                    end
                                    %end
                                end
                                
                                %%
                                xx =    vm(temprxnpos(m))*tempgeneprobs(l); % flux times probability
                                
                                %%
                                if  v(temprxnpos(m)) < 0
                                    
                                    %%Update the flux bounds
                                    
                                    tem = max([lbf(temprxnpos(m)),xx,lbg(temprxnpos(m))]);  %make sure we arent violating the original bounds; also get the lowest value if there were multiple modifications for the rxn
                                    lbg(temprxnpos(m)) = min([tem,-thresh]);   % prevents the solver from crashing
                                    
                                    ub11(1*numel(ubg) + temprxnpos(m)) = 1000;
                                    
%                                     weights11(1*numel(ubg) + temprxnpos(m)) = (-1*KAPPA/abs(vm(temprxnpos(m))))*abs(f0);   % v0 f0 are the wild type values..
%                                     weights11(1*numel(ubg) + temprxnpos(m))
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);
                                    
%                                     weights11(1*numel(ubg) + temprxnpos(m)) = min([(KAPPA*(-1)*abs(f0))/(abs(vv)), weights11(1*numel(ubg) + temprxnpos(m)) ]);
                                    weights11(1*numel(ubg) + temprxnpos(m)) = min([((-1)*KAPPA*abs(f0))/(abs(vv)), (-1*KAPPA*abs(f0)/abs(vm(temprxnpos(m))))]);
                                elseif v(temprxnpos(m)) > 0
                                    
                                    tem = min([xx,ubf(temprxnpos(m)),ubg(temprxnpos(m))]);
                                    ubg(temprxnpos(m)) = max(tem,thresh);
                                    
                                    ub11(2*numel(ubg) + temprxnpos(m)) = 1000;
                                    
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);
                                    
                                    weights11(2*numel(ubg) + temprxnpos(m)) = min([((-1)*KAPPA*abs(f0))/abs(vv), weights11(2*numel(ubg) + temprxnpos(m)) ]);  % new weights based on KAPPA, normalized with growth rate
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    dxdt0 = [zeros(size(S,1),1);lbg;ubg];
%     lb11(lb11==ub11) = ub11(lb11 == ub11) - thresh; % prevents solver from crashing
    param.itlim = 1000000;
    %  optimizeCbModel
    [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
    
    %f(ci) = v00(ci,weights);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%
    coun  = 1;        lbh = lbg; ubh  = ubg;
    
    % clean up section for TFs that did not compute properly
    while ((status1(ci) ~= 5) || (v00(ci,logical(weights)) < 0))
        disp(ci)
        lbh(lbh ~= lbf) = lbh(lbh ~= lbf) - mthresh;
        ubh(ubh ~= ubf) = ubh(ubh ~= ubf) + mthresh;
        dxdt0 = [zeros(size(S,1),1);lbh;ubh];
        [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
        coun = coun + 1;
        if (coun > 50),
            % if its impossible to estimate, then
            % check the unweighted g.r and the one with max weight prom -
            % if very less difference use it - any way warn the user about
            % the problem at the iteration number - ci;
            [v3,f3,status3] = glpk(-weights,S,dxdt,lbg,ubg,ctype);
            [v30,f30,status30] = glpk(-weights00,A,dxdt0,lb11,ub11,ctype1);
            %if abs((f3- v30(find(weights)))/abs(f3)) < 0.1
            f00(ci) = -f3;
            v00(ci,logical(weights)) = -f3;
            %else
            disp(' problem in'); disp(ci);break;  % if that doesnt work,  display a warning
            %end
            
            disp('check'); disp(ci); break;
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    f(ci) = v00(ci,logical(weights));

%     lbg_st(ci,:) = lbg;
%     ubg_st(ci,:) = ubg;
%     
%     lb_st(ci,:) = lbf;
%     ub_st(ci,:) = ubf;
%     
    [v_ko(ci,:),f1] = glpk(-weights,S,dxdt,lbg,ubg,ctype);
    
    f_ko(ci) = -f1;
    
%     ktime = toc;
%     waitpar = [num2str(ceil(ci/numel(tfnames)*100)),'% complete. Time taken:',num2str(ceil(ktime)),' secs'];
    %waitbar(ci/numel(tfnames),hw,waitpar);
    
    
%     ff00(scou,ci) = v00(ci,logical(weights));
    %disp(ff00(scou,ci))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if datathreshflag
%         if all(lost_xn(k))   % if none of the probabilities of a gene can be estimated, then ->
%             disp('Error: xns cant be estimated')
%             v00(ci,:) = NaN;
%             f_ko(ci) = NaN;
%             v_ko(ci,:) = NaN;
%             f00(ci) = NaN;
%             %break;
%         end
%     end
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear tempgenepos tempgeneprobs temprxnpos k
    %disp(bnumstobekoed(ci))
end

% f_ko = -f1';
% v_ko = v_ko;
%
% f00_ko(:,scou) = v00(:,find(weights));
% v00_ko = v00;
%
% f = f00_ko;
% v = v00_ko;
v = v00;

%lostxns(:,scou) = lost_xn;
lostxns = [];

end


function [probtfgene lost_xn datathreshflag] = probtfgenecalc(expression,expressionid,regulated,regulator,datathreshoff,datathreshOE,DATATHRESHVAL,DATATHRESHVALOE,litevidence,prob_prior)
% i need to find the reactions that each gene has influence on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expressionid = expressionid;
tic
lost_xn = false(size(regulated));
% remove_inte.ractions1 = false(size(regulated));
disp('finding probabilities')
% cou = 1;cou1 = 1;cou3 = 1;
%data= knnimpute(datbackup);
data = expression;
data = knnimpute(data);
data = quantilenorm(data);  %its already normalized..
data1 = data;
%kappavals = [0,0.0001,0.001,0.05,0.1,0.25,0.33,0.5,0.75,1];
%datathreshvals = [0,0.01,0.05,0.1,0.2,0.25,0.33,0.4,0.5,0.75,1];
%datathresh = quantile(data(:),datathreshvals(scou));
if isempty(datathreshoff)
     if (~exist('DATATHRESHVAL','var')) | (isempty(DATATHRESHVAL))
        DATATHRESHVAL = 0.33; 
     end
    datathreshoff = quantile(data(:),DATATHRESHVAL);
end
% datathresh = quantile(data(:),DATATHRESHVAL);

if isempty(datathreshOE)
     if (~exist('DATATHRESHVALOE','var')) | (isempty(DATATHRESHVALOE))
        DATATHRESHVALOE = 0.9; 
     end
    datathreshOE = quantile(data(:),DATATHRESHVALOE);
end

if datathreshoff < 0,
    data(data1>=datathreshoff) = 1;
    data(data1 < datathreshoff) = 0;
else
    data(data1 < datathreshoff) = 0;
    data(data1>=datathreshoff) = 1;
end

% overexpression is 2
data(data1>=datathreshOE) = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try probtfgene = ones(size(regulated));
catch M1
    probtfgene = 1; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = cellfun(@(x) find(strcmp(x,expressionid)),regulated,'UniformOutput',false);
l = cellfun(@(x) find(strcmp(x,expressionid)),regulator,'UniformOutput',false);

for  i = 1:length(regulated)
%         k = find(ismember(expressionid,regulated(i)));
%         l = find(ismember(expressionid,regulator(i)));
    if ~isempty(k{i}) & ~isempty(l{i})
        tarexp = data1(k{i},:);
        tarbin = data(k{i},:);
        tfbin = data(l{i},:);

        try ksh = kstest2(tarexp(tfbin == 2),tarexp(tfbin== 0));
            if  ksh

                prob1 = sum(tarbin(tfbin == 2) > 0)/numel(tarbin(tfbin==2));
                probtfgene(i) = prob1;
                %   this formula also gives the same answer  - (sum(~tfbin(tarbin == 1))/numel(tfbin(tarbin==1))) * (sum(tarbin)/sum(~tfbin))

%                 else

%                     probtfgene(i) = 1;  % no effect

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

probtfgene = probtfgene(:);

toc

datathreshflag = 0;

if ~isempty(litevidence)
    probtfgene(litevidence) = prob_prior(litevidence);  % u could set those interactions that u think have strong literature evidence to have predefined
end                                                            % probabilities

toc


end
