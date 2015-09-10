

% Adjust Transport Reactions:
iSM810noExchange = changeRxnBounds(iSM810,iSM810exchangeRxns(:,1),0,'u');

%% Griffin
iSM810G = iSM810noExchange;
iSM810G = changeRxnBounds(iSM810G,'R800',1,'u'); % ammonia
iSM810G = changeRxnBounds(iSM810G,'R804',1,'u'); % oxygen
iSM810G = changeRxnBounds(iSM810G,'R805',1000,'u'); %CO2
iSM810G = changeRxnBounds(iSM810G,'R812',1,'u'); % glycerol
iSM810G = changeRxnBounds(iSM810G,'R822',1,'u'); % Asparagine
iSM810G = changeRxnBounds(iSM810G,'R838',1,'u'); % Phosphate
iSM810G = changeRxnBounds(iSM810G,'R882',1,'u'); % Phosphate
iSM810G = changeRxnBounds(iSM810G,'R841',1,'u'); % Sulfate
iSM810G = changeRxnBounds(iSM810G,'R851',1,'u'); % Citrate
iSM810G = changeRxnBounds(iSM810G,'R904',1,'u'); % Citrate
iSM810G = changeRxnBounds(iSM810G,'R858',1,'u'); % Ethanol
iSM810G = changeRxnBounds(iSM810G,'R924',1,'u'); % Iron (FE3)

%% Sauton
iSM810S = iSM810G;
iSM810S = changeRxnBounds(iSM810S,'R858',0,'u'); % Ethanol

% carbon source models (based on Sauton for non-carbon metabolites)
tmp = changeRxnBounds(iSM810S,'R812',0,'u');
tmp = changeRxnBounds(tmp,'R851',0,'u');
tmp = changeRxnBounds(tmp,'R904',0,'u');
% tmp = changeRxnBounds(tmp,'R924',10^-4,'u');

iSM810Sac = changeRxnBounds(tmp,'R808',1,'u'); % acetate
iSM810Sac = changeRxnBounds(iSM810Sac,'R844',1,'u'); % acetate

iSM810Sakg = changeRxnBounds(tmp,'R847',1,'u'); % 2-oxoglutarate
iSM810Sala = changeRxnBounds(tmp,'R820',1,'u'); % L-alanine
iSM810Sasn = changeRxnBounds(tmp,'R822',1,'u'); % Asparagine

iSM810Scaproate = changeRxnBounds(tmp,'R908',1,'u'); % caproic acid (hexanoate)

iSM810Schol = changeRxnBounds(tmp,'R932',1,'u'); % cholesterol

iSM810Scit = changeRxnBounds(tmp,'R851',1,'u'); % citrate
iSM810Scit = changeRxnBounds(iSM810Scit,'R904',1,'u');

iSM810Sdala = changeRxnBounds(tmp,'R853',1,'u'); % D-alanine

iSM810Sf6p = changeRxnBounds(tmp,'R899',1,'u'); % D-fructose-6-phosphate

iSM810Sg6p = changeRxnBounds(tmp,'bDG6P_Transport',1,'u'); % D-glucose-6-phosphate

iSM810Sglc = changeRxnBounds(tmp,'R863',1,'u'); % D-glucose

iSM810Sgln = changeRxnBounds(tmp,'R829',1,'u'); % l-glutamine

iSM810Sglu = changeRxnBounds(tmp,'R830',1,'u'); % l-glutamate
iSM810Sglu = changeRxnBounds(iSM810Sglu,'R864',1,'u'); % l-glutamate

iSM810Sgly = changeRxnBounds(tmp,'R866',1,'u'); % glycine

iSM810Sllac = changeRxnBounds(tmp,'R876',1,'u'); % l-lactate

iSM810Smal = changeRxnBounds(tmp,'R879',1,'u'); % l-malic acid
iSM810Smal = changeRxnBounds(iSM810Smal,'R880',1,'u'); % l-malic acid

iSM810Soleate = changeRxnBounds(tmp,'R918',1,'u'); % oleate (9-octadecenoate)

iSM810Spalmitate = changeRxnBounds(tmp,'R916',1,'u'); % palmitate (hexadecanoate)

iSM810Spropionate = changeRxnBounds(tmp,'R907',1,'u'); % propionate

iSM810Spyruvate = changeRxnBounds(tmp,'R885',1,'u'); % pyruvate

iSM810Stre = changeRxnBounds(tmp,'R902',1,'u'); % trehalose

% nitrogen sources
tmp = changeRxnBounds(iSM810S,'R800',0,'u');
tmp = changeRxnBounds(tmp,'R822',0,'u');

iSM810SnitALA = changeRxnBounds(tmp,'R820',1,'u'); % Alanine

iSM810SnitASN = changeRxnBounds(tmp,'R822',1,'u'); % Asparagine

iSM810SnitGLN = changeRxnBounds(tmp,'R829',1,'u'); % Glutamine

iSM810SnitGLU = changeRxnBounds(tmp,'R830',1,'u'); % l-glutamate
iSM810SnitGLU = changeRxnBounds(iSM810SnitGLU,'R864',1,'u'); % l-glutamate

iSM810SnitSER = changeRxnBounds(tmp,'R886',1,'u'); % serine



%% 7H9
iSM810_7H9 = iSM810noExchange;
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R800',1,'u'); % ammonia
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R804',1,'u'); % oxygen
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R805',1000,'u'); %CO2
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R812',1,'u'); % glycerol
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R830',1,'u'); % glutamate
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R864',1,'u'); % glutamate
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R838',1,'u'); % Phosphate
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R882',1,'u'); % Phosphate
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R841',1,'u'); % Sulfate
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R851',1,'u'); % Citrate
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R904',1,'u'); % Citrate
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R858',1,'u'); % Ethanol
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R863',1,'u'); % Glucose
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R924',1,'u'); % Iron (FE3)
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R925',1,'u'); % Biotin

% to account for BSA:
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R820',10^-4,'u'); % ALA
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R821',10^-4,'u'); % ARG
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R822',10^-4,'u'); % ASN
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R823',10^-4,'u'); % ASP
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R848',10^-4,'u'); % ASP
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R849',10^-4,'u'); % ASP
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R850',10^-4,'u'); % ASP
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R826',10^-4,'u'); % CYS
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R829',10^-4,'u'); % GLN
% iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R830',10^-4,'u'); % GLU
% iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R864',10^-4,'u'); % GLU
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R866',10^-4,'u'); % GLY
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R831',10^-4,'u'); % HIS
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R868',10^-4,'u'); % HIS
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R869',10^-4,'u'); % HIS

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R832',10^-4,'u'); % ILE
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R870',10^-4,'u'); % ILE
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R871',10^-4,'u'); % ILE

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R833',10^-4,'u'); % LEU
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R874',10^-4,'u'); % LEU
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R875',10^-4,'u'); % LEU

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R834',10^-4,'u'); % LYS
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R877',10^-4,'u'); % LYS
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R878',10^-4,'u'); % LYS

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R835',10^-4,'u'); % MET
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R881',10^-4,'u'); % PHE

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R839',10^-4,'u'); % PRO
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R883',10^-4,'u'); % PRO
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R884',10^-4,'u'); % PRO

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R886',10^-4,'u'); % SER
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R887',10^-4,'u'); % SER

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R842',10^-4,'u'); % THR
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R891',10^-4,'u'); % THR
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R892',10^-4,'u'); % THR

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R893',10^-4,'u'); % TRP
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R894',10^-4,'u'); % TYR

iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R843',10^-4,'u'); % VAL
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R897',10^-4,'u'); % VAL
iSM810_7H9 = changeRxnBounds(iSM810_7H9,'R898',10^-4,'u'); % VAL
