%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit the initiation-mediated model with constant probability
%
%
% Inputs:
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% Author: Gregory Roth
%
%   original version: 29.04.2022,
%   last version: 29.04.2022%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../models/')
addpath('../gillespie/')
%data
expnamelist={'pPD12'};
conditionlist={'fe'};

%%
cl=clock;
rng('shuffle');
%% selct the fit strategy

nbrun=1000; %number of run for the search
%% load the data
load(strcat('../data/data_for_fit/data_',expnamelist{1},'_',conditionlist{1}));
datatofit=dataall(2,:);
nb_rna_ini=round(sum(dataall(:,1)));
nb_rna=round(max(sum(dataall)));
nb_nuc=round(dataall(1,1))+(nb_rna-nb_rna_ini);
 
sigma=1;
N=length(timepts);
MLF_tot=@(x)N*log(sigma*sqrt(2*pi))+sum((datatofit-rna_traj_probaribo_flux([kon,koff,x(1),export,ini],nb_rna,nb_nuc,timepts)).^2)/(2*sigma^2);    
%parameter limits
lb=1e-05;%     x =[deltaon,deltaoff]
ub=9e-05;
%initial condition for the parameters
P0=4.0883e-05;
%% search
ms = MultiStart('UseParallel',false,'Display','iter');
problem = createOptimProblem('fmincon',...
     'objective',@(x)MLF_tot(x),...
     'x0',P0,'lb',lb,'ub',ub);
[x,fval,eflag,output,manymins] = run(ms,problem,nbrun);

%% calculate the bestfit parameters
bestfit_parameter=[kon,koff,x(1),export,ini];

save('bestfit_param_temp','bestfit_parameter');
%save
filename=strcat('./bestfits/bestfit_probaribo_flux_',expnamelist{1},'_',erase(num2str(cl([1,2,3,4,5]))," "),'_',conditionlist{1});
save(filename,'x','fval','eflag','output','manymins','nbrun','condition','datatofit','timepts','bestfit_parameter','expname','dataall');

%% remove paths
rmpath('../models/')
rmpath('../gillespie/')



