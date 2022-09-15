%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit the nmd model
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
addpath('/../models/stochastic')
addpath('../gillespie/')

%data
expnamelist={'pPD17'};
conditionlist={'fe'};

%bestfit
load(strcat('./bestfits/bestfit_pPD12_fe.mat'),'bestfit_parameter')


%%
cl=clock;

%% selct the fit strategy

nbrun=50; %number of run for the search
%% load the data
load(strcat('../data/data_for_fit/data_',expnamelist{1},'_',conditionlist{1}));
datatofit=dataall(2,:);
nb_rna_ini=round(sum(dataall(1:2,1)));
nb_rna=round(max(sum(dataall([1,2],:))));
nb_nuc=round(dataall(1,1))+(nb_rna-nb_rna_ini);


sigma=1;
N=length(timepts);
MLF_tot=@(x)N*log(sigma*sqrt(2*pi))+sum((datatofit-rna_traj_nmd_noHistory([bestfit_parameter,x(1),.2],nb_rna,nb_nuc,timepts)).^2)/(2*sigma^2);     
%parameter limits
lb=0;%     x =[deltaon,deltaoff]
ub=.15;
%initial condition for the parameters
P0=rand(1)*.15;

%% search
ms = MultiStart('UseParallel',false,'Display','iter');
problem = createOptimProblem('fmincon',...
     'objective',@(x)MLF_tot(x),...
     'x0',P0,'lb',lb,'ub',ub);
[x,fval,eflag,output,manymins] = run(ms,problem,nbrun);

%% calculate the bestfit parameters
bestfit_parameter=[bestfit_parameter,x(1),.2];

%save
filename=strcat('./bestfits/bestfit_nmd_',expnamelist{1},'_',erase(num2str(cl([1,2,3,4,5]))," "),'_',conditionlist{1});
save(filename,'x','fval','eflag','output','manymins','nbrun','condition','timepts','bestfit_parameter','datatofit','expname','dataall');



%% remove paths
rmpath('../models/stochastic/')
rmpath('../gillespie/')



