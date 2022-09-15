%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit the ribosome flux model
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
addpath('../models/stochastic/')
addpath('../gillespie/')
%data
expnamelist={'pPD12'};
conditionlist={'ctrl'};

%%
cl=clock;

%% selct the fit strategy

nbrun=20; %number of run for the search
%% load the data
load(strcat('../../data/data_for_fit/data_',expnamelist{1},'_',conditionlist{1}));
datatofit=dataall(2,:);
nb_rna_ini=round(sum(dataall(:,1)));
nb_rna=round(max(sum(dataall)));
nb_nuc=round(dataall(1,1))+(nb_rna-nb_rna_ini);
    
sigma=1;
N=length(timepts);
MLF_tot=@(x)N*log(sigma*sqrt(2*pi))+sum((datatofit-rna_traj_wflux([kon,koff,x(1),export,ini],nb_rna,nb_nuc,timepts)).^2)/(2*sigma^2);     
%parameter limits
lb=0.001;%     x =[deltaon,deltaoff]
ub=.006;
%initial condition for the parameters
P0=.0036;


%% search
ms = MultiStart('UseParallel',false,'Display','iter');
problem = createOptimProblem('fmincon',...
     'objective',@(x)MLF_tot(x),...
     'x0',P0,'lb',lb,'ub',ub); 
[x,fval,eflag,output,manymins] = run(ms,problem,nbrun);

%% calculate the bestfit parameters 
bestfit_parameter=[kon,koff,x(1),export,ini];

%save
filename=strcat('./bestfits/bestfit_',expnamelist{1},'_',erase(num2str(cl([1,2,3,4,5]))," "),'_',conditionlist{1});
save(filename,'x','fval','eflag','output','manymins','nbrun','condition','datatofit','timepts','bestfit_parameter','expname','dataall');


%% remove paths
rmpath('../models/stochastic')
rmpath('../gillespie/')



