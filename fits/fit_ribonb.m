%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit the ribosome number model
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
addpath('../delay_telegraph_model/')

%data
expnamelist={'pPD12'};
conditionlist={'fe'};
%%
cl=clock;
rr=rand(3);
nbrun=50; %number of run for the search
%% load the data
load(strcat('../../data/data_for_fit/data_',expnamelist{1},'_',conditionlist{1}));
datatofit=dataall(2,:);
nb_rna_ini=round(sum(dataall([1,2],1)));
nb_rna=round(max(sum(dataall([1,2],:))));
nb_nuc=round(dataall(1,1))+(nb_rna-nb_rna_ini);

%ATTENTION the ORF is different than from the suntag experiment
t_elongation = 77.5./3600;
    
%calculating the steady-state ribo distribution
load(strcat('../../delay_telegraph_model/data_for_fit/dataToFit.mat'),'dataVal_fe','maxRna_ire')
dataVal=dataVal_fe;
distfit_fe= FSP_delaytelegraph([kon,koff,ini,t_elongation],maxRna_ire);
%match the length of the distribution with the one of the data
distfit_fe(length(dataVal)+1:end)=[];
ribo_nbs=0:1:length(distfit_fe)-1;

sigma=1;
N=length(timepts);

MLF_tot=@(x)N*log(sigma*sqrt(2*pi))+sum((datatofit-rna_traj_ribonb_noHistory([kon,koff,x(1),export,ini,t_elongation],nb_rna,nb_nuc,distfit_fe,ribo_nbs,timepts)).^2)/(2*sigma^2);

%parameter limits
lb=0;%     x =[deltaon,deltaoff]
ub=.4;
%initial condition for the parameters
P0=rand(1)*.05;

%% search
ms = MultiStart('UseParallel',false,'Display','iter');
problem = createOptimProblem('fmincon',...
     'objective',@(x)MLF_tot(x),...
     'x0',P0,'lb',lb,'ub',ub);

[x,fval,eflag,output,manymins] = run(ms,problem,nbrun);

%% calculate the bestfit parameters
bestfit_parameter=[kon,koff,x(1),export,ini,t_elongation];
%save
filename=strcat('bestfits/bestfit_ribonb_noHistory_',expnamelist{1},'_',erase(num2str(cl([1,2,3,4,5]))," "),'_',conditionlist{1});
save(filename,'x','fval','eflag','output','manymins','nbrun','condition','datatofit','timepts','bestfit_parameter','expname','dataall');



%% remove paths
rmpath('../models/')
rmpath('../gillespie/')
rmpath('../delay_telegraph_model/')

