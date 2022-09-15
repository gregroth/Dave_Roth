%% routine for fitting the delayed two-state model to data 
% Inputs:
            % trmean= vector with mean transcriptional levels
            % cpmean= vector with EP contact probability corresponding to trmean
            % dataVal=cell with the histogram values for each distribution
            % binsize= vector with the size of the histogram bins
            % maxRna= vector with the max number of mRNA per cell considered in the histograms
            % cpdist= vector with EP contact probability corresponding to distribution data
            % selectedclones = ector of 0 and 1 of size= # clone with distribution; 1= we use this clone for the fit
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none
%
% Author: Gregory Roth
%
%   original version: 19.08.2022,
%   last version: 19.08.2022%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% Shuffle the seed
cl=clock;
rng('shuffle')
%% data used for the fitting
    % distribution of mRNA data
    dataVal=[];
    binsize=[];
    maxRna=[];

%% load the data
load('./data_for_fit/dataToFit.mat')
load('./data_for_fit/histograms.mat')
load('./data_for_fit/histogram_corrected.mat')
%% load the bestfit parameter

%% data for specific condition
condition=2 ;
if condition==1
    dataVal=dataVal_ctrl;
    maxRna=maxRna_ire;
    N=sum(histogram_ctrl(:,2));
elseif condition==2
    dataVal=dataVal_fe;
    maxRna=maxRna_ire;
    N=sum(histogram_fe(:,2));
elseif condition==3
    dataVal=dataVal_nonstemloop;
    maxRna=maxRna_nonstemloop;
    N=sum(histogram_nonstemloop(:,2));
elseif condition==4
    dataVal=dataVal_stemloop;
    maxRna=maxRna_stemloop;
    N=sum(histogram_stemloop(:,2));
elseif condition==5
    dataVal=dataVal_xx;
    maxRna=maxRna_xx;
    N=sum(histogram_xx(:,2));
end
%% Fitting parameters
%number of runs
nbrun=2000;
%free parameters (kon,koff,mu)
freeparam=[1,1,1];
%fixed parameters
%fixedparam=[1.5279    1.5868];
%fixedparam=[0.1022    0.0491    4.8581];
fixedparam=[];
%selected clones for the fit, belongs to 1,2,...,7
%selectedclones=[]; %must be a vector of o and 1 of of size= # clone with distribution
%constrain bounds 
lb = [0,0,0];
ub = [10,10,50];
%initial point
x0temp=[0.0427    3.6310   10.0483];
x0temp=rand(1,3);

initialpt=x0temp(freeparam==1);

%% multistart fit
[x,fval,eflag,output,manymins]=fit_fct(dataVal,binsize,maxRna,nbrun,freeparam,fixedparam,lb,ub,initialpt);
fval=N*fval;
%save the bestfit data
Xtemp=arrayfun(@(x)x.X,manymins,'UniformOutput',false);
Xlist=cell2mat(Xtemp);
Xmat=reshape(Xlist,sum(freeparam),length(manymins))';
freeparamstr=num2str(freeparam);
filesavename=strcat('test_bestfit_condition_',num2str(condition),'_',freeparamstr(~isspace(freeparamstr)),'_',num2str(ceil(100*rand)),'_',erase(num2str(cl([1,2,3,4,5]))," "));
bestfit_param=zeros(1,3);
bestfit_param(freeparam==1)=x;
bestfit_param(freeparam==0)=fixedparam;
save(strcat('./bestfit_parameters/',filesavename),'x','fval','eflag','output','manymins','Xmat','bestfit_param','freeparam','fixedparam','condition')
%save the the same information with current name
save(strcat('./bestfit_parameters/current_bestfit.mat'),'x','fval','eflag','output','manymins','Xmat','bestfit_param','freeparam','fixedparam','condition')

toc

