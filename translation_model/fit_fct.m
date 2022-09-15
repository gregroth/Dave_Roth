%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit the delay-two-state model to the ribosome distribution
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

function [x,fval,eflag,output,manymins]=fit_fct(dataVal,binsize,maxRna,nbrun,freeparam,fixedparam,lb,ub,initialpt)
        %constrain bounds 
        lb = lb(freeparam==1);
        ub = ub(freeparam==1);

        %likelihood functions
        MLF = @(x)ML_dist(freeparam,fixedparam,x,maxRna,dataVal,binsize);


        problem = createOptimProblem('fmincon',...
             'objective',@(x)MLF(x),...
             'x0',initialpt,'lb',lb,'ub',ub);

        %ms = MultiStart('PlotFcns',@gsplotbestf);
        ms = MultiStart('UseParallel',false,'Display','iter');
        %set up parallel processing
        %parpool
         % euler = parcluster('local');
          %pool = parpool(euler,24);
        [x,fval,eflag,output,manymins] = run(ms,problem,nbrun);
end
