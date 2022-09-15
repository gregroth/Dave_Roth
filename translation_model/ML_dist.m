%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Gregory Roth
%
%   original version: 29.04.2022,
%   last version: 29.04.2022%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=ML_dist(freeparam,fixedparam,x,maxRna,dataVal,binsize)
    totparam=zeros(1,3);
    totparam(freeparam==1)=x;
    totparam(freeparam==0)=fixedparam;
    kon=totparam(1);
    koff=totparam(2);
    mu=totparam(3);
    out=ML_delayedtelegraph_hist([kon,koff,mu],maxRna,dataVal,binsize);
end
