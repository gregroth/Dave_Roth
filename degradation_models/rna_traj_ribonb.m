%% approximate the mean degradation dynamics of mRNA for the ribosome number model
%variables: y=[nuclear/cyto,off,on,nb ribo,lifestate]
%
%parameters:
%   - kon = rate at wich an RNA switch on translation   
%   - koff = rate at wich an RNA switch off translation     
%   - delta = magnitude of increment in the rate of degradation
%   - mu = export rate
%   - ini = initiation rate when ON
%
% Inputs:
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% Author: Gregory Roth
%
%   original version: 14.03.2022,
%   last version: 14.03.2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rna=rna_traj_ribonb(param,nb_rna,nb_nuc,tpts)
tmax=15;

%% simulations
nb_sim=15;
for k=1:nb_sim
initialcondition_nuc=[1,0,0,0,0];
    for i=1:nb_nuc
        [t,x]=single_rna_ribonb(param,tmax,initialcondition_nuc);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    for i=nb_nuc+1:nb_rna
        [t,x]=single_rna_ribonb(param,tmax,initialcondition_nuc);
        xprint=print_traj(x,t,tpts+1);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    rnatot(k,:)=size(x5,1)-sum(x5)-sum(x1);
end
rna=mean(rnatot);

end