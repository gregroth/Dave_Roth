%% approximate the mean degradation dynamics of mRNA for translation state dependent degradation model
%variables: y=[nuclear/cyto,off,on,lifestate]
%
%parameters:
%   - kon = rate at wich an RNA switch on translation   
%   - koff = rate at wich an RNA switch off translation     
%   - delta = degradation rate when ON
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
function rna=rna_traj_transstate(param,nb_rna,nb_nuc,tpts)
tmax=15;
%% simulations
nb_sim=15;
for k=1:nb_sim
initialcondition_nuc=[1,0,0,0];
initialcondition_rna=[0,1,0,0];
    for i=1:nb_nuc
        [t,x]=single_rna_trans_state(param,tmax,initialcondition_nuc);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x4(i,:)=xprint(:,4)';
    end
    for i=nb_nuc+1:nb_rna
        [t,x]=single_rna_trans_state(param,tmax,initialcondition_rna);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x4(i,:)=xprint(:,4)';
    end
    rnatot(k,:)=size(x4,1)-sum(x4)-sum(x1);
end
rna=mean(rnatot);

end