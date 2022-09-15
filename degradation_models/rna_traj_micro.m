%% approximate the mean degradation dynamics of mRNA for the micro-flux model
%variables: y=[nuclear/cyto,off,on,nb ribo,lifestate]
%
%parameters:
%   - kon = rate at wich an RNA switch on translation   
%   - koff = rate at wich an RNA switch off translation     
%   - delta = magnitude of increment in the rate of degradation
%   - delta_micro_nontrans = miRNA mediated degradation rate when OFF
%   - delta_micro_trans = miRNA mediated degradation rate when ON
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
function rna=rna_traj_micro(param,nb_rna,nb_nuc,tpts)
%% parameters
kon = param(1);%0.01;      
koff = param(2);%1;                     
delta =param(3);
export = param(4);
ini = param(5);
delta_micro=param(6);
delta_micro_trans=param(7);


tmax=15;


%% simulations
nb_sim=20;
for k=1:nb_sim
initialcondition_cyto=[0,1,0,0,0];
initialcondition_nuc=[1,0,0,0,0];
    for i=1:nb_nuc
        [t,x]=single_rna_wflux_micro_trans(kon,koff,delta,export,ini,delta_micro,delta_micro_trans,tmax,initialcondition_nuc);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    for i=nb_nuc+1:nb_rna
        [t,x]=single_rna_wflux_micro_trans(kon,koff,delta,export,ini,delta_micro,delta_micro_trans,tmax,initialcondition_cyto);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    rnatot(k,:)=size(x5,1)-sum(x5)-sum(x1);
end
rna=mean(rnatot);

end