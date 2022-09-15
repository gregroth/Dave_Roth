%% approximate the mean degradation dynamics of mRNA for the NMD-flux model
%variables: y=[nuclear/cyto,off,on,nb ribo,lifestate]
%
%parameters:
%   - kon = rate at wich an RNA switch on translation   
%   - koff = rate at wich an RNA switch off translation     
%   - delta = magnitude of increment in the rate of degradation
%   - prob_nmd = probability to initiate NMD after an initiation event
%   - resistant = proportion of mRNA that are resistant to NMD
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
function rna=rna_traj_nmd(param,nb_rna,nb_nuc,tpts)
%% parameters
kon = param(1);%0.01;      
koff = param(2);%1;                     
delta =param(3);
export = param(4);
ini = param(5);
prob_nmd=param(6);
resistent_frac=param(7);
tmax=15;

nb_cyto=nb_rna-nb_nuc;
nb_rna_resistent_nuc=round(resistent_frac*nb_nuc);
nb_rna_resistent_cyto=round(resistent_frac*nb_cyto);

      
initialcondition_nuc=[1,0,0,0,0];
initialcondition_cyto=[0,1,0,0,0];
nb_sim=15;
for k=1:nb_sim
    for i=1:nb_rna_resistent_nuc
        [t,x]=single_rna_wflux(kon,koff,delta,export,ini,tmax,initialcondition_nuc);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    for i=nb_rna_resistent_nuc+1:nb_nuc
        [t,x]=single_rna_wflux_nmd(kon,koff,delta,export,ini,prob_nmd,tmax,initialcondition_nuc);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    for i=nb_nuc+1:nb_nuc+nb_rna_resistent_cyto
        [t,x]=single_rna_wflux(kon,koff,delta,export,ini,tmax,initialcondition_cyto);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    for i=nb_nuc+nb_rna_resistent_cyto+1:nb_rna
        [t,x]=single_rna_wflux_nmd(kon,koff,delta,export,ini,prob_nmd,tmax,initialcondition_cyto);
        xprint=print_traj(x,t,tpts);
        x1(i,:)=xprint(:,1)';
        x5(i,:)=xprint(:,5)';
    end
    
    rnatot(k,:)=size(x5,1)-sum(x5)-sum(x1);
end
rna=mean(rnatot);

end