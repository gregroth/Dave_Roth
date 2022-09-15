function [t,x]=single_rna_trans_state(param,tmax,initialcondition)
%   Gregory Roth, 2022
%   FMI, gregory.roth@fmi.ch

%% Reaction network:
%   1. export :          export    --mu--> export - 1
%   2. translation ON                 OFF     --kon--> ON
%   3. translation OFF                 ON     --koff--> OFF
%   5. degradation ON:       ON  --> deg + ON


%% Rate constants
p.kon = param(1);%0.01;      
p.koff = param(2);%1;  
p.delta = param(4); 
p.export = param(5);

%% Initial state
tspan = [0, tmax]; %seconds
x0    = initialcondition;     %trans on, trans off

%% Specify reaction network
prop = @propensities;
stoich_matrix = [-1 1 0 0  %export
                 0 -1 1 0 %translation ON 
                 0 1 -1 0 %translation OFF 
                 0 0 0 1];    %degradation (proportional to the # of ribo that have passed through the rna)

%% Run simulation
[t,x]=directMethod(stoich_matrix, prop, tspan, x0, p);
end


function a = propensities(x, p)
% Return reaction propensities given current state x
nuclear_rna = x(1);
trans_off = x(2);
trans_on    = x(3);
degraded_rna= x(4);
a=[p.export*nuclear_rna.*(1-degraded_rna);       %export
     p.kon*trans_off.*(1-degraded_rna).*(1-nuclear_rna);       %translation ON 
     p.koff*trans_on.*(1-degraded_rna).*(1-nuclear_rna);%translation OFF
     p.delta.*trans_on.*(1-degraded_rna).*(1-nuclear_rna)]; %degradation 
end
