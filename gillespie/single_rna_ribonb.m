function [t,x]=single_rna_ribonb(parameters,tmax,initialcondition)
%   Gregory Roth, 2022
%   FMI, gregory.roth@fmi.ch

%% Reaction network:
%   1. export :          export    --mu--> export - 1
%   2. translation ON                 OFF     --kon--> ON
%   3. translation OFF                 ON     --koff--> OFF
%   4. initiation event     ribo --> ribo +1 triggers termination with
%   delay
%   5. degradation ON:       live --> deg



%% Rate constants
p.kon = parameters(1);     
p.koff = parameters(2);                     
p.delta = parameters(3); 
p.export = parameters(4);
p.ini = parameters(5);
t_elongation=parameters(6);
%% Initial state
tspan=[0,tmax]; 
x0=initialcondition;

%% Specify reaction network
prop = @propensities;
stoich_matrix = [-1 1 0 0 0     %export
                 0 -1 1 0 0     %translation ON 
                 0 1 -1 0 0     %translation OFF 
                 0 0 0 1 0      %initiation event
                 0 0 0 0 1];    %degradation (proportional to the # of ribo that have passed through the rna)

%% Run simulation
stoich_delay=[0,0,0,-1,0];
[t,x] = directMethod_delay(stoich_matrix, prop, tspan, x0,4,stoich_delay,t_elongation, p);
end


function a = propensities(x, p)
% Return reaction propensities given current state x
nuclear_rna = x(1);
trans_off = x(2);
trans_on    = x(3);
ribo_nbr = x(4);
degraded_rna= x(5);
a=[p.export*nuclear_rna.*(1-degraded_rna);       %export
     p.kon*trans_off.*(1-degraded_rna).*(1-nuclear_rna);       %translation ON 
     p.koff*trans_on.*(1-degraded_rna).*(1-nuclear_rna);       %translation OFF 
     p.ini.*trans_on.*(1-degraded_rna).*(1-nuclear_rna);       %initiation event
     p.delta.*max(0,ribo_nbr).*(1-degraded_rna).*(1-nuclear_rna)]; %degradation 
end
