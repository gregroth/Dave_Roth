%% This function calculates an approximation of the probability distribution
%of the delayed 2-state model. The algorithm follows the stationary Finite State Projectio
%approximation method provided by Ankit Gupta, Jan Mikelson and Mustafa and
%adapted for delay equation by Fu,Tineke,Grima 2022

%Input
%
%           - parm=[kon,koff,kini,delta=1], parameter vector
%           - maxRna = maximum number of RNA molecules
%           - state space = (0,0),...,(0,maxRna+1),(1,0),....,(1,maxRna+1)
%           and 0=OFF and 1=ON
%Output
%
%           -  dist.join= vector containing the distribution of the state
%              of the promoter and the number of RNA
%           -  dist.rna= vector containing the marginal distribution of the
%              number of RNA
%           -  dist.mean= mean of the # RNA
%           -  dist.var= variance of the # RNA

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

function dist = FSP_delaytelegraph(parameter,maxRna) 
    kon=parameter(1);
    koff=parameter(2);
    mu=parameter(3);
    tau=parameter(4);
    
    B=mu.*(diag(-[ones(1,maxRna),0],0)+diag(ones(1,maxRna),-1));
    A1=kon.*eye(maxRna+1);
    A2=koff.*eye(maxRna+1);
   
    A=sparse([-A1,A2;A1,-A2+B]);
    
    %calculate the history
    p0=[zeros(1,maxRna+1),1,zeros(1,maxRna)]';
    tildeP=expm(A.*tau)*p0;
    tildePdiff=tildeP-[0,tildeP(1:end-1)']';
    
    %calculate the delay terms
    h1=mu*kon/(kon+koff);
    
    %calculate the non homogenous term
    nh=h1*tildePdiff;
    
    %calculate the solution
    A(end,:)=[];
    A(:,end)=[];
    nh(end)=[];
    sol = A \ (-nh);
    
    %remove infinite value
    sol(isinf(sol))=0;
    
    %calculate the ditribution of PolII
    dist=sol(1:maxRna+1)+[sol(maxRna+2:end)',0]';
    dist(dist<0)=0;

end