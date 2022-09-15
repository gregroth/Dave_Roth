function [ t, x ] = directMethod( stoich_matrix, propensity_fcn, tspan, x0,...
                                  params)
%   2022-09-15. This code is adapted from the functon direcMethod from Nezar Abdennur,
%   Created: 2012-01-19, www.sysbiolab.uottawa.ca
                            
                             
output_length = 1000000;
num_variable = size(stoich_matrix, 2);
T = zeros(output_length, 1);
X = zeros(output_length, num_variable);
T(1)     = tspan(1);
X(1,:)   = x0;
reac_count = 1;

while T(reac_count) < tspan(2)        
    %reaction propensities
    a = propensity_fcn(X(reac_count,:), params);
    
    % time at next reaction (tau)
    a0 = sum(a);
    r = rand(1,2);
    tau = -log(r(1))/a0;
    % next reaction
    mu = find((cumsum(a) >= r(2)*a0), 1,'first');
    
  
    if reac_count + 1 > output_length
        t = T(1:reac_count);
        x = X(1:reac_count,:);
        warning('output length is too large');
        return;
    end
    
    % Update time and make reaction mu
    T(reac_count+1)=T(reac_count)+tau;
    X(reac_count+1,:)=X(reac_count,:)+stoich_matrix(mu,:);    
    reac_count = reac_count + 1;
end  
%return simulation trajectory
t = T(1:reac_count);
x = X(1:reac_count,:);
if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,:) = X(reac_count-1,:);
end    
end

