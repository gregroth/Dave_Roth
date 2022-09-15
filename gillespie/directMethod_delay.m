function [ t, x ] = directMethod_delay(stoich_matrix, propensity_fcn, tspan, x0,...
                                  delay_reaction,stoich_delay,delay_time,params)
%% Attention this is a modified version of the Gillespie algorithm to take into account reaction with delay
%   Gregory Roth, 2022
%   FMI, gregory.roth@fmi.ch
%%
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
        %check if there is a delay reaction happening in the time interval
        %[T(rxn_count),T(rxn_count)   + tau]

        [valmin,kmin]=min(termination_times);
        if valmin<=T(rxn_count)+tau %if there is one, we execute it and forget about mu
            T(rxn_count+1)   = valmin;
            termination_times(kmin)=[]; %remove the delayed reaction that happened
            X(rxn_count+1,:) = X(rxn_count,:) + stoich_delay;    
            rxn_count = rxn_count + 1;
        elseif mu==delay_reaction %check if mu is the delayed reaction
            termination_times(end+1)=T(rxn_count)+tau+delay_time;
            T(rxn_count+1)   = T(rxn_count)   + tau;
            X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);    
            rxn_count = rxn_count + 1;
        else
            T(rxn_count+1)   = T(rxn_count)   + tau;
            X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);    
            rxn_count = rxn_count + 1;
        end


    end  
    %return simulation trajectory
    t = T(1:reac_count);
    x = X(1:reac_count,:);
    if t(end) > tspan(2)
        t(end) = tspan(2);
        x(end,:) = X(reac_count-1,:);
    end    
end



