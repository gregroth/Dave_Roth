function [ t, x ] = directMethod_probaribo_flux( stoich_matrix, propensity_fcn, tspan, x0,...
                                  params,ini_reaction,prob_delta)
%% ATTENTION this is a modified version of the directMethod function specially written for the initiation mediated degradation pathway
% When a ribo passes through, the degradation pathway can be activated with a 
% probability p_delta that depends on the number of initiation events
%%
    
    output_length = 1000000;
    num_variable = size(stoich_matrix, 2);
    T = zeros(output_length, 1);
    X = zeros(output_length, num_variable);
    T(1)     = tspan(1);
    X(1,:)   = x0;
    reac_count = 1;
    %degradation following nmd pathway
    stoich_nmnd=zeros(1,num_species);
    stoich_nmnd(end)=1;

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

        % Update time and carry out reaction mu

        if mu==ini_reaction %check if a ribosome has been loaded
            r_delta=rand(1);
            if r_delta<=1-exp(-delta*X(rxn_count,4)) %degradation is initiated with proba fct of nb of initiation events
                T(rxn_count+1)   = T(rxn_count)   + tau;
                X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:) +stoich_nmnd;    
                rxn_count = rxn_count + 1;
            else
                T(rxn_count+1)   = T(rxn_count)   + tau;
                X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);    
                rxn_count = rxn_count + 1;
            end
        else
            T(rxn_count+1)   = T(rxn_count)   + tau;
            X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);    
            rxn_count = rxn_count + 1;
        end
    end  

    %return simulation trajectory
    t=T(1:rxn_count);
    x=X(1:rxn_count,:);
    if t(end)>tspan(2)
        t(end)=tspan(2);
        x(end,:)=X(rxn_count-1,:);
    end    

end

