%% This function calculates the trajectories at specific time points tp
%   Gregory Roth, 2022
%   FMI, gregory.roth@fmi.ch
%%

function xp=print_traj(x,t,tp)
    xp=zeros(length(tp),size(x,2));
    for i=1:length(tp)
        currentk=find(t<=tp(i),1,'last');
        xp(i,:)=x(currentk,:);
    end
end
