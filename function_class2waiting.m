function [f] = function_class2waiting(level,phase,c)
    %this function is prob of being in state (a,b)
    serverAvailable = max(c-phase,0);
    if level >= serverAvailable % This function was wrong before because I had > instead of >=
        f = 1;
    else
        f=0;
    end
end

