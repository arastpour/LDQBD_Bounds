function [ p ] = function_ErlangA_truncation_level( lambda, mu, gamma, c, epsilon_h )
%ErlangA_truncation_level: Computes a truncation level p that guarantees
% that the probability in the truncated upper tail is less than epsilon_h, for 
% an Erlang A system with parameters:
% lambda: arrival rate
% mu: service rate
% gamma: abandonment rate
% c: number of servers

%Initialization: 
p=1; l=0; u=0; q_l=1; q_u=1; converge=0; Delta2=1;

while ~converge,
    if (q_u > q_l)||(l==c),
        u=u+1;
        b=lambda/(c*mu+u*gamma);
        q_u=b*q_u;
        % Normalize
        Sum=1+q_u;
        q_u=q_u/Sum;q_l=q_l/Sum;
        if b<1, Delta2=b*q_u/(1-b); end
    else
        l=l+1;
        a=(c-l+1)*mu/lambda;
        q_l=a*q_l;
        % Normalize;
        Sum=1+q_l;
        q_u=q_u/Sum;q_l=q_l/Sum; 
    end
    if  Delta2 <epsilon_h, converge=1; end
end
p=c+u;

end

