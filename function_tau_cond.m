function tau_condition = function_tau_cond (k,maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c)
    Q_k_1 = function_Q(1, k, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); %Q_1^(k)
    
    Q_k_plus1_1 = function_Q(1, k+1, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); %Q_1^(k+1)
    Q_k_plus1_2 = function_Q(2, k+1, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); %Q_2^(k+1)  
    
    Q_k_plus2_1 = function_Q(1, k+2, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); %Q_1^(k+2)
    Q_k_plus2_2 = function_Q(2, k+2, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); %Q_2^(k+2)
    
    tau_condition = true(1);
    for a = 1:maxPhase+1
        for b = 1:maxPhase+1
            med_cond = ((Q_k_plus1_2(b,b)/(Q_k_1(a,a) * Q_k_plus1_1(b,b))) >= (Q_k_plus2_2(b,b)/(Q_k_plus1_1(a,a) * Q_k_plus2_1(b,b))));
            tau_condition = tau_condition && med_cond;
        end
    end
end

