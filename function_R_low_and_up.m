function [R_lowerbar, R_upperbar, m] = function_R_low_and_up(epsilon_r,ell, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, m_initial)
    Msize = maxPhase+1;
    continue_r = true(1);
    m = m_initial;
    while continue_r
        W1 = zeros((m-1)*Msize, 2*Msize); 
        W1(1:Msize, 1:Msize)                           = function_Q(2, ell+1,   maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); 
        W1((m-2)*Msize+1:(m-1)*Msize, Msize+1:2*Msize) = function_Q(0, ell+m-1, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); 
        V = -function_BlockTridiagonalInverse2(ell, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, m); % V calcualtes the rows 1:Msize of inv(W2)
        T = V(1:Msize, 1:Msize); %Tij = expected time in transient state (ell+1,j) given the syetem starts from (ell+1,i), before visiting levels ell and ell+m 
        Pmat = V * sparse(W1);         
        P1 = Pmat(1:Msize, 1:Msize); %P1ij = prob of getting absorbed to (ell,j) given the system starts from (ell+1,i)
        P2 = Pmat(1:Msize, Msize+1:2*Msize); %P1ij = prob of getting absorbed to (ell+m,j) given the system starts from (ell+1,i)
        
        Delta = zeros(size(T));
        P_abs_ell = sum(P1,2); %element i is the probability of getting absorbed to level ell, given the system started from (ell+1,i)
        P_abs_lPlusm = sum(P2,2); %element i is the probability of getting absorbed to level ell + m, given the system started from (ell+1,i)

        Delta_mainDiag = diag(T).* P_abs_lPlusm./P_abs_ell; %it is column vec
        Delta_body =  P_abs_lPlusm * Delta_mainDiag';
        Delta = Delta_body - diag(diag(Delta_body)) + diag(Delta_mainDiag);

        if all(all(Delta < epsilon_r))
            continue_r = false(1);
        else
            m = m+1;
        end
    end
    Delta(abs(Delta)<eps) = 0;
    R_lowerbar = lambda2 * T;
    R_upperbar = lambda2 *(T + Delta);
end

