function Qmatrix = function_Q(matNumber, ell, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c)
%Matrix Q0
    Q0 = lambda2 * eye(maxPhase+1);
    Q1 = zeros(maxPhase+1, maxPhase+1);
    Q2 = zeros(maxPhase+1, maxPhase+1);    
    
    if ell > 0
        for i=1:maxPhase+1 %Matrix Q2 
           ii = i -1; % class 1 count
           m = min (ii, c); % count of servers serving class 1
           e = c - m; % count of servers not serving class 1
           n = min (ell, e); %count of servers serving class 2
           maxx = max(ell - n, 0);% count of class 2 not getting service
           Q2(i, i) = n * mu2 + maxx * gamma2;
        end
    end
    
    %Matrix Q1. Note that we need to have Q0 and Q2 to calculate the diagonal elements of Q1
%     for  i = 1:maxPhase %above the diagonal
%         Q1(i, i+1) = (lambda1);
%     end
    Q1 = Q1 + diag(lambda1*ones(1,maxPhase),1);
        
    for i = 2:maxPhase+1 %below the diagonal
        ii = i -1; % class 1 count
        m = min (ii, c);  % count of servers serving class 1
        %e = c - m;
        %n = min (ell, e);
        maxx = max(ii - m, 0);
        Q1(i, i-1) = m * mu1 + maxx * gamma1; 
    end        
        
    for i = 1:maxPhase+1 %diagonal
        Q1(i, i) = -sum(Q0(i,:)) -sum(Q1(i,:)) -sum(Q2(i,:));
    end  
 
    %Function output:
    if matNumber == 0 %Matrix Q0 
        Qmatrix = Q0;
    elseif matNumber == 1 %Matrix Q1 
        Qmatrix = Q1;     
    elseif matNumber == 2 %Matrix Q2 
        Qmatrix = Q2;
    end



    
end







