
clearvars
formatOut = 'yyyy_mm_dd';
runDate = datestr(now,formatOut);
format longEng

%parameter inputs
mu1 = 1;
mu2 = 1; %10
gamma = 1;
c = 100;
rho = 0.95; 
propFromClass1 = 0.1; %it is $d$ in the paper
lambda1 =    propFromClass1 * rho * c * mu1; %to calculate lambda 1 and 2 basedon other parameters
lambda2 = (1-propFromClass1)* rho * c * mu2;
gamma1 = gamma; gamma2 = gamma;
epsilon_h = 10^(-6); %it is $\epsilon_p$ in the paper
epsilon_r = 10^(-10); 
error = 0.01; %it is $\epsilon_s$ in the paper
epsilon_aggrProb = 0.001; % It is $\epsilon_a$ in the paper
paramV(1,:) = horzcat (1, rho, propFromClass1, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, epsilon_r, epsilon_h);

cutOffPoint =60;%It is to preassign the size of x and R matrices. It does not actually need to be precise (as long as it is larger than 1). It is only to make sure that the size of upper_x_ellPrime_Matrix, lower_x_ellPrime_Matrix, upper_Pi_ellPrime_Matrix, and lower_Pi_ellPrime_Matrix do not change within the code. In case it is too small, the size of these matrices will start grwoing at each loop which will slow down the code. 
tau = inf; %it is $\ell^\tau$ in the paper
offeredLoad = (lambda1/mu1 + lambda2/mu2)/c; 
if offeredLoad >=1 && gamma1 == 0 && gamma2 == 0
    disp('The stability condition is not met because offeredLoad >=1 && gamma1 == 0 && gamma2 == 0');
    keyboard
end 
maxPhase = function_ErlangA_truncation_level(lambda1, mu1, gamma1, c, epsilon_h); %it is $p^\ast$ in the paper
lower_x_ellPrime_Matrix = zeros(cutOffPoint+1,maxPhase+1); %each row for a level, starting from 0,1,...,ellPrime. Each column is for a phase, starting from 0,1,...,maxPhase
upper_x_ellPrime_Matrix = inf * ones(cutOffPoint+1,maxPhase+1); 
lower_PI_ellPrime_Matrix = zeros(cutOffPoint+1,maxPhase+1); 
upper_PI_ellPrime_Matrix = ones(cutOffPoint+1,maxPhase+1);
indicatorV = zeros(cutOffPoint+1,1); %the ell-th element will indicate if the prob vector associated with level ell-1 satisfies the error tolerance 
indicatorM = zeros(cutOffPoint+1,cutOffPoint+1); %I use this matirx just for debugging. At iteration k, I put indicatorV associated with this iteration in column k of this matrix
startTime = tic;
k = 0;    
m_initial = 4;
[lower_R0,upper_R0,m]=function_R_low_and_up(epsilon_r,0, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, m_initial);
if m == -1, textToReport = 'Matrix W2 was close to singular'; function_writeToExcel(textToReport,lineNumber,runDate,paramV(problem,:),namesV); keyboard, end
med_logic = min(min(upper_R0 - lower_R0)) < - eps || any(any(lower_R0 < -eps)); % I use eps instead of 0 to account for the floating-point relative error
if med_logic
    disp (['ERROR: Matrix W2 is close to singular;  min(min(upper_R0 - lower_R0)) = ', num2str(min(min(upper_R0 - lower_R0))), '     min(min(lower_R0)) = ', min(min(lower_R0))] )%Note that depending on specs of your machine, you may get this message on some machines and not on others.
    keyboard
end
lower_R0(lower_R0<eps) = 0;
Q_1_0 = function_Q(1, 0, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); %$Q_1^(0)$
Q_2_1 = function_Q(2, 1, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); %$Q_2^(1)$    
[lower_x0, upper_x0, varDrop] = function_x0(Q_1_0, Q_2_1, lower_R0, upper_R0, lambda1, mu1,gamma1,c,maxPhase,lambda2, mu2, gamma2, m, epsilon_r);
lower_x0(lower_x0<eps) = 0;
lower_c = sum(lower_x0);    upper_c1 = sum(upper_x0); upper_c2 = inf;  upper_c = upper_c1 + upper_c2;
lower_x_ellPrime_Matrix(1,:) = lower_x0;        upper_x_ellPrime_Matrix(1,:) = upper_x0;
lower_PI_ellPrime_Matrix = lower_x_ellPrime_Matrix ./upper_c;
upper_PI_ellPrime_Matrix = upper_x_ellPrime_Matrix ./lower_c;

if min(upper_PI_ellPrime_Matrix - lower_PI_ellPrime_Matrix,[],'all')< -eps || any (lower_PI_ellPrime_Matrix < -eps,'all')
    disp('ERROR: the PI_0 inequality does not hold.')
    errorCount = errorCount +1;
end
lower_PI_ellPrime_Matrix(lower_PI_ellPrime_Matrix<eps) = 0;
lower_Rk_minus1 = lower_R0;           upper_Rk_minus1 = upper_R0;
lower_xk_minus1 = lower_x0;           upper_xk_minus1 = upper_x0;
IterTime = toc(startTime);

lowerBoundSumV = -1;
lower_xk_plus_1 = []; upper_xk_plus_1 = [];  lower_Rk_plus_1 = []; upper_Rk_plus_1 = [];
while  sum(lowerBoundSumV) < 1 - epsilon_aggrProb 
    iterationTIme = tic;
    k = k+1;
    if k < tau
        lower_xk = lower_xk_minus1 * lower_Rk_minus1;
        upper_xk = upper_xk_minus1 * upper_Rk_minus1;
        m_initial = m;
        [lower_Rk,upper_Rk,m]=function_R_low_and_up(epsilon_r,k, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, m_initial);
        med_logic = min(min(upper_Rk - lower_Rk)) < - eps || any(any(lower_Rk < -eps)); % I use eps instead of 0
        if med_logic
            disp (['min(min(upper_Rk - lower_Rk)) = ', num2str(min(min(upper_Rk - lower_Rk))), '     min(min(lower_Rk)) = ', min(min(lower_Rk))] )
            upper_PI_ellPrime = 0; lower_PI_ellPrime = 0;
            keyboard
        end
        lower_Rk(lower_Rk<eps) = 0;
    end

    tau_condition = function_tau_cond (k,maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); % this function returns True if the condition is met
    if (tau == inf && k >= c) && tau_condition
        tau = k;
    end
    lower_c = lower_c + sum(lower_xk);        
    upper_c1 = upper_c1 + sum(upper_xk);
    upper_c2 = inf;
    if  k >= tau
        m_initial = m;
        [lower_Rk_plus_1,upper_Rk_plus_1,m]=function_R_low_and_up(epsilon_r,k+1, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, m_initial);  
        lower_Rk_plus_1(lower_Rk_plus_1<eps) = 0;
        lower_xk_plus_1 = lower_xk * lower_Rk;    upper_xk_plus_1 = upper_xk * upper_Rk;

        eigenValues = eig(upper_Rk_plus_1);  
        if min(eigenValues) > 0 && max(eigenValues) < 1
            medianMatrix = eye(maxPhase+1) - upper_Rk_plus_1;
            upper_c2 = sum(upper_xk_plus_1 / medianMatrix);   %upper_xk_plus_1 / medianMatrix is equivalet to upper_xk_plus_1 * medianMatrix^(-1) ;       
        end
    end

    upper_c = ((upper_c1 + upper_c2));
    lower_x_ellPrime_Matrix(k+1,:) = lower_xk; % I record the bounds for all x_ellPrime (unnormalized probs), but I will use only those that satisfy the torlerance error for calculating performance measures
    upper_x_ellPrime_Matrix(k+1,:) = upper_xk;
    lower_PI_ellPrime_Matrix = lower_x_ellPrime_Matrix ./upper_c;
    upper_PI_ellPrime_Matrix = upper_x_ellPrime_Matrix ./lower_c;

    if min(upper_PI_ellPrime_Matrix - lower_PI_ellPrime_Matrix,[],'all') < -eps || any(lower_PI_ellPrime_Matrix < -eps,'all')
        disp('ERROR: the PI_0 inequality does not hold. 2')
        keyboard
    end    
    lower_PI_ellPrime_Matrix(lower_PI_ellPrime_Matrix<eps) = 0;

    diffMat = abs(upper_PI_ellPrime_Matrix - lower_PI_ellPrime_Matrix) < error; % it is the matrix of differences between lower and upper bounds for all levels.
    indicatorV = all(diffMat,2); %if all elements of a row of diffMat is 1 (meaning that pi(level, phase), for all phase values, is within the tolerence error) then the indicator becomes 1. Otherwise it is 0.
    firstZeroInd = find(~indicatorV,1); %to find the index of the first element in indicatorV that is zero. In other words, for all levels below this index, the tolerance level is satisfied. NOTE: firstZeroInd - 1 is the index of the last non-zero element in indicatorV. And element i in indicatorV is associated with level i-1. For example, the first element is for level 0. Therefore, firstZeroInd - 2 is the highest level that we have its probs within the tolerance error
    if (~isempty(firstZeroInd)  && firstZeroInd > 1 )   %sometimes, it is possible that, for example, pi for levels 3-10 is within the tolerence error, but pi for 0-2 are not. firstZeroInd gives me the first level for which the tolerence level is NOT satisfied. I calculate the sum of probs for all levels below that level and check if the sum is larger than 1 - epsilon_aggrProb. If it is, then I go out of the while loop. 
        %all(indicatorV) means that all levels satisfy the prob.
        %tolerence level and there are NO zeros in the index vector.
        %(~isempty(firstZeroInd)  && firstZeroInd > 1 ) means that
        %there is at least one level that does not satisfy the tol. and
        %that level is not level 0.
        lowerBoundSumV = sum(lower_PI_ellPrime_Matrix(1:firstZeroInd-1,:),2);
    elseif all(indicatorV)
        lowerBoundSumV = sum(lower_PI_ellPrime_Matrix,2);
    else
        lowerBoundSumV = -1;
    end        
    lower_xk_minus1 = lower_xk;           upper_xk_minus1 = upper_xk;
    lower_Rk_minus1 = lower_Rk;           upper_Rk_minus1 = upper_Rk; 
    lower_xk        = lower_xk_plus_1;    upper_xk        = upper_xk_plus_1;
    lower_Rk        = lower_Rk_plus_1;    upper_Rk        = upper_Rk_plus_1;       

    disp (['k = ', num2str(k),  ', sum(lowerBoundSumV) = ', num2str(sum(lowerBoundSumV)),...
        ', elapsed time (sec) = ', num2str(toc(startTime)), ', Iteration time = ', num2str(IterTime)])        
    IterTime = toc(iterationTIme);

    tookTooLong = false(1);
    hours = 5;
    if toc(startTime)>3600*hours
        lowerBoundSumV = 1; %it couses that the while loop terminates
        tookTooLong = true(1);
    end
end
    
if tookTooLong
    disp(['It has been running for more than ', num2str(hours),'. So I am going to terminate the code and going to go to the next problem.'])
    keyboard
else
    ttime = toc(startTime);
    if isempty(firstZeroInd)
        upperBmat = upper_PI_ellPrime_Matrix;
        lowerBmat = lower_PI_ellPrime_Matrix;
    else
        upperBmat = upper_PI_ellPrime_Matrix(1:firstZeroInd-1,:);
        lowerBmat = lower_PI_ellPrime_Matrix(1:firstZeroInd-1,:);
    end
    f = zeros (size(lowerBmat));
    for level = 1:size(lowerBmat,1)
        for phase = 1:size(lowerBmat,2)
              ell = level -1;
              p = phase - 1;
              f(level,phase) = function_class2waiting(ell,p,c);

        end
    end
    averageProbs = (lowerBmat + upperBmat)/2;
    performanceMeasure = sum(sum(f .* averageProbs,2));
    disp(['performanceMeasure = ', num2str(performanceMeasure), ', total elapsed time (sec) = ',  num2str(toc(startTime))])
    calculateBounds = 'true'; % As discussed in the paper, we can calculate the bounds for some performance measures and we cannot do it for some others. We manually set calculateBounds to true and false to indicate whether we want to calculate the bounds or not.
    if calculateBounds
        boundMat =  (upperBmat - lowerBmat)/2;
        performanceMeasureBound = sum(sum(f .* boundMat,2));
        performanceMeasureErrorLowerBound = - performanceMeasureBound;
        performanceMeasureErrorUpperBound = performanceMeasureBound + epsilon_aggrProb;
        disp(['performanceMeasureErrorLowerBound = ', num2str(performanceMeasureErrorLowerBound), ', performanceMeasureErrorUpperBound = ', num2str(performanceMeasureErrorUpperBound)])
    end
end
    


