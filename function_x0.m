function [lower_x0, upper_x0,varDrop] = function_x0(Q_1_0, Q_2_1, lower_R0, upper_R0, lambda1, mu1,gamma1,c,maxPhase,lambda2, mu2, gamma2,m, epsilon_r)
%imshow('Help3.jpg')                     
    %In these lines, I find the index of the element, for which the
    %distance between its prob.^-1 and the average of the distances is
    %minimum. I used to drop the last line all the time but it caused
    %difficutlties because the output bounds on x0 were from VERY different
    %orders of magnitude.
    %upper_PI_ellPrime = function_MMCK(lambda1, mu1,gamma1,c,maxPhase);
    %med1 = upper_PI_ellPrime.^-1;
    %med2 = mean(med1);
    %med3 = abs(med1 - med2*ones(size(med1)));
    %varDrop = find(med3 == min(med3));
    varDrop=-4;%I set it to -4 because my step size below is 5 and I want varDrop to start from -4+5 = 1 and then continue  with 6, 11,...
    varDropNeedsUpdate_Lower= true(1); %it is to check if the current varDrop satisfies CONDITIONN1 below, which isrelated to lowerbound calcs
    varDropNeedsUpdate_Upper = true(1); %it is to check if the current varDrop satisfies CONDITIONN2 below, which isrelated to upperbound calcs
    while varDropNeedsUpdate_Lower || varDropNeedsUpdate_Upper
        varDrop = min(varDrop + 5,maxPhase);
        varDropNeedsUpdate_Lower = false(1); %If varDrop is fine, then these falses will cause that the code gets out of the whilse loop. Otherwise, it will comeback and will go throught the while loop with the updated varDrop 
        varDropNeedsUpdate_Upper = false(1);
        lower_z = lower_R0 * Q_2_1;
        lower_k = Q_1_0 + lower_z;

        lower_k_cut = lower_k;
        lower_k_cut (varDrop, :) = [];
        lower_k_cut (:, varDrop) = [];    

        lower_b_cut = - lower_k(varDrop, :);
        lower_b_cut (:, varDrop) = [];



        lastwarn('') % it is to reset lastwarn
        warnmsgSingularityTest = 0; msgidSingularityTest = 0;
        lower_k_cut(abs(lower_k_cut)<eps) = 0;
        lower_x0 = lower_b_cut / lower_k_cut;
        [warnmsgSingularityTest, msgidSingularityTest] = lastwarn;
         CONDITIONN1 =  strcmp(msgidSingularityTest,'MATLAB:singularMatrix') || min(lower_x0) < -eps || strcmp(msgidSingularityTest,'MATLAB:nearlySingularMatrix');  %check if upper_k_cut is invertible. And if upper_k_cut is inverse positive. check the eps fucntion to see why I am usin this instead of 0.
        if  CONDITIONN1
            varDropNeedsUpdate_Lower= true(1);
            if varDrop == maxPhase
               disp('the equation systems cannot be solved due to numerical issues.');
               keyboard
            end            
        end
        lower_x0(abs(lower_x0)<eps) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~varDropNeedsUpdate_Lower %it is not to waste time incase varDropNeedsUpdate_Lower is true. In this case, I want to go back to the while loop, update varDrop and do the calcas from the beginning
            upper_z = upper_R0 * Q_2_1;
            upper_k = Q_1_0 + upper_z;

            upper_k_cut = upper_k;
            upper_k_cut (varDrop, :) = []; 
            upper_k_cut (:, varDrop) = [];  

            upper_b_cut = - upper_k(varDrop, :);
            upper_b_cut (:, varDrop) = [];     

            %we conjecture that upper_k_cut is inverse positive. Here, I check it,
            %and if it does not hold, I recalculate upper_R0 with a larger m. I
            %increase m until the condition holds.

            lastwarn('') % it is to reset lastwarn
            warnmsgSingularityTest = 0; msgidSingularityTest = 0;
            upper_k_cut(abs(upper_k_cut)<eps)=0;
            upper_x0 = upper_b_cut / upper_k_cut;
            [warnmsgSingularityTest, msgidSingularityTest] = lastwarn;
            CONDITIONN2 =  min(upper_x0 - lower_x0)< -eps  || strcmp(msgidSingularityTest,'MATLAB:singularMatrix') || min(upper_x0) < -eps || strcmp(msgidSingularityTest,'MATLAB:nearlySingularMatrix'); %check if upper_k_cut is invertible. And if upper_k_cut is inverse positive. check the eps fucntion to see why I am usin this instead of 0.
            diff2V = [100 , 100];
            while CONDITIONN2
                disp(['The conjecture does not hold! fucntion_x0---m = ',num2str(m), '   min(upper_x0 - lower_x0) = ', num2str(min(upper_x0 - lower_x0)), '   maxPhase = ', num2str(maxPhase)])
                m = m+200;
                %upper_R0     = function_R_up (m,         0, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c);
                m_initial = m;                  
                [lower_R0, upper_R0, m] = function_R_low_and_up(epsilon_r, 0, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, m_initial);
                upper_z = upper_R0 * Q_2_1;
                upper_k = Q_1_0 + upper_z;
                upper_k_cut = upper_k;
                upper_k_cut (varDrop, :) = [];
                upper_k_cut (:, varDrop) = [];
                lastwarn('') % it is to reset lastwarn
                warnmsgSingularityTest = 0; msgidSingularityTest = 0;
                upper_k_cut(abs(upper_k_cut)<eps)=0;
                upper_x0 = upper_b_cut / upper_k_cut;
                [warnmsgSingularityTest, msgidSingularityTest] = lastwarn;            
                CONDITIONN2 =  min(upper_x0 - lower_x0)< -eps  || strcmp(msgidSingularityTest,'MATLAB:singularMatrix') || min(upper_x0) < -eps || strcmp(msgidSingularityTest,'MATLAB:nearlySingularMatrix'); %check if upper_k_cut is invertible. And if upper_k_cut is inverse positive. check the eps fucntion to see why I am usin this instead of 0.
                %CONDITIONN2 = min(upper_x0 - lower_x0)< -eps  || strcmp(msgidSingularityTest,'MATLAB:singularMatrix') || min(upper_x0) < -eps); %check if upper_k_cut is invertible. And if upper_k_cut is inverse positive. check the eps fucntion to see why I am usin this instead of 0.
                diff1V = diff2V;
                diff2V = [min(upper_x0 - lower_x0) ,  min(upper_x0)];
                diff = abs(diff2V - diff1V);
                if max(diff) < eps %here, I check if adding to m actually helps. If not, I go back and update varDrop
                    varDropNeedsUpdate_Upper = true(1);
                    CONDITIONN2 = false(0);%to jump out of this while loop, go back up, update varDrop and do everything again.
                    if varDrop == maxPhase
                       disp('UPPER: I have exhausted all possible values for varDrop but it seems it did not work!')
                       there is something wrong here
                    end
                end
            end
        end
    end
    %upper_x0(abs(upper_x0)<eps) = 0;
    lower_x0 = cat(2, lower_x0(1:varDrop-1), 1, lower_x0(varDrop:end));
    upper_x0 = cat(2, upper_x0(1:varDrop-1), 1, upper_x0(varDrop:end));
end

