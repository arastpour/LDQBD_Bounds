%inplementation of (1)-(6) in Inversion of Block Tridiagonal Matrices Egil Kristoffer Gorm Hansen & Rasmus Koefoed Jakobsen Department of Computer Science – January 14th, 2010 Supervisor: Stig Skelboe 
function g1 = function_BlockTridiagonalInverse2(ell, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c, m)
    %I will calculate g1, which includes blocks g1j of matrix G  =A^-1  
    Msize = maxPhase + 1;    
    N = m-1;% The original matrix that we want to calculate its inverse is of size N*Msize by N*Msize    
    a11 = function_Q(1, ell+1, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c);
    d11L = a11;
    
    %calculate d11R
    CR = zeros(Msize, N*Msize);% I will put c_iR,   i = N,...,2 in this matrix. NOTE: for convenience, I set the size of DR and CR to be (Msize, N*Msize) and not (Msize, (N-1)*Msize). Basically, I assume that c1R = d11R = 0. But I update d11R later.
    %DR = zeros(Msize, N*Msize);% I will put d_iiR, i = N,...,2 in this matrix
    aNN       = function_Q(1, ell+N  , maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); 
    aNminus1N = function_Q(0, ell+N-1, maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); 
    dNNR = aNN;
    %DR(:,(N-1)*Msize+1 : N*Msize) = aNN; %dNNR
    CR(:,(N-1)*Msize+1 : N*Msize) = -aNminus1N / dNNR; %cNNR
    for i = N-1:-1:2
        a_i_minus1_i = function_Q(0, ell+i-1  , maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); 
        a_ii         = function_Q(1, ell+i    , maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c);   
        a_i_plus1_i  = function_Q(2, ell+i+1  , maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); 
        c_i_plus1_R  = CR(:,i*Msize+1 : (i+1)*Msize);
        d_ii_R = a_ii + c_i_plus1_R * a_i_plus1_i;
        c_i_R  = - a_i_minus1_i / d_ii_R;
        CR(:,(i-1)*Msize+1 : i*Msize) = c_i_R;
        %DR(:,(i-1)*Msize+1 : i*Msize) = d_ii_R;
    end
    c_2_R = CR(:,Msize+1 : 2*Msize);
    a_2_1  = function_Q(2, ell+2  , maxPhase, lambda1, lambda2, mu1, mu2, gamma1, gamma2, c); 
    d11R = a11 + c_2_R * a_2_1;    
    g11 = inv(-a11 + d11L + d11R);
    g1 = zeros(Msize, N*Msize);
    g1(:,1:Msize)=g11;
    for j = 2:N
        c_j_R = CR(:,(j-1)*Msize+1 : j*Msize);
        g1_j_minus1 = g1(:,(j-2)*Msize+1:(j-1)*Msize);
        g1(:,(j-1)*Msize+1 : j*Msize) = g1_j_minus1 * c_j_R;
    end
end
% to check:
% A = zeros( N*Msize , Msize);
% A (1:2*Msize,:) = vertcat( a11 , a_2_1 );
% g1 * A should be an Identity matrix of size Msize

% B = zeros( N*Msize , Msize);
% B ((N-2)*Msize+1:N*Msize,:) = vertcat( aNminus1N , aNN );
%g1 * B should be zero matrix of size Msize