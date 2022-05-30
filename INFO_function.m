function [I,M] = INFO_function(SList,X,p)
% Computes the mutual information I_X(S) and the MMSE matrix  M_X(s) for a
% discrete random vector supporte on the columns of X with pmf p

d = size(X,1);          % dimension of vector
K = size(X,2);          % number of atoms
numS = size(SList,3);   % number of SNR matrices

% Expectations are approximated using Monte Carlo integration. There is an
% outloop and and innerloop. Optimal allocation depends on the problem
% dimensions. If d and K are small, then is it best to put all of the
% itertions in the inner loop. 
numMC= 5;%1e2;             % number of trials (outer loop)
T = 1e3;                % number of trials (inner loop)


%initiallize MI and MMSE
I = zeros(1,numS);
M = zeros(d,d,numS);

for nS=1:numS
    S = SList(:,:,nS);
    S2 = sqrtm(S);
       
    Ch = zeros(d,d);
    Dh = 0;
    for itr = 1:numMC
        W = randn(d,T);
    
        % this four loop computed expectation w.r.t. px
        for k=1:K
            [Ci,Di] = CD_function(S2*X(:,k)*ones(1,T)+W,S,S2,X,p);
            %[Ci,Di] = CD_function_old(S2*X(:,k)+W(:,1),S,S2,X,p); % this
            %function only works for a single sample at a time
            Ch = Ch  + (p(k)/numMC)*Ci;
            Dh = Dh  + (p(k)/numMC)*Di;
        end
    end  
    M(:,:,nS) = Ch;
    I(nS) = Dh;
end

end


function [C,D] = CD_function(y,S,S2,X,p)
% Compute the average C(y) and D(y) for multiple instances of y

T = size(y,2);  % number of samples to be computed in parallel
K = size(X,2);  % number of atoms in distribution of X

% the vector u it the exponent of the unormalized probababilities. 
u = X'*S2*y - (1/2)*diag(X'*S*X)*ones(1,T); % u is KxT
L = max(u,[],1);                            % L is KxT
q = (p*ones(1,T)).*exp(u-ones(K,1)*L);      % q is KxT

R0 = sum(q,1);                              % R0 is 1xT 
Gavg =  sum(sum(u.*q,1)./R0)/T;             % Gavg is 1x1
D = Gavg - sum(L)/T - sum(log(R0))/T;       % D is 1x1 (averge over T)

q = q./R0;                         % normalize to get posterior prob.
Q = diag(sum(q,2))/T - q*q'/T;     % compute averge of probabilities
C = X*Q*X';                        % estimatoe of conditional covariance
end

function [C,D] = CD_function_old(y,S,S2,X,p)
% This function computes C(y) and D(y) for a single instance of y

u = X'*S2*y - (1/2)*diag(X'*S*X);
L = max(u);
q = p.*exp(u-L);

R0 = sum(q);
G = u'*q;
R1 = X*q;
R2 = X*diag(q)*X';

D = G/R0 - L - log(R0);
C = R2/R0 - R1*R1'/R0^2;
end