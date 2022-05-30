function [I,M] = INFO_function_GMM(sList,aXList,muXList,tauX2List)
% This function compute the mutual information I_X(s) and the MMSE M_X(s)
% for a Gaussian mixture with paramters aXList, muXList, and tauX2List

% Find number of mixtures and make sure that all distrbution paramters are
% column vectors. 
d = length(aXList);    
aXList = reshape(aXList, [d,1]);
muXList = reshape(muXList, [d,1]);
tauX2List = reshape(tauX2List, [d,1]);
tauXList = sqrt(tauX2List);


phi = @(u) 1/sqrt(2*pi)*exp(-u.^2/2);

I = zeros(size(sList));
I2 = zeros(size(sList));

gamma = 1e-7;
sList2 = sList*(1+gamma);

% do integral seperately for each mixture. One could combine the integrals
% but this can lead to numeical issues. 
for k=1:length(aXList)
    
    % parameters of the Y vector
    B = (muXList-muXList(k))*sqrt(sList).*(1 + tauX2List*sList).^(-1/2);
    A = (1 + tauX2List(k)*ones(size(tauX2List))*sList).^(1/2).*(1 + tauX2List*sList).^(-1/2);
    C = -(1/2)*log(1 + tauX2List*sList) + 1/2 +log(aXList)*ones(size(sList));
    
    P = @(u) -1/2*(u*A - B).^2 + C;
    pmax = @(u) max(P(u));
    g = @(u) pmax(u) + log(sum(exp(P(u)-ones(size(aXList))*pmax(u))));    
    
    % parameters of the Y vector
    B2 = (muXList-muXList(k))*sqrt(sList2).*(1 + tauX2List*sList2).^(-1/2);
    A2 = (1 + tauX2List(k)*ones(size(tauX2List))*sList2).^(1/2).*(1 + tauX2List*sList2).^(-1/2);
    C2 = -(1/2)*log(1 + tauX2List*sList2) + 1/2 +log(aXList)*ones(size(sList2));
    
    I = I - aXList(k)*integral(@(u) phi(u).*g(u),-inf,inf,'ArrayValued',true,'RelTol',1e-10);

    P2 = @(u) -1/2*(u*A2 - B2).^2 + C2;
    pmax2 = @(u) max(P2(u));
    g2 = @(u) pmax2(u) + log(sum(exp(P2(u)-ones(size(aXList))*pmax2(u))));  
    
    I2 = I2 - aXList(k)*integral(@(u) phi(u).*g2(u),-inf,inf,'ArrayValued',true,'RelTol',1e-10);
    
end
n2 = 1;
M = 2*(I2-I)./(gamma*sList);
end
