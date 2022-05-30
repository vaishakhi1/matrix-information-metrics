function [I,M] = INFO_function_2GMM(s,aXList,mucell,sz)
% This function compute the mutual information I_X(s) and the MMSE M_X(s)
% for a Gaussian mixture with paramters aXList, muXList, and tauX2List

% Find number of mixtures and make sure that all distrbution paramters are
% column vectors.
%reshape(aXList, [d,1]);
% muXList = reshape(muXList, [d,1]);
% tauX2List = reshape(tauX2List, [d,1]);
% tauXList = sqrt(tauX2List);
% 
% 
% d = real(eig(s));
% dpos = sum(d<0);
% dreal = isreal(d);
% if dreal==0 || dpos~= 0
%     I = 0;
% M = zeros(2);
%     
% else
%     
%     
%     Q = ((a-b)*eye(sz) + b*ones(sz,sz));
%     qtilde = sqrtm(Q);
%     
%     %compute exp
%     for i = 1:3
%         % tcell{i} = stilde*stilde*tcell1{i};
%         mucell{i} = qtilde*mucell1{i};
%     end
%     cell2mat(mucell);
%     EX = sum(repmat(aXList,sz,1).*cell2mat(mucell),2);
%     CovX = zeros(sz);
%     for i = 1:3
%         CovX = CovX+aXList(i)*(mucell{i}-EX)*(mucell{i}-EX)';
%     end
%     [V,D] = eig(CovX);
%     if (D(1,1) <= 0 || D(2,2) <= 0)
%         stophere = 1;
%     end
%     stdX = (V*sqrt(D)*V');
%     for i = 1:3
%         Munew(:,i) = stdX^(-1)*(mucell{i});
%     end
%     mucell = mat2cell(Munew, [sz],[1,1,1]);
    
    % Q = [a b; b a];
%    qtilde = sqrtm(Q);
%     lQ = length(Q);
    basis = eye(sz);
    gam = 1e-6;
    
    %CovX = [aXList;aXList].*(cell2mat(mucell1) - EX)
    
    
    % u1 = [sort(-logspace(-8,3,500)) logspace(-8,3,500)];
    % u2 = u1;
    
    stilde = sqrtm(s);%sqrtm(qtilde*s*qtilde');
    for i = 1:3
        %    tcell{i} = stilde*stilde*tcell1{i};

        mucell1{i} = stilde*mucell{i};
    end
    
    ngroups = length(aXList);
    nbases = length(basis);
    HY = sum(Entropy_3comm_corrected(aXList,mucell1));
    delHY = zeros(nbases,nbases);
    
    HYX = 0.5*log(det(2*pi*exp(1)*eye(2)));
    
    I = real(HY);
    
    
    for i = 1:nbases
        for j = 1:nbases
            ei = basis(:,i);
            ej = basis(:,j);
            %delrt = sqrt(gam)*(abs(s).*(ei*ej'));
            del = 0.5*(gam)*((ei*ej')) + 0.5*(gam)*((ej*ei'));
            stilde = sqrtm(s+del);% + (delrt);
            %stilde = sqrtm(s);%sqrtm(qtilde*sgam*qtilde');
            for l = 1:3
                %tcell{l} = stilde*stilde*tcell1{l};
                mucell1{l} = stilde*mucell{l};
            end
            ytemp(:,2*(i-1)+j) = Entropy_3comm_corrected(aXList,mucell1);
            delHY(i,j) = sum(Entropy_3comm_corrected(aXList,mucell1));
            M(i,j) = 2*(delHY(i,j) - I)./(del(i,j));
        end
        
    end
    M = [1 0.5; 0.5 1].*M;
    %M = 0.5*(M+M');
    %M = qtilde*M*qtilde';
    %M = zeros(2);
    % HYX = 0.5*log(det(2*pi*exp(1)*eye));
    %
    % I = HY - HYX;
    %M = 2*abs(delHY - HY)/gam;
end



% M1 = 2*(Ia-Ib)./(gamma);
% M2 = 2*(Ia-Ic)./gamma;

% monte carlo approach.
