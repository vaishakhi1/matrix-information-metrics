function HY = Entropy_3comm_corrected(aXList,mucell)
u1 = [sort(-logspace(-2,2,30)) logspace(-2,2,8)];
u2 = u1;
 ngroups = length(aXList);
 HY = zeros(ngroups,1);
% mudiff = zeros(1,2,9);
%  for l = 1:3
%      for k = 1:3
%          mudiff(:,:,(l-1)*3+k) = mucell{k} - mucell{l};
%          all((l-1)*3+k) = aXList(k);
%      end
% end
 
%    A = @(u1,u2,i,k) -(0.5* sum( (repmat([u1;u2],1,1,9) - mudiff).^2 ));
 
for l=1:3

        clear mudiff
     for k = 1:3
         mudiff(:,:,k) = mucell{k} - mucell{l};
         all(:,:,k) = aXList(k);
     end
     
     A = @(u1,u2) bsxfun(@plus,bsxfun(@plus, -(0.5* sum( (bsxfun(@minus,[u1;u2], mudiff)).^2 )), log(all)),1);% + log(all) + 1;
         %t = (eye(2)+tcell{i1}+tcell{k1});
         %m = mucell{k1} - mucell{i1};
    pmax = @(u1,u2) max(A(u1,u2),[],3); 

%         A = @(u1,u2,i,k) log(realmin+sqrt(2*pi)*mvnpdf([u1;u2]',(mucell{k} - mucell{i})',(eye(2)+tcell{i}+tcell{k}))');%-diag(0.5*(([u1;u2]-m)*t*([u1;u2]-m)))
%     A = @(u1,u2,i,k) -(0.5* (    (u1 - mucell{k}(1) + mucell{i}(1)).*(u1 - mucell{k}(1) + mucell{i}(1)) ...
%         + (0+tcell{i}(1,2)+tcell{k}(1,2))*(u2 - mucell{k}(2) + mucell{i}(2)).*(u1 - mucell{k}(1) + mucell{i}(1))...
%         + (0+tcell{i}(2,1)+tcell{k}(2,1))*(u1 - mucell{k}(1) + mucell{i}(1)).*(u2 - mucell{k}(2) + mucell{i}(2))...
%         + (1+tcell{i}(2,2)+tcell{k}(2,2))*(u2 - mucell{k}(2) + mucell{i}(2)).*(u2 - mucell{k}(2) + mucell{i}(2))    ));

   phi = @(u1,u2) 1/(2*pi)*exp(-0.5*sum([u1;u2].^2));

    
    %pmax = @(u1,u2) max([ A(u1,u2,1,2)', A(u1,u2,1,1)', A(u1,u2,1,3)',  A(u1,u2,2,2)', A(u1,u2,2,1)', A(u1,u2,2,3)',  A(u1,u2,3,2)', A(u1,u2,3,1)', A(u1,u2,3,3)'],[],2)'; 
    
       p2 = @(u1,u2) exp(bsxfun(@minus, A(u1,u2), pmax(u1,u2)));

    py =  @(u1,u2) pmax(u1,u2) + log(sum(p2(u1,u2),3));
                
                %pz = @(u1,u2,i,k) (aXList(i)*aXList(k))*p1(u1,u2,i,k).*(pmax(u1,u2) + log(py(u1,u2))); 
      HY(l) = -integral2(@(u1,u2) aXList(l)*phi(u1,u2).*py(u1,u2),-inf,inf,-inf,inf,'RelTol',1e-10,'method','iterated');
%      blah = 1
 end
end

