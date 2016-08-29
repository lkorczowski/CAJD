function [C,Ct] = convertDat_BAJD(X,Cf)

N = size(X,1);
T = size(X,2);
K = size(X,3);

% get Ct
Ct = zeros(N,N,K);
trat=0;
for k=1:K
    Ct(:,:,k) = cov(X(:,:,k)');
    trat=trat+trace(Ct(:,:,k));
end
traf=0;
for f=1:size(Cf,3)
   traf=traf+trace(Cf(:,:,f));
end
% Ct = mean(Ct,3);

C = cat(3,trat*Cf/traf,Ct);