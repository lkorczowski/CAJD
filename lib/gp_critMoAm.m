function crit = gp_critMoAm(A,B)

% crit = gp_critMoAm(A,B)
%
% function computing the Moreau-Amari index for matrices A and B. This is
% usually used to estimate the performances of (joint) blind source
% separation algorithms when both the mixing and unmixing matrices are
% known.
%
% Inputs:
%   - A: NxNxM array containing true mixing matrices used to simulate data
%   - B: NxNxM array containing estimated unmixing matrices
%
% where N and M correspond to the number of channels and of datasets
% respectively.
%
% Output:
%   - crit: 1xM vector containing the Moreau-Amari index for each dataset
%
% reference: SELF-ADAPTIVE SOURCE SEPARATION BY DIRECT AND RECURSIVE
%            NETWORKS. O. Macchi, E. Moreau
%            Proc. International Conference on Digital Signal Processing
%            (DSP93), 1993.
%
% *** History: 18-Feb-2015
% *** Author: Florent BOUCHARD, GIPSA-Lab, 2015

crit=NaN;
 if any([isempty(A) isempty(B)])
%     return;
 else
%     % get variables
     N = size(A,2);
 M = size(A,3);
 end
%     
%     % normalize columns of A and B to avoid scale indeterminancy issues
%     A = gp_multinormc(A);
%     B = gp_multinormc(B);
%     
%     % compute product B'A
%     P = multiprod(multitransp(B),A);
%     
%     % rows
%     sr = sum(squeeze(sum(abs(P),1)./max(abs(P),[],1)-1));
%     % columns
%     sc = sum(squeeze(sum(abs(P),2)./max(abs(P),[],2)-1));
%     
%     crit = (sr+sc)/2/N/(N-1);
% end
% end

%% old code using for loops (without dependencies)
% % initialize output structure
% crit = zeros(1,M);
M = size(A,3);
 for m=1:M
    % normalize columns of A and B to avoid scale indeterminancy issues
    for i=1:N
        A(:,i,m) = A(:,i,m)/sqrt(A(:,i,m)'*A(:,i,m));
        B(:,i,m) = B(:,i,m)/sqrt(B(:,i,m)'*B(:,i,m));
    end
    % compute product B'A
    P = B(:,:,m)'*A(:,:,m);
    s = 0;
    % rows
    for i=1:N
        s = s + sum(abs(P(i,:)))/max(abs(P(i,:))) - 1;
    end
    % columns
    for i=1:N
        s = s + sum(abs(P(:,i)))/max(abs(P(:,i))) - 1;
    end
    crit(m) = s/2/N/(N-1);
 end
