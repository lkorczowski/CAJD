function [B,C,info,options] = gp_AJD_GPT(C,options)
% [B,D,C,info,options] = gp_AJD_GPT(C,options)
% Compute approximate joint diagonalization of matrices C. 
%
% INPUTS
% -----------------
%   C : are assumed to be symetrical matrices [N x N x K]
%   options* is a structure with
%     .B0 is the initial B. default: eye(N)
%     .eps is the tolerance. default: 1e-12
%     .itMax is the maximum of iterations. default: 100
%     .A is the real left handed mixing matrix (if available). Only for performance
%     criterion
%
% OUTPUTS
% -----------------
%   B : the left handed unmixing matrix
%   C : the product B'*C*B
%
% see also: gp_CAJD_GPT, gp_BAJD_GPT
%
% *** History: 09-Fev-2016
% *** Author: Florent BOUCHARD & Louis KORCZOWSKI, GIPSA-Lab, 2016
% *** Reference: MINING THE BILINEAR STRUCTURE OF DATA WITH APPROXIMATE
%                JOINT DIAGONALIZATION
%                L. Korczowski, F. Bouchard, C. Jutten, M. Congedo
% *** Contact: louis.korczowski@gmail.com
% *** Licence: GNU GPLv3



% get variables from X and C
N = size(C,1); % number of channels
K = size(C,3); % number of matrices

% initialize parameters according to options input
%%% default
B0 = eye(N);
opt = struct('B0',B0,'eps',1e-12,'itMax',100,'A',[]);
%%% get final options
if nargin<2
    options = gp_getOpt([],opt);
else
    options = gp_getOpt(options,opt);
end

%%% get parameters
B0      = options.B0;
eps     = options.eps;
itMax   = options.itMax;
A       = options.A;

% init unmixing matrix
B = B0;

% init C according to B
Ctmp = zeros(N,N,K);
for k=1:K
    Ctmp(:,:,k) = B'*C(:,:,k)*B;
end

% init info structure
info = struct('it',0,'conv',nan,'B',B,'critOFF',gp_critOFF(C),'critA',gp_critMoAm(A,B));

% Main loop
conv   = 1;
it     = 0;
while it<itMax && conv>eps
    it   = it+1;
    conv = 0;
    for i=1:N
        for j=1:N
        if  j~=i
            % get beta
            beta  = - sum(Ctmp(i,j,:).*Ctmp(j,j,:))/sum(Ctmp(j,j,:).*Ctmp(j,j,:));
            % update B
            B(:,i) = B(:,i) + beta*B(:,j);
            % update C
            tmp1 = Ctmp(i,i,:); tmp2 = Ctmp(i,j,:); tmp3 = Ctmp(j,j,:);
            Ctmp(i,:,:) = Ctmp(i,:,:) + beta*Ctmp(j,:,:);
            Ctmp(:,i,:) = Ctmp(:,i,:) + beta*Ctmp(:,j,:);
            Ctmp(i,i,:) = tmp1 + 2*beta*tmp2 + beta^2*tmp3; % diagonal elements
            % update conv
            conv = conv + beta^2;
        end
        end
    end
    % get conv
    conv = sqrt(conv/2/N/(N-1));
    % fill info struct
    info(it+1).it      = it;
    info(it+1).conv    = conv;
    info(it+1).B       = B;
    info(it+1).critOFF = gp_critOFF(Ctmp);
    info(it+1).critA   = gp_critMoAm(A,B);
end
C = Ctmp;