function [B,D,S,info,options] = gp_BAJD_GPT(X,options)
% [B,D,S,info,options] = gp_CAJD_GPT(X,C,options)
% Compute bilinear joint diagonalization of matrices X and C. 
%
% INPUTS
% -----------------
%   X : are assumed to be non-symetrical matrices (possibly rectangular) [N x T x K]
%   options* is a structure with
%     .B0 is the initial B. default: eye(N)
%     .D0 is the initial D. default: randn(T,N)
%     .eps is the tolerance. default: 1e-12
%     .itMax is the maximum of iterations. default: 100
%     .A is the real left handed mixing matrix (if available). Only for performance
%     criterion
%     .E is the real right handed mixing matrix (if available). Only for performance
%     criterion
%
% OUTPUTS
% -----------------
%   B : the left handed unmixing matrix
%   D : the right handed unmixing matrix
%   S : the product B'*X*D 
%
% see also: gp_AJD_GPT, gp_CAJD_GPT
%
% *** History: 09-Fev-2016
% *** Author: Florent BOUCHARD & Louis KORCZOWSKI, GIPSA-Lab, 2016
% *** Reference: MINING THE BILINEAR STRUCTURE OF DATA WITH APPROXIMATE
%                JOINT DIAGONALIZATION
%                L. Korczowski, F. Bouchard, C. Jutten, M. Congedo
% *** Contact: louis.korczowski@gmail.com
% *** Licence: GNU GPLv3



% get variables from X and C
N = size(X,1); % number of channels
T = size(X,2); % number of samples
K = size(X,3); % number of X matrices

% initialize parameters according to options input
%%% default
B0 = eye(N);
D0 = orth(randn(T,N));
opt = struct('B0',B0,'D0',D0,'eps',1e-12,'itMax',100,'A',[],'E',[]);
%%% get final options
if nargin<2
    options = gp_getOpt([],opt);
else
    options = gp_getOpt(options,opt);
end

%%% get parameters
B0      = options.B0;
D0      = options.D0;
eps     = options.eps;
itMax   = options.itMax;
A       = options.A;
E       = options.E;

% init unmixing matrices
B = B0;
D = D0;

% init S according to B and D
S = zeros(N,N,K);
for k=1:K
    S(:,:,k) = B'*X(:,:,k)*D;
end

% init info structure
info = struct('it',0,'conv',nan,'B',B,'D',D,'critOFFS',gp_critOFF(S),'critA',gp_critMoAm(A,B),'critE',gp_critMoAm(E,D));

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
            beta = - sum(S(i,j,:).*S(j,j,:))/sum(S(j,j,:).*S(j,j,:));
            % get gamma
            gamma = - sum(S(j,i,:).*S(j,j,:))/sum(S(j,j,:).*S(j,j,:));
            % update B
            B(:,i) = B(:,i) + beta*B(:,j);
            % update D
            D(:,i) = D(:,i) + gamma*D(:,j);
            % update S
            S(i,:,:) = S(i,:,:) + beta*S(j,:,:);
            S(:,i,:) = S(:,i,:) + gamma*S(:,j,:);
            % update conv
            conv = conv + beta^2 + gamma^2;
        end
        end
    end
    % get conv
    conv = sqrt(conv/2/N/(N-1));
    % fill info struct
    info(it+1).it       = it;
    info(it+1).conv     = conv;
    info(it+1).B        = B;
    info(it+1).D        = D;
    info(it+1).critOFF = gp_critOFF(S);
    info(it+1).critA    = gp_critMoAm(A,B);
    info(it+1).critE    = gp_critMoAm(E,D);
end