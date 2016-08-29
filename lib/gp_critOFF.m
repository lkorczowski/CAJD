function crit = gp_critOFF(M)
% crit = gp_critOFF(M)
% function computing a non diagonality measure of the set of matrices M.
% Matrix representations are supposed to be over the two first dimensions.
%
% Input:
%   - M: multi-dimensional array of at least two dimensions
%
% Output:
%   - crit: value of the non diagonality measure
%
% *** History: 05-Jan-2015
% *** Author: Florent BOUCHARD, GIPSA-Lab, 2015

  sz   = size(M);
 nDim = length(sz);
 if nDim<2,       error('input are not matrices');      end;
 % if sz(1)~=sz(2), error('input matrices are not square'); end;
 
% reshape data so it is 3-D
 p = prod(sz)/sz(1)/sz(2);
 M = reshape(M,[sz(1),sz(2),p]);
 
 % new version (with dependencies)

% 
% % get diag and off parts of M
% [~,M_d,M_o] = gp_multidiag(M);
% 
% % compute squares of frobenius norms of matrices in M_d and M_o
% n2_d = multitrace(multiprod(M_d,multitransp(M_d)));
% n2_o = multitrace(multiprod(M_o,multitransp(M_o)));
% 
% % compute output
% tmp = n2_o./n2_d;
% crit = sum(tmp)/p;
% % normalize to obtain measure between 0 and 1
% crit = crit*min(sz(1:2))/(sz(1)*sz(2)-min(sz(1:2)));
% end

%% old version with for loops (without dependencies)
% % compute crit
% sz   = size(M);
cOFF = 0;
for i=1:p
    tmp_d = diag(diag(M(:,:,i)));
    tmp_o = M(:,:,i) - tmp_d;
    
    cOFF = cOFF + trace(tmp_o*tmp_o')/trace(tmp_d*tmp_d')/p;

end
     crit = cOFF*min(sz(1:2))/(sz(1)*sz(2)-min(sz(1:2)));
