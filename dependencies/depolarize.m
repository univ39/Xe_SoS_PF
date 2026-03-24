function [d_depol] = depolarize(d,flipvec,TR2T1)
%DEPOLARIZE takes into account flip angle train and corrects for
%depolarization in hyperpolarized acquisition, in dim=2
%   INPUTS
%             d   data to be corrected (#excitations,#points per exc)
%       flipvec   vector of flip angle
%         TR2T1   TR of each excitation/T1
%   OUTPUTS
%       d_depol   depolarized data
%
%   See depolarize_corr
%   02/2025 Kylie Yeung
if ~exist("flipvec","var");flipvec=[];end
if isempty(flipvec);flipvec=0;end

if ~exist("TR2T1","var");TR2T1=[];end
if isempty(TR2T1);TR2T1=0;end

if isscalar(flipvec);flipvec=repmat(flipvec,1,size(d,1));end

cos_factor=cumprod(cosd(flipvec).*exp(-TR2T1));
switch find(size(d) == length(flipvec))
    case 1
        factor=repmat([1 cos_factor(1:end-1)]',1,size(d,2));
    case 2
        factor=repmat([1 cos_factor(1:end-1)],size(d,1),1);
end
d_depol=d.*factor;
end

