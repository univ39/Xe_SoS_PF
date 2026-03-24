function [d_corr] = depolarize_corr(d,flipvec,TR2T1)
%DEPOLARIZE_CORR takes into account flip angle train and corrects for
%depolarization in hyperpolarized acquisition, in dim=2
%   INPUTS
%             d   data to be corrected
%       flipvec   vector of flip angle
%   OUTPUTS
%        d_corr   corrected data
%
%   07/2024 Kylie Yeung
if isscalar(flipvec);flipvec=repmat(flipvec,1,size(d,1));end
if ~exist("TR2T1","var");TR2T1=[];end
if isempty(TR2T1);TR2T1=0;end

cos_factor=cumprod(cosd(flipvec).*exp(-TR2T1));
switch find(size(d) == length(flipvec))
    case 1
        b1_factor=repmat([1 cos_factor(1:end-1)]',1,size(d,2));
    case 2
        b1_factor=repmat([1 cos_factor(1:end-1)],size(d,1),1);
end
d_corr=d./b1_factor;
end

