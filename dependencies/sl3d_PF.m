%% Forward modelling of a pseudo3d trajectory
% Kylie Yeung 05/2025

function [Image,Image_zf,d_ramped_depol,F,E, BILD] = shepp_logan3d_PF(wf,FA,T1,F,E,nz_homodyne)
%% get variables
mtx_acq=wf.mtx;mtx_reco=mtx_acq;
npts_arm=size(wf.ks,2);
narms=size(wf.phi,2);
npts=npts_arm*narms;
nz=wf.mtx(3); % number of z encodes
n_fullctr=nz-nz_homodyne*2; % number of z encodes fully sampled in the center

k=wf.k;
if ~exist('FA','var');FA=[];end
if ~exist('F','var');F=[];end
if ~exist('E','var');E=[];end
if ~exist('T1','var')
    TR2T1=[];
else
    TR2T1=wf.t(npts_arm)./T1;
end

%% get encoding and reconstruction matrices
z = k(:,:,3);[~,i,~]=unique(z);
Kz = single(z(sort(i)))./0.5*pi./mtx_reco(3).*mtx_acq(3);
[F_z,k,dim,E_z] = getFz(k,mtx_acq, mtx_reco,npts);

if isempty(E);E=getE(k,mtx_acq,mtx_reco,dim,narms,npts);end
if isempty(F);F = getF(E,[],'cholesky');end

%% Generate a discretized Shepp-Logan test image
BILD=flipud(phantom3d('Modified Shepp-Logan',mtx_reco));
BILD=reshape(BILD,[mtx_reco(1)*mtx_reco(2) mtx_reco(3)]);

%% Forward encode the image into k-space data
d=E*(E_z*BILD.').';

%% undersample data in Z
ind=ones(1,mtx_reco(3));ind(2:2:2*nz_homodyne)=0;
d=bsxfun(@times,d,ind);

%% ramping and depolarizing
d_tmp=reshape(d,[npts_arm narms.*mtx_acq(3)]); % [#pts per exc, #exc]
if ~isempty(FA)
    d_tmp=bsxfun(@times, d_tmp, sind(FA)); % [#pts per exc, #exc]
end
d_ramped_depol=depolarize(d_tmp.',FA,TR2T1);

%% Partial Fourier
dxy=F*(reshape(d_ramped_depol.',[npts mtx_acq(3)]));
Kz_sorted=cat(2,Kz(:,1:2:end),flip(Kz(:,2:2:end),2));
dxy_sorted=cat(2,dxy(:,1:2:end),flip(dxy(:,2:2:end),2));
ZM=single(-(nz/2):nz/2-1);
E_z=exp(-1i*(Kz_sorted(:)*ZM(:).'));
% [im_homodyne] = PinvItPartialFourier(dxy_sorted.', E_z, 1, 0, 15, 0);
[im_homodyne] = homodyne_zPE(dxy_sorted.', [],15,E_z);

%% reconstruction
Image_zf=reshape((F_z*dxy.').',mtx_acq);
Image=reshape(im_homodyne.',mtx_acq);
end