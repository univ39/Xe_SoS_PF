function [im, ksp, MSEs] = PinvItPartialFourier_xyz(data, Exy, Fxy, Ez, mtx,its, zf, n, dir, refim)
%Iterative reconstruction for Cartesian data in 1 dimension
%
%       data                                data e.g. [data_xy, data__z, coils, passes]
%       its                                   number of iterations
%       zf                                    number of zero-filled outer z-encodes (at start and end of z) to ignore
%       n                                      number of z-encodes toremove and reconstruct
%       dir                                   remove from start(0) or end(1)
%
% from https://mrrl.ucla.edu/file/59559/M229_Lecture14_PartialFourier.pdf

if ~exist("refim","var");refim=[];end

Fz=pinv(Ez,0.5);
nz=size(data, 2);

if dir == 0
    idx1 = 1+nz-n-zf:nz;
    idx2 = 1:n+zf;
    idx3 = 1:nz-n-zf;
else
    idx1 = 1:n+zf;
    idx2 = 1+nz-n-zf:nz;
    idx3 = 1+n+zf:nz;
end

data_pk=data;
data_pk(:,idx1,:,:) = 0;

data_center = data_pk;
data_center(:,idx2,:,:) = 0;

im_ph = quickRecon(data_center);

im_init = quickRecon(data_pk);
im_init = abs(im_init).*exp(1i*angle(im_ph));

tmp_k =  quickEncode(im_init);

MSEs=zeros(1,its);
for ii=1:its
    tmp_k(:,idx3,:,:) = data_pk(:,idx3,:,:);
    if zf > 1
        tmp_k(:,[1:zf end-zf+1:end],:,:) = 0;
    end
    tmp_im = quickRecon(tmp_k);
    
    tmp_im = abs(tmp_im).*exp(1i*angle(im_ph));
    tmp_k = quickEncode(tmp_im);
    
    im=abs(tmp_im);
    if ~isempty(refim)
        MSEs(:,ii)=mse(im(:)./max(im(:)),refim(:)./max(refim(:)));
    end

    if zf > 1
        tmp_k(:,[1:zf end-zf+1:end],:,:) = 0;
    end
end

im = tmp_im;
ksp = tmp_k;

function im=quickRecon(data)
    im = reshape((Fz*(Fxy*data).').',mtx);
end

function data=quickEncode(im)
    data = (Ez*(Exy*reshape(im,[],mtx(3))).').';
end

end

