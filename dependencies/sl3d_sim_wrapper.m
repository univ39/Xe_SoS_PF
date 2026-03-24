%% Forward modelling of a pseudo3d trajectory
% Kylie Yeung 05/2025

function [Image,MSE,PSNR,SSIM] = sl3d_sim_wrapper(wf,FA,depol_corr,zPEmode,params,PF)
%% get variables
[mtx_acq,mtx_reco,npts_arm,narms,npts,k,nexc] = getKeyWFVars(wf);
TR2T1=(params.TR)./(params.T1);
if ~exist('PF','var');PF=[];end

%% get encoding and reconstruction matrices
[F_z,k,dim,E_z] = getFz(k,mtx_acq, mtx_reco,npts.*mtx_reco(3));

%% Generate a discretized Shepp-Logan test image
[BILD_complex] = phantom3d_complex(mtx_reco);
BILD_reshape=reshape(BILD_complex,[mtx_reco(1)*mtx_reco(2) mtx_reco(3)]);

%% Forward encode the image into k-space data
d=params.E_xy*(E_z*BILD_reshape.').';

%% depolarize
if isempty(PF)
    switch zPEmode
        case 'subsequent'
            data_clean=reshape(d,[npts_arm narms.*size(d,2)]); % [#pts per exc, #exc]
            data=noise_depol(data_clean);
            data_reshaped=reshape(data,[npts mtx_acq(3)]).';
        case 'intlv'
            data_clean=reshape(permute(reshape(d,[npts_arm narms size(d,2)]),[1 3 2]),npts_arm,[]); % [#pts per exc, #exc]
            data=noise_depol(data_clean);
            data_reshaped=reshape(permute(reshape(data,[npts_arm mtx_acq(3) narms]),[2 1 3]),mtx_acq(3),[]);
        case 'hybrid'
            data_clean=reshape(d,[npts_arm narms.*size(d,2)]); % [#pts per exc, #exc]
            data_clean=cat(2,data_clean(:,1:end/2),flip(data_clean(:,end/2+1:end),2));
            data=noise_depol(data_clean);
            data=cat(2,data(:,1:end/2),flip(data(:,end/2+1:end),2));
            data_reshaped=reshape(data,[npts mtx_acq(3)]).';
        case 'variable_density'
    end
    Image=reshape(params.F_xy*(F_z*data_reshaped).',mtx_acq);

else
    %% reconstruction
    mtx_acq(3)=PF.n_fullctr+PF.nz_homodyne;
        switch zPEmode
        case 'subsequent'
            data_clean=reshape(d,[npts_arm narms.*size(d,2)]); % [#pts per exc, #exc]
            if params.depol_signal
                data=noise_depol(data_clean);
            else
                data_clean=data_clean.*sind(FA(1).*(params.B1)); % [#pts per exc, #exc]
                [data] = noise_data(data_clean,[],params.noise);
            end
            data_reshaped=reshape(data,[npts mtx_acq(3)]);
        end
        if strcmp(PF.dir,'io')    
            d_full = zeros(npts,mtx_reco(3));
            d_full(:,1:PF.n_fullctr)=data_reshaped(:,1:PF.n_fullctr);
            d_full(:,PF.n_fullctr+1:2:end)=data_reshaped(:,PF.n_fullctr+1:end);
            
            % with Phase correction
            data_sorted=cat(2,flip(d_full(:,1:2:end),2),d_full(:,2:2:end));
            z = PF.full_k(:,:,3);[~,i,~]=unique(z);
            Kz = single(z(sort(i)))./0.5*pi;
            Kz_sorted=cat(2,flip(Kz(:,1:2:end),2),Kz(:,2:2:end));
            ZM = single(-(mtx_reco(3)-1)/2:(mtx_reco(3)-1)/2);        
            E_z = exp(1i*(-Kz_sorted(:)*ZM(:).'));
            [Image] = PinvItPartialFourier_xyz(data_sorted, params.E_xy, params.F_xy, E_z, mtx_reco,10, 0, PF.nz_homodyne, 1, []);
    
        elseif strcmp(PF.dir,'oi')
            d_full = zeros(npts,mtx_reco(3));
            d_full(:,1:2:PF.nz_homodyne*2)=data_reshaped(:,1:PF.nz_homodyne);
            d_full(:,end-PF.n_fullctr:end)=data_reshaped(:,end-PF.n_fullctr:end);
            % with Phase correction
            data_sorted=cat(2,d_full(:,1:2:end),flip(d_full(:,2:2:end),2));
            z = PF.full_k(:,:,3);[~,i,~]=unique(z);
            Kz = single(z(sort(i)))./0.5*pi;
            Kz_sorted=cat(2,Kz(:,1:2:end),flip(Kz(:,2:2:end),2));
            ZM = single(-(mtx_reco(3)-1)/2:(mtx_reco(3)-1)/2);        
            E_z = exp(1i*(-Kz_sorted(:)*ZM(:).'));
            [Image] = PinvItPartialFourier_xyz(data_sorted, params.E_xy, params.F_xy, E_z, mtx_reco,10, 0, PF.nz_homodyne, 0, []);
        else
            error('PF direction unknown')
        end
        % [F_z,k,dim,E_z] = getFz(PF.full_k,mtx_acq, mtx_reco,npts.*mtx_reco(3));
        % Image = reshape((F_z*d_full.').',mtx_reco);
        % % figure;mat2montage(Image);
end

% figure;nexttile;mat2montage(BILD_complex);nexttile;mat2montage(Image);

%% get metrics
[MSE,PSNR,SSIM]=getImMetrics(BILD_complex,Image);

%% subfunctions
    function data=noise_depol(data_clean)
        % scale by sine
        data_clean=bsxfun(@times, data_clean, sind(FA.*(params.B1))); % [#pts per exc, #exc]
        data=depolarize(data_clean,FA.*(params.B1),TR2T1);

        % add noise
        [data] = noise_data(data,[],params.noise);
        
        % depolarize corr
        if depol_corr
            data=depolarize_corr(data,FA,TR2T1);
        end
    end

    function [MSE,PSNR,SSIM]=getImMetrics(ref_im,recon_im)    
        ref_im=ref_im./max(ref_im(:));
        recon_im=recon_im./max(recon_im(:));
    
        MSE=immse(single(gather(abs(recon_im))),single(abs(ref_im)));
        PSNR=psnr(single(gather(abs(recon_im))),single(abs(ref_im)));
        SSIM=ssim(single(gather(abs(recon_im))),single(abs(ref_im)));
    end

end