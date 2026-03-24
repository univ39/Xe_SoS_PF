%% Numerical phantom simulations for stack-of-spirals with Partial Fourier
% This script runs the simulations, which can take a long time depending on
% compute resources. The simulated data from the paper can also be directly
% plotted using figure_sim_results.m.
close all;
addpath('dependencies/')
clearvars -except F_xy E_xy
%% load in waveform
wfn='waveforms/spiral4D_out_129xe_fov400_250_mtx80_51_arms16_16_kdt4_gmax32_smax119_dur3p2_numecho1_GA0_zPE4.mat';TR=0.0156;
% wfn='waveforms/spiral4D_out_129xe_fov400_250_mtx80_51_arms12_12_kdt4_gmax33_smax129_dur4_numecho1_GA0_zPE4.mat'; TR=0.0078;

load(wfn);
[mtx_acq,mtx_reco,npts_arm,narms,npts,k,nexc] = getKeyWFVars(wf);
nz_homodyne=15; n_fullctr=mtx_reco(3)-nz_homodyne.*2; nexc_pf=(mtx_reco(3)-nz_homodyne).*narms;

%% get encoding and reconstruction matrices
if ~exist('E_xy','var'); E_xy = getE(k(:,1:npts,1:2),mtx_acq,mtx_reco,2); end
if ~exist('F_xy','var'); F_xy = getF(E_xy,[],'cholesky');end

%% design different Z encoding schemes
wf_io=wf;
wf_oi=wf;wf_oi.k(:,:,3)=flip(wf.k(:,:,3));

% get Kz for PF
z = wf.k(:,:,3);
tmp=ones(1,51);tmp(n_fullctr+2:2:end)=0;
tmp2=repelem(tmp,1,npts);
wf_pf=wf;wf_pf.k=wf.k(:,1:npts_arm*nexc_pf,:);

wf_io_pf=wf_pf;wf_io_pf.k(:,:,3)=z(:,tmp2==1);
PF_io.dir='io';PF_io.n_fullctr=n_fullctr;PF_io.nz_homodyne=nz_homodyne;
PF_io.full_k=wf_io.k;

wf_oi_pf=wf_pf;wf_oi_pf.k(:,:,3)=flip(z(:,tmp2==1));
PF_oi.dir='oi';PF_oi.n_fullctr=n_fullctr;PF_oi.nz_homodyne=nz_homodyne;
PF_oi.full_k=wf_oi.k;

figure;plot(CFA_pf);hold on;plot(CFA);plot(RFA_pf);plot(RFA);hold off
legend('CFA (undersampled)','CFA','RFA(undersampled)','RFA','location','southeast')
set(findall(gca, 'Type', 'Line'), 'LineWidth', 2)
xlabel('N_{exc}')
ylabel('FA [{^o}]')
params.T1=23;

%% Parameters to iterate over
Noise_range=10.^[-3:0.2:2];
B1_range=0.25:0.05:1.75; % global

%% initialize
iterations=25;
io=initialize_structs(iterations,Noise_range,B1_range);
io_pf=initialize_structs(iterations,Noise_range,B1_range);

io_ramped=initialize_structs(iterations,Noise_range,B1_range);
oi_ramped=initialize_structs(iterations,Noise_range,B1_range);
io_ramped_pf=initialize_structs(iterations,Noise_range,B1_range);
oi_ramped_pf=initialize_structs(iterations,Noise_range,B1_range);

%% simulate
idx=1;
zPEmode='subsequent';
for ni=1%:iterations
    for nnoise=1%:length(Noise_range)
        for nb1=1%:length(B1_range)
            progressBar(idx,iterations*length(Noise_range)*length(B1_range),'Simulating...')
            idx=idx+1;

            params.noise=Noise_range(nnoise);
            params.B1=B1_range(nb1);
            params.F_xy=F_xy;
            params.E_xy=E_xy;
            params.TR=TR;
            params.depol_signal=true;

            % % constant flip angle
            [~,MSE,PSNR,SSIM] = sl3d_sim_wrapper(wf_io,CFA,0,zPEmode,params);
            io=store_in_struct(io,MSE,PSNR,SSIM,nnoise,nb1,ni);
            % constant flip angle and PF
            [~,MSE,PSNR,SSIM] = sl3d_sim_wrapper(wf_io_pf,CFA_pf,0,zPEmode,params,PF_io);
            io_pf=store_in_struct(io_pf,MSE,PSNR,SSIM,nnoise,nb1,ni);


            % % ramped flip angle
            [~,MSE,PSNR,SSIM] = sl3d_sim_wrapper(wf_io,RFA,0,zPEmode,params);
            io_ramped=store_in_struct(io_ramped,MSE,PSNR,SSIM,nnoise,nb1,ni);
            [~,MSE,PSNR,SSIM] = sl3d_sim_wrapper(wf_io,RFA,0,zPEmode,params);
            oi_ramped=store_in_struct(oi_ramped,MSE,PSNR,SSIM,nnoise,nb1,ni);

            % ramped flip angle and PF
            [~,MSE,PSNR,SSIM] = sl3d_sim_wrapper(wf_io_pf,RFA_pf,0,zPEmode,params,PF_io);
            io_ramped_pf=store_in_struct(io_ramped_pf,MSE,PSNR,SSIM,nnoise,nb1,ni);
            [~,MSE,PSNR,SSIM] = sl3d_sim_wrapper(wf_oi_pf,RFA_pf,0,zPEmode,params,PF_oi);
            oi_ramped_pf=store_in_struct(oi_ramped_pf,MSE,PSNR,SSIM,nnoise,nb1,ni);
        end
    end
end

clearvars E_xy F_xy
save(sprintf('simresults_%s.mat',string(datetime('today','Format','yyyyMMdd'))));



%% subfunction
function out=initialize_structs(iterations,Noise_range,B1_range)
out=struct();
out.MSE=zeros(length(Noise_range),length(B1_range),iterations);
out.PSNR=zeros(length(Noise_range),length(B1_range),iterations);
out.SSIM=zeros(length(Noise_range),length(B1_range),iterations);
end

function out=store_in_struct(out,MSE,PSNR,SSIM,nsnr,nb1,ni)
out.MSE(nsnr,nb1,ni)=MSE;
out.PSNR(nsnr,nb1,ni)=PSNR;
out.SSIM(nsnr,nb1,ni)=SSIM;
end