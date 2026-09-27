%%-------------------------------------------------------------------------
%% load data
%%-------------------------------------------------------------------------


set(0,'DefaultFigureWindowStyle','docked')

addpath data/
addpath utils/


load phs_tissue;            % tissue phase from transversal orientation (in ppm, normalized by gyro*TE*B0)
load phs_wrap;              % raw phase from trans orient without processing (in radians)
load phs_unwrap;            % phase from trans orient processed with Laplacian unwrapping and BET masking (in radians)


load spatial_res;           % voxel size
load msk;                   % brain mask => obtained by eroding the BET mask by 5 voxels (by setting peel=5 in LBV)


load magn;                  % magnitude from transversal orientation
load magn_raw;              % raw magnitude without masking from transversal orientation
load mp_rage;               % mprage


load chi_33;                % chi_33 from STI solution (in ppm)
load chi_cosmos;            % COSMOS from 12 orientations (in ppm)     
 

N = size(msk);


imagesc3d2(chi_cosmos, N/2, 1, [90,90,-90], [-0.10,0.14], [], 'Cosmos')
imagesc3d2(chi_33, N/2, 2, [90,90,-90], [-0.10,0.14], [], '\chi_{33}')


imagesc3d2(msk, N/2, 3, [90,90,-90], [0,1], [], 'Mask')


imagesc3d2(phs_tissue, N/2, 4, [90,90,-90], [-0.05,0.05], [], 'Input Phase')
imagesc3d2(phs_wrap, N/2, 5, [90,90,-90], [-pi,pi], [], 'Raw Phase')
imagesc3d2(phs_unwrap, N/2, 6, [90,90,-90], [-pi,pi], [], 'Unwrapped Phase')


imagesc3d2(magn, N/2, 7, [90,90,-90], [0,0.5], [], 'Magnitude')
imagesc3d2(magn_raw, N/2, 8, [90,90,-90], [0,0.5], [], 'Raw Magnitude')
imagesc3d2(mp_rage, N/2, 9, [90,90,-90], [0,0.9], [], 'MPRAGE')



%%-------------------------------------------------------------------------
%% create dipole kernel
%%-------------------------------------------------------------------------


[ky,kx,kz] = meshgrid(-N(1)/2:N(1)/2-1, -N(2)/2:N(2)/2-1, -N(3)/2:N(3)/2-1);

kx = (kx / max(abs(kx(:)))) / spatial_res(1);
ky = (ky / max(abs(ky(:)))) / spatial_res(2);
kz = (kz / max(abs(kz(:)))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;


R_tot = eye(3);     % orientation matrix for transverse acquisition


kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ (k2 + eps) );    


kernel_disp = fftshift(mean(abs(kernel), 4));
mosaic(squeeze(kernel_disp(1+end/2,:,:)), 1, 1, 11, '', [0,2/3])

 

%%-------------------------------------------------------------------------
%% TKD recon
%%-------------------------------------------------------------------------


thre_tkd = 0.19;      % TKD threshold parameter


kernel_inv = zeros(N);
kernel_inv( abs(kernel) > thre_tkd ) = 1 ./ kernel(abs(kernel) > thre_tkd);


chi_tkd = real( ifftn( fftn(phs_tissue) .* kernel_inv ) ) .* msk; 


imagesc3d2(chi_tkd, N/2, 21, [90,90,-90], [-0.10,0.14], [], 'TKD')



%%-------------------------------------------------------------------------
%% closed-form L2 recon
%%-------------------------------------------------------------------------

[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;

DtD = abs(kernel).^2;


reg_param = 9e-2;   % gradient regularization parameter


chi_L2 = real( ifftn(conj(kernel) .* fftn(phs_tissue) ./ (DtD + reg_param * EtE))) .* msk;


imagesc3d2(chi_L2, N/2, 22, [90,90,-90], [-0.10,0.14], [], 'CF L2')



%%-------------------------------------------------------------------------
%% compute performance metrics for TKD and CFL2 recons
%%-------------------------------------------------------------------------
 

% RMSE (lower the better, given as percentage, RMSE=0% is the best possible):
rmse_tkd_cosmos = compute_rmse(chi_tkd, chi_cosmos);
rmse_tkd_chi33 = compute_rmse(chi_tkd, chi_33);


rmse_L2_cosmos = compute_rmse(chi_L2, chi_cosmos);
rmse_L2_chi33 = compute_rmse(chi_L2, chi_33);


disp(['TKD RMSE: ', num2str(rmse_tkd_cosmos), ' % relative to COSMOS and ', num2str(rmse_tkd_chi33), ' % relative to Chi_33'])
disp(['CFL2 RMSE: ', num2str(rmse_L2_cosmos), ' % relative to COSMOS and ', num2str(rmse_L2_chi33), ' % relative to Chi_33'])



% HFEN (lower the better, given as percentage, HFEN=0% is the best possible):
hfen_tkd_cosmos = compute_hfen(chi_tkd, chi_cosmos);
hfen_tkd_chi33 = compute_hfen(chi_tkd, chi_33);


hfen_L2_cosmos = compute_hfen(chi_L2, chi_cosmos);
hfen_L2_chi33 = compute_hfen(chi_L2, chi_33);
 

disp(['TKD HFEN: ', num2str(hfen_tkd_cosmos), ' % relative to COSMOS and ', num2str(hfen_tkd_chi33), ' % relative to Chi_33'])
disp(['CFL2 HFEN: ', num2str(hfen_L2_cosmos), ' % relative to COSMOS and ', num2str(hfen_L2_chi33), ' % relative to Chi_33'])



% SSIM (higher the better, between [0,1], SSIM=1 is the best possible):
ssim_tkd_cosmos = compute_ssim(chi_tkd, chi_cosmos);
ssim_tkd_chi33 = compute_ssim(chi_tkd, chi_33);


ssim_L2_cosmos = compute_ssim(chi_L2, chi_cosmos);
ssim_L2_chi33 = compute_ssim(chi_L2, chi_33);


disp(['TKD SSIM: ', num2str(ssim_tkd_cosmos), ' % relative to COSMOS and ', num2str(ssim_tkd_chi33), ' % relative to Chi_33'])
disp(['CFL2 SSIM: ', num2str(ssim_L2_cosmos), ' % relative to COSMOS and ', num2str(ssim_L2_chi33), ' % relative to Chi_33'])



