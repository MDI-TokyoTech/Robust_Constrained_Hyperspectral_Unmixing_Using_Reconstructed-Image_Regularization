%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the demo file of the method proposed in the
% following reference:
% 
% K. Naganuma and S. Ono
% ``Towards Robust Hyperspectral Unmixing: Mixed Noise Modeling and Image-Domain Regularization''
%
% Update history:
% October 13, 2023: v1.0 
%
% Copyright (c) 2023 Kazuki Naganuma and Shunsuke Ono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Adding path
addpath('./methods')
addpath('./sub_func')

%% Setting parameters
% size
n1 = 64; % the number of vertical pixels
n2 = 64; % the number of horizontal pixels
m = 10; % the number of spectral signatures included in an endmember library
l = 224; % the number of bands
num_endmembers_exist = 4;

% data type
data_flag = 1; % 1 for Gaussian, 2 for Legendre

% use GPU or not
use_GPU = 1; % 0 if you do not use GPU, 1 if you use GPU

% noise
sigma = 0.1;
Sp_rate = 0.05;

% image-domain regularization
type_IDR = 'HTV';
% type_IDR = 'SSTV';
% type_IDR = 'HSSTV';

%% Generating datasets
params.n1 = n1;
params.n2 = n2;
params.m = m;
params.num_endmembers_exist = num_endmembers_exist;
params.data_flag = data_flag;
params.use_GPU = use_GPU;
params.sigma = sigma;
params.Sp_rate = Sp_rate;

data = make_synthdata(params);

HSI_noisy = data.HSI_noisy;
HSI_org = data.HSI_org;
X_true = data.X_true;
X_true_exp = data.X_true_exp;
abundance_exp = data.abundance_exp;
abundance = data.abundance;
idxs_exist = data.idxs_exist;
A = data.A;
A_exist = data.A_exist;
names_endmembers = data.names_endmembers;

%% Preparation of unmixing
% hyperparameters in the objective function
lambda1 = 2;
lambda2 = 2;
lambda3 = 1;

% Observed HS image
Y = cube2mat(HSI_noisy, n1, n2, l);

% parameters for noise
para_ep = 0.95;
para_et = 0.9;
ep = para_ep*sigma*sqrt(n1*n2*l*(1 - Sp_rate));
eta = para_et*(Sp_rate*0.5*n1*n2*l);

%% Unmixing
params.ep = ep;
params.eta = eta;
params.lambda1 = lambda1;
params.lambda2 = lambda2;
params.lambda3 = lambda3;

sizes.n1 = n1;
sizes.n2 = n2;
sizes.m  = m;
sizes.l  = l;

if strcmp(type_IDR, 'HTV')
    results = RHUIDR_HTV(A,Y,params,sizes);
elseif strcmp(type_IDR, 'SSTV')
    results = RHUIDR_SSTV(A,Y,params,sizes);
elseif strcmp(type_IDR, 'HSSTV')
    results = RHUIDR_SSTV(A,Y,params,sizes);
else
    disp('The image-domain regularization cannot be selected')
end

X_est = results.X;

%% Calculating metrics
U_est = A*X_est; 
HSI_rec = reshape(U_est', n1, n2, l);
abundance_est_exist = reshape(X_est(idxs_exist,:)',n1,n2,num_endmembers_exist);

val_SRE = SRE(abundance_est_exist, abundance);
val_RMSE = RMSE(abundance_est_exist, abundance);
val_Ps = Ps(abundance_est_exist, abundance, 0.316);
val_PSNR = PSNR_array(HSI_rec, HSI_org, 1);
val_SSIM = ssim_index3d(255*HSI_rec, 255*HSI_org, [1 1 1], HSI_rec>0);

%% Displaying metrics
disp('****************** Noisy HSI *************************')
val_PSNR_noisy = PSNR_array(HSI_noisy, HSI_org, 1);
disp(append('PSNR of noisy image : ', num2str(val_PSNR_noisy)));
disp('****************** Unmixing results ******************')
disp(append('SRE  : ', num2str(val_SRE)));
disp(append('RMSE : ', num2str(val_RMSE)));
disp(append('Ps   : ', num2str(val_Ps)));
disp(append('PSNR : ', num2str(val_PSNR)));
disp(append('SSIM : ', num2str(val_SSIM)));

%% Plotting results
num_band_show = 100;
fig = figure('Name', append('RHUIDR-', type_IDR));
fig.Position(2) = 100;
fig.Position(3) = 400*num_endmembers_exist;

subplot(2, num_endmembers_exist + 1, 1)
imshow(HSI_noisy(:, :, num_band_show), [0 1]);
title('Noisy HS image ');

subplot(2, num_endmembers_exist + 1, num_endmembers_exist + 2);
imshow(HSI_rec(:, :, num_band_show), [0 1]);
title('Reconstructed HS image');

for idx_abun = 1:num_endmembers_exist
    subplot(2, num_endmembers_exist + 1, idx_abun + 1)
    imshow(abundance(:, :, idx_abun), [0 1], Colormap=parula);
    title(append('True abundnance ', num2str(idx_abun)));

    subplot(2, num_endmembers_exist + 1, num_endmembers_exist + idx_abun + 2);
    imshow(abundance_est_exist(:, :, idx_abun), [0 1], Colormap=parula);
    title(append('Estimated abundnance ', num2str(idx_abun)));
end

%% functions
function result=PSNR_array(recon,org,skip)
 [m, n,~]=size(org);
    vec_PSNRs = 10*log10(m*n./sum((recon - org).^2, [1, 2]));
    result = mean(vec_PSNRs, "all");

end

function [result] = SRE(recon, org)
    result = 10*log10(sum(org.^2, "all")/sum((recon - org).^2, "all"));
end


function [result] = RMSE(recon, org)
    [m,n,k] = size(org);
    result = sqrt(sum((recon - org).^2, "all")/(n*m*k));
end

function [Ps] = Ps(recon, org, thres)
    [m,n,~] = size(org);
    tmp = sum((recon - org).^2, 3)./sum(org.^2, 3);
    map_prob = double(tmp <= thres);
    Ps = sum(map_prob, "all")/(m*n);
end
%% code ends here.
