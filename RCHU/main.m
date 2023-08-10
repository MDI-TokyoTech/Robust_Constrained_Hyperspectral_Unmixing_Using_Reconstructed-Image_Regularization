%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the demo file of the method proposed in the
% following reference:
% 
% K. Naganuma, Y. Nagamatsu, and S. Ono
% ``Robust Constrained Hyperspectral Unmixing Using Reconstructed-Image Regularization.''
%
% Update history:
% Augast 7, 2023: v1.0 
%
% Copyright (c) 2023 Kazuki Naganuma, Yuki Nagamatsu, and Shunsuke Ono
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

% reconstructed-image regularization
type_RIR = 'HTV';
% type_RIR = 'SSTV';
% type_RIR = 'HSSTV';

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

if strcmp(type_RIR, 'HTV')
    results = RCHU_HTV(A,Y,params,sizes);
elseif strcmp(type_RIR, 'SSTV')
    results = RCHU_SSTV(A,Y,params,sizes);
elseif strcmp(type_RIR, 'HSSTV')
    results = RCHU_SSTV(A,Y,params,sizes);
else
    disp('The reconstructed-image regularization cannot be selected')
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
fig = figure('Name', append('RCHU-', type_RIR));
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
function [cpsnr,psnr]=PSNR_array(recon,org,skip)
org=org(skip+1:end-skip,skip+1:end-skip,:);
recon=(recon(skip+1:end-skip,skip+1:end-skip,:));

  [m, n,~]=size(org);

    sse=squeeze(sum(sum((org-recon).^2))); %square sum of error  
    mse=sse./(m*n);  %mean square error of each band.
    maxval=squeeze(max(max(org)));
    psnr= 10*log10( (maxval.^2) ./mse);
    cpsnr=mean(psnr);

end

function [result] = SRE(recon, org)
    result = 10*log10(sum(org.^2, "all")/sum((recon - org).^2, "all"));
end


function [result] = RMSE(recon, org)
    [m,n,k] = size(org);
    result = sqrt(sum((recon - org).^2, "all")/(n*m*k));
end

function [Ps] = Ps(recon, org, thres)
    [m,n,k] = size(org);
    re_org = (reshape(org,m*n,k))';
    re_recon = (reshape(recon,m*n,k))';
    j = 0;
    for i=1:k
        for l = 1:m*n
        X_recon = abs(re_recon(i,l))*norm(re_recon(i,l));
            X_diff = abs(re_org(i,l)-re_recon(i,l))*abs(re_org(i,l)-re_recon(i,l));
            if (X_diff/X_recon) <= thres
                j = j+1;
            end
        end
    end
    Ps = j/(m*n*k);
end
%% code ends here.