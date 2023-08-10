function data = make_synthdata(params)

%% Adding path
addpath(genpath('./utils'));

%% Setting parameters
if ~isfield(params, 'num_endmembers_exist')
    num_endmembers_exist = 4;
else 
    num_endmembers_exist = params.num_endmembers_exist;
end
if ~isfield(params, 'm')
    m = 10; % the number of spectral signatures included in the endmember library
else
    m = params.m;
end
if ~isfield(params, 'n1')
    n1 = 64; % the number of vertical pixels
else
    n1 = params.n1;
end
if ~isfield(params, 'n2')
    n2 = 64; % the number of horizontal pixels
else
    n2 = params.n2;
end
if ~isfield(params, 'data_flag')
    data_flag = 1;
else
    data_flag = params.data_flag;
end
if ~isfield(params, 'use_GPU')
    use_GPU = 0;
else
    use_GPU = params.use_GPU;
end

%% Loading an spectral library
l = 224; % the number of bands
load('./Endmembers/Spectra_224.mat', 'Libs_endmembers');
num_Libs = numel(Libs_endmembers);


%% Preparing the endmember library for experiments
A = zeros(l, m); % endmember library
idxs_tmp = randperm(num_Libs, m);
names_endmembers = cell(1, m);

for idx_tmp = 1:m
    idx_Lib = idxs_tmp(idx_tmp);
    names_endmembers{idx_tmp} = Libs_endmembers{idx_Lib}.name;
    A(:, idx_tmp) = transpose(Libs_endmembers{idx_Lib}.reflectance);
end

if use_GPU == 1
    A = gpuArray(A);
end

% idx = idxs_exist
idxs_exist = randperm(m, num_endmembers_exist);
idxs_exist = sort(idxs_exist, 'ascend');
A_exist = A(:, idxs_exist);

clear idx_tmp idxs_tmp idx_Lib;

%% Generating abundance

if data_flag == 1
    abundance = getAbundanciesSampleGaussianFields(num_endmembers_exist,n1,n2);
elseif data_flag == 2
    abundance = getAbundanciesSampleLegendre(num_endmembers_exist,n1,1,1);
else
    disp(append('The data type ', num2str(data_flag), 'is wrong'))
end

abundance_exp = zeros(n1, n2, m);
abundance_exp(:, :, idxs_exist) = abundance;

if use_GPU == 1
    abundance = gpuArray(abundance);
    abundance_exp = gpuArray(abundance_exp);
end

X_true = cube2mat(abundance, n1, n2, num_endmembers_exist);
X_true_exp = cube2mat(abundance_exp, n1, n2, m);


%% setting noise condition
if ~isfield(params, 'sigma')
    sigma = 0.10;
else
    sigma = params.sigma;
end
if ~isfield(params, 'Sp_rate')
    Sp_rate = 0.05;
else
    Sp_rate = params.Sp_rate;
end

HSI_org = mixing(A_exist, abundance);

HSI_noisy = addGauNoise(HSI_org, sigma);
HSI_noisy = addStNoise(HSI_noisy);
HSI_noisy = addSpNoise(HSI_noisy, Sp_rate);

if use_GPU == 1
    HSI_noisy = gpuArray(HSI_noisy);
end

%% Returing variables
data.HSI_noisy = HSI_noisy;
data.HSI_org = HSI_org;
data.X_true = X_true;
data.X_true_exp = X_true_exp;
data.abundance_exp = abundance_exp;
data.abundance = abundance;
data.idxs_exist = idxs_exist;
data.A = A;
data.A_exist = A_exist;
data.names_endmembers = names_endmembers;

end

%% functions
function[Y] = addGauNoise(X, sigma)
     Y = X + sigma*randn(size(X));
end

function [X] = addSpNoise(X, Sp_rate)
     Sp = 0.5*ones(size(X));
     Sp = imnoise(Sp, 'salt & pepper', Sp_rate);
     X(Sp==0) = 0;
     X(Sp==1) = 1;
end

function HSI_striped = addStNoise(HSI_input)
intensity = 0.2;
rate_of_stripe = 0.05;
sigma_of_stripe = 0.05;
[n1, n2, n3] = size(HSI_input);
stripe1 = 2*(imnoise(0.5*ones(1, n2, n3), "salt & pepper", rate_of_stripe) - 0.5).*rand(1, n2, n3)*intensity.*ones(n1, n2, n3);
stripe2 = sigma_of_stripe*randn(1, n2, n3).*ones(n1, n2, n3);

stripe = stripe1 + stripe2;
stripe = intensity*stripe/(max(stripe, [], "all"));

HSI_striped = HSI_input + stripe;

end