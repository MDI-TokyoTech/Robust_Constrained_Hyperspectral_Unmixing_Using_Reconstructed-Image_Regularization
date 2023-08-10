%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the file of the method proposed in the
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

function results = RCHU_SSTV(A, Y, params, sizes)

ep = params.ep;
eta = params.eta;
lambda1 = params.lambda1;
lambda2 = params.lambda2;
lambda3 = params.lambda3;

n1 = sizes.n1;
n2 = sizes.n2;
m  = sizes.m;
l  = sizes.l;

if ~isfield(params, 'tol')
    tol = 1e-5;
else 
    tol = params.tol;
end
if ~isfield(params, 'num_iter_max')
    num_iter_max = 100000;
else 
    num_iter_max = params.num_iter_max;
end
if ~isfield(params, 'use_GPU')
    use_GPU = 100000;
else 
    use_GPU = params.use_GPU;
end

%% Setting functions
Db = @(z) z(:,:,[2:end, 1])-z;
Dbt = @(z) z(:,:,[end,1:end-1]) - z;
Dv = @(z) z([2:end, 1],:,:) - z;
Dvt = @(z) z([end,1:end-1],:,:) - z;
Dh = @(z) z(:,[2:end, 1],:)-z;
Dht = @(z) z(:,[end,1:end-1],:) - z;
Dhsi = @(z) cat(4, Dv(Db(z)), Dh(Db(z)));
Dhsit = @(z) Dbt(Dvt(z(:,:,:,1))) + Dbt(Dht(z(:,:,:,2)));
Mat2HSI_hsi = @(z) reshape(z',n1,n2,l);
HSI2Mat_hsi = @(z) reshape(z,n1*n2,l)';
K = @(z) Dhsi(Mat2HSI_hsi(z));
Kt = @(z) HSI2Mat_hsi(Dhsit(z));

D3 = @(z) Dv(Mat2HSI_hsi(z));
D3t = @(z) HSI2Mat_hsi(Dvt(z));

Dabun = @(z) cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:)-z);
Dabunt = @(z) [-z(1,:,:,1); - z(2:end-1,:,:,1) + z(1:end-2,:,:,1); z(end-1,:,:,1)] ...
            +[-z(:,1,:,2), - z(:,2:end-1,:,2) + z(:,1:end-2,:,2), z(:,end-1,:,2)];
Mat2HSI_abun = @(z) reshape(z',n1,n2,m);
HSI2Mat_abun = @(z) reshape(z,n1*n2,m)';

D = @(z) Dabun(Mat2HSI_abun(z));
Dt = @(z) HSI2Mat_abun(Dabunt(z));

%% Initializings variables
if use_GPU == 1
    X = gpuArray(zeros(size(A'*Y)));
    S = gpuArray(zeros(size(Y)));
    L = gpuArray(zeros(size(Y)));
    Z1 = gpuArray(X);
    Z2 = gpuArray(D(X));
    Z3 = gpuArray(K(A*X));
    Z4 = gpuArray(Y);
    Z5 = gpuArray(D3(A*X));
else
    X = zeros(size(A'*Y));
    S = zeros(size(Y));
    L = zeros(size(Y));
    Z1 = X;
    Z2 = D(X);
    Z3 = K(A*X);
    Z4 = Y;
    Z5 = D3(A*X);
end

%% Calculating stepsizes
sig_A_max = max(svd(gather(A), 'econ'));
gamma_1 = 1/(9 + 33*sig_A_max^2);
gamma_2 = 1;
gamma_3 = 1/5;
gamma_4 = 1/3;

%% Start algorithm

disp('****************** Algorithm starts ******************')

for i = 1:num_iter_max
    
    Xpre = X;
    Spre = S;
    Lpre = L;

    % Updating X : steps 4, 5, and 6
    X = X - gamma_1 * (Z1 + Dt(Z2) + A'*(Kt(Z3)) + A'*Z4);
    X = max(0,X); % sum-to-oneを除いて非負制約のみに

    % Updating S : steps 7 and 8
    S = S - gamma_2 * Z4;
    S = ProjFastL1Ball(S,eta); 
    
    % Updating L : steps 9 and 10
    L = L - gamma_3 * (Z4 + D3t(Z5));
    L = ProxL1norm(L,gamma_3);
    
    % Updating Z1 : steps 11 and 12
    Z1 = Z1 + gamma_4 * (2*X - Xpre);
    Z1temp = ProxRowGroupL12norm(Z1./gamma_4,lambda1./gamma_4);
    Z1 = Z1 - gamma_4 * Z1temp;    
    
    % Updating Z2 : steps 13 and 14
    Z2 = Z2 + gamma_4 * D(2*X - Xpre);
    Z2temp = ProxL1norm(Z2/gamma_4,lambda2/gamma_4); 
    Z2 = Z2 - gamma_4 * Z2temp;

    % Updating Z3 : steps 15 and 16
    Z3 = Z3 + gamma_4 * K(A*(2*X - Xpre));
    Z3temp = ProxL1norm(Z3/gamma_4,lambda3/gamma_4); 
    Z3 = Z3 - gamma_4 * Z3temp;
    
    % Updating Z4 : steps 17, 18, 19, and 20
    Z4 = Z4 + gamma_4 * (A*(2*X - Xpre) + (2*S - Spre) + (2*L - Lpre));
    Z4temp = ProjL2Ball(Z4./gamma_4,Y,ep); 
    Z4 = Z4 - gamma_4 * Z4temp;

    % Updating Z5 : steps 21 and 22
    Z5 = Z5 + gamma_4 * D3(2*L - Lpre);
    Z5temp = 0;
    Z5 = Z5 - gamma_4 * Z5temp;
    
    % stopping condition
    Res = X-Xpre;
    error = norm(Res(:),2)/norm(X(:),2);
      
    if error < tol
        break;
    end
   
end

disp('****************** Algorithm ends ********************')

results.X = gather(X);
results.S = gather(S);
results.L = gather(L);

