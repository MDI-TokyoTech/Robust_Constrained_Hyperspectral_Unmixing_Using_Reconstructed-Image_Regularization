function[X] = ProxRowGroupL12norm(X,gamma)

thresh = (sqrt(sum(X.^2,2)).^(-1))*gamma; % グループ（行）ごとのしきい値 γ/||x||_2の計算
thresh(thresh>1) = 1; % しきい値が１を超えるところを１に置き換え
onemat = ones(size(thresh)); % thresh と同じサイズのオール１ベクトル
coef = repmat(onemat - thresh,1,size(X,2)); % max{1 - γ/||x||_2,0}部に相当

X = coef.*X; % 計算結果