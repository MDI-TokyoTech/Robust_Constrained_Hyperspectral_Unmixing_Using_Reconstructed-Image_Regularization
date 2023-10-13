function[X] = ProxTVL12norm(X,gamma)

thresh = (sqrt(sum(sum(X.^2,4),3)).^(-1))*gamma; % グループ（差分の方向数＋バンド数）ごとのしきい値 γ/||x||_2の計算
thresh(thresh>1) = 1; % しきい値が１を超えるところを１に置き換え
onemat = ones(size(thresh)); % thresh と同じサイズのオール１ベクトル
coef = repmat(onemat - thresh,1,1,size(X,3),size(X,4)); % max{1 - γ/||x||_2,0}部に相当

X = coef.*X; % 計算結果