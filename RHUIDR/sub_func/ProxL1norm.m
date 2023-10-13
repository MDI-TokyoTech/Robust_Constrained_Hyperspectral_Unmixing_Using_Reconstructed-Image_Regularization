function[X] = ProxL1norm(X,gamma)

X = sign(X) .* max( abs(X) - gamma, 0 ); 