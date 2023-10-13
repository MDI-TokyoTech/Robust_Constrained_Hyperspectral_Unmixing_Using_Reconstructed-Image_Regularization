function [lambda, aux]  = generation(n, m, tipo, theta1, theta2)
% [lambda, aux]  = generation(n, m, tipo, theta1, theta2)
%
%    Copyright 2007,2008 Grupo Inteligencia Computacional, Universidad del País Vasco / Euskal Herriko Unibertsitatea (UPV/EHU)
%
%    Website: http://www.ehu.es/ccwintco
%
%    This file is part of the HYperspectral Data Retrieval and Analysis tools (HYDRA) program.
%
%    HYDRA is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    HYDRA is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with HYDRA.  If not, see <http://www.gnu.org/licenses/>.

%% Generation
cova = getCovariances(n, m, tipo, theta1, theta2);
nb = n;
tt = blockmatrix (nb, 1, m, 1);
for j=1:nb,
   tt(j) = cova((j-1)*m+1:(j-1)*m+m);
end
vc = zeros(2*n, 1);
cc = blockmatrix (nb, 1, 2*m, 1);
for j = 1:nb,
   	aux = tt(j);
   	vc(1:m) = aux';
   	vc(m+1) = 0;
   	vc(m+2:2*m) = aux(m:-1:2);
  	cc(j) = vc';
end
vc = blockmatrix (2*n, 1, 2*m, 1);
for j = 1:n,
   vc(j) = cc(j);
end
vc(n+1) = zeros(2*m, 1);
for j = n:-1:2,
   vc(n+(n-j)+2) = cc(j);
end
aux = vc(:);
aux1 = reshape(aux, 2*m, 2*n);
lambda = real(fft2(aux1));
lambda = reshape (lambda, 1, 4*n*m);
