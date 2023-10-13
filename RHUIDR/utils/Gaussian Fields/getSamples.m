function [z1, z2] = getSamples(lambda, n, m, nh)
% [z1, z2] = getSamples(lambda, n, m, nh)
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

%% Samples
ndata = sqrt(size(lambda, 2)/4);
norma = sqrt(4*ndata*ndata);

kk = [];
for j =0:ndata-1,
   kk = [kk; 2*ndata*j:1:2*ndata*j+(ndata-1)];
end
kk = kk + 1;
ind_k = reshape(kk, ndata, ndata);

z1 = zeros (n, m, nh); z2 = zeros (n, m, nh);
a1 = randn(nh, 4*ndata*ndata);
a2 = randn(nh, 4*ndata*ndata);
for i=1:nh,
	w1 = a1(i, :).*sqrt(lambda);
	w2 = a2(i, :).*sqrt(lambda);
   	w = fft2(w1+i*w2);
   %%%% Not sure about normalization
  	w = w / norma;
   	b1 = real(w(ind_k)); 
	b2 = imag(w(ind_k));
	z1 (:, :, i) = b1(1:n, 1:m);
	z2 (:, :, i) = b2(1:n, 1:m);
end
