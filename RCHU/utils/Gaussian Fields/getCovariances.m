function cova = getCovariances(n, m, type, theta1, theta2)
% cova = getCovariances(n, m, type, theta1, theta2)
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

%% Covariances
ndata2 = n*m;
[x y] = meshgrid(0:1:n-1, 0:1:m-1);
xx = reshape(x, ndata2, 1); yy = reshape(y, ndata2, 1);
cova = zeros (1, ndata2);
d = sqrt((xx-xx(1)).*(xx-xx(1))+(yy-yy(1)).*(yy-yy(1)));
sv = covarianceType(d, type, theta1, theta2);
cova(:) = sv';
