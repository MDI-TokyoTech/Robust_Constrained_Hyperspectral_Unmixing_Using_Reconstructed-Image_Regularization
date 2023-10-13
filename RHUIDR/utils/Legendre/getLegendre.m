function [P] = getLegendre(npoints,k1,k2,theta1,theta2)
% [P] = getLegendre(npoints,k1,k2,theta1,theta2)
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
%    You should have received a copy of the GNU General Public License
%    along with HYDRA.  If not, see <http://www.gnu.org/licenses/>.
%
% Computes the Legendre polynomies method. The resulting matrix is computed by:
% theta1*legendre(k1,rank)'*theta2*legendre(k2,rank) where rank is
% an acumulative distribution from (0,1) along a vector of npoints.
%
% Input:
%   - npoints: spatial dimensions of the images, the images have the same
%   width and height.
%   - k1,k2: degrees of the Lengendre polinomies.
%   - theta1,theta2: weights for each Legendre polynomy
%
% Output:
%   - P: a 3-dimensional matrix (width x height x ne).

%% Legendre
rank = ((1:npoints)-1)/(npoints-1);

p1 = legendre(k1,rank);
p1 = (p1(1,:));

p2 = legendre(k2,rank);
p2 = (p2(1,:));

P = ((theta1*p1')*(theta2*p2));
