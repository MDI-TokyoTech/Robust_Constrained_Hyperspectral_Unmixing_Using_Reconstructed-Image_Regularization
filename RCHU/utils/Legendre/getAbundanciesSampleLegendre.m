function [abundanciesSample] = getAbundanciesSampleLegendre(ne,npoints,minOrder,maxOrder,maxCoef)
% [abundanciesSample] = getAbundanciesSampleLegendre(ne,npoints,minOrder,maxOrder,maxCoef)
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
% Return an abundancies sample by Legendre polynomies method.
%
% Input:
%   - ne: number of images in the sample (default: 5).
%   - npoints: spatial dimensions of the images, the images have the same
%   width and height (default: 256).
%   - minOrder,maxOrder: minimum and maximum orders of the Legendre
%   polynomies (defaults: 1, 10).
%   - maxCoef: to weight the Legendre polynomies (default: 100).
%
% Output:
%   - abundanciesSample: a 3-dimensional matrix (width x height x ne).

%% Check parameters
switch nargin
    case 0
        ne = 5; npoints = 256; minOrder = 1; maxOrder = 10; maxCoef = 100;
    case 1
        npoints = 256; minOrder = 1; maxOrder = 10; maxCoef = 100;
    case 2
        minOrder = 1; maxOrder = 10; maxCoef = 100;
    case 3
        maxOrder = 10; maxCoef = 100;
    case 4
        maxCoef = 100;
end
if ne <= 0
    warning('The number of images in sample has to be greater than zero. Default used (ne = 5).\n');
    ne = 5;
end
if npoints <= 0
    warning('The number of points in the images has to be greater than zero. Default used (npoints = 256).\n');
    npoints = 256;
end
if minOrder <= 0
    warning('The minimum Legendre order has to be greater than zero. Default used (minOrder = 1).\n');
    minOrder = 1;
end
if maxOrder <= 0
    warning('The maximum Legendre order has to be greater than zero. Default used (maxOrder = 10).\n');
    maxOrder = 10;
end
if maxOrder <= minOrder
    warning('The maximum order has to be greater than minimum order. minOrder + 1 used.\n');
    maxOrder = minOrder + 1;
end
if maxCoef <= 0
    warning('The maximum coeficient has to be greater than zero. Default used (maxCoef = 100).\n');
    maxCoef = 100;
end


%% Legendre
tic;
abundanciesSample = zeros(npoints,npoints,ne);

for i = 1:ne
    Q = zeros(npoints);
    for k1 = minOrder:maxOrder
        for k2 = minOrder:maxOrder
            theta1 = (rand-0.5)*2*maxCoef;
            theta2 = (rand-0.5)*2*maxCoef;
            P = getLegendre(npoints,k1,k2,theta1,theta2);
            Q = Q + P;
        end
    end
    abundanciesSample(:,:,i) = normalize(Q);
end
abundanciesSample = renormalize(abundanciesSample);

t = toc;
fprintf('Sample completed in %g seconds.\n',t);