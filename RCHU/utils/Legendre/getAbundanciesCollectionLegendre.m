function [abundanciesCollection] = getAbundanciesCollectionLegendre(n,ne,npoints,minOrder,maxOrder,maxCoef)
% [abundanciesCollection] = getAbundanciesCollectionLegendre(n,ne,npoints,minOrder,maxOrder,maxCoef)
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
%
% Returns a collection of synthetic samples of abundancy images by Legendre
% polynomies method.
%
% Input:
%   - n: number of samples to build (> 0, defaults are 10).
%   - ne: number of images in each sample (>= 2, default is 5).
%   - npoints: spatial dimensions of the images (> 0, default is 256).
%   - minOrder,maxOrder: minimum and maximum orders of the Legendre
%   polynomies (defaults: 1, 10).
%   - maxCoef: to weight the Legendre polynomies (default: 100).
%
% Output:
%   - abundanciesCollection: a 4-dimensional matrix (n x width x height x ne).

%% Check parameters
% deafults
switch nargin
    case 0
        n = 10; ne = 5; npoints = 256; minOrder = 1; maxOrder = 10; maxCoef = 100;
    case 1
        ne = 5; npoints = 256; minOrder = 1; maxOrder = 10; maxCoef = 100;
    case 2
        npoints = 256; minOrder = 1; maxOrder = 10; maxCoef = 100;
    case 3
        minOrder = 1; maxOrder = 10; maxCoef = 100;
    case 4
        maxOrder = 10; maxCoef = 100;
    case 5
        maxCoef = 100;
end
% Check n
if n <= 0
	error('n has to be greater than zero');
end
% More checks made in getAbundanciesSampleLegendre;

%% Build samples
% Initialize collection
abundanciesCollection = zeros(n,npoints,npoints,ne);
% Build
for i=1:n
    fprintf('Synthesizing sample: %g.\n',i);
    abundanciesCollection(i,:,:,:) = getAbundanciesSampleLegendre(ne, npoints, minOrder, maxOrder, maxCoef);
end
fprintf('Collection completed.\n');
